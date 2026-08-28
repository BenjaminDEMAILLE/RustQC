//! Per-contig depth of coverage accumulation.
//!
//! One [`DepthAccum`] covers one contig. Aligned blocks are recorded as
//! increments in a delta array the length of the contig, and a prefix sum at
//! the end turns that into per-base depth in a single linear pass.
//!
//! # Upstream semantics
//!
//! The filters and the CIGAR walk reproduce mosdepth 0.3.14 run without
//! `--fast-mode`, whose help text describes that flag as "dont look at
//! internal cigar operations or correct mate overlaps". Default mode
//! therefore does both, and so does this module:
//!
//! - records carrying any bit of [`MOSDEPTH_DEFAULT_EXCLUDE`] are skipped
//!   (mosdepth's `-F` default of 1796);
//! - records with `MAPQ` below the cutoff are skipped (mosdepth's `-Q`,
//!   default 0);
//! - `M`, `=` and `X` cover the reference, `D` and `N` advance without
//!   covering, and `I`, `S`, `H` and `P` do not advance at all.

use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use crate::common::bam_flags::*;

/// Bit mask matching mosdepth's `-F` default: `UNMAP | SECONDARY | QCFAIL | DUP`.
pub const MOSDEPTH_DEFAULT_EXCLUDE: u16 = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;

/// Accumulates per-base depth for a single contig.
#[derive(Debug)]
pub struct DepthAccum {
    /// Delta array of length `contig_len + 1`; a `+1` at a block start and a
    /// `-1` one past its end, summed into depth by [`DepthAccum::into_depths`].
    deltas: Vec<i32>,
    /// Contig length in bases.
    len: usize,
    /// Records with `MAPQ` strictly below this value are ignored.
    mapq_cut: u8,
    /// Records carrying any of these flag bits are ignored.
    exclude_flags: u16,
}

impl DepthAccum {
    /// Allocate for one contig of `length` bases.
    pub fn new(length: u64, mapq_cut: u8, exclude_flags: u16) -> Self {
        let len = length as usize;
        Self {
            deltas: vec![0i32; len + 1],
            len,
            mapq_cut,
            exclude_flags,
        }
    }

    /// Add one record's aligned blocks. Records failing the filters are ignored.
    pub fn process_read(&mut self, record: &bam::Record) {
        if !self.passes_filters(record) {
            return;
        }
        let mut pos = record.pos();
        for op in record.cigar().iter() {
            match op {
                // Reference-consuming and query-consuming: covers the reference.
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    let n = i64::from(*n);
                    self.add_block(pos, pos + n);
                    pos += n;
                }
                // Reference-consuming only: advances without covering.
                Cigar::Del(n) | Cigar::RefSkip(n) => pos += i64::from(*n),
                // Neither reference-consuming nor covering.
                Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
    }

    /// Consume the delta array and return per-base depth for the contig.
    pub fn into_depths(self) -> Vec<u32> {
        let mut depths = Vec::with_capacity(self.len);
        let mut running = 0i32;
        for delta in self.deltas.iter().take(self.len) {
            running += delta;
            // `running` cannot go negative: every `-1` is emitted only after
            // its matching `+1`, and both are clamped to the same range.
            depths.push(running.max(0) as u32);
        }
        depths
    }

    /// Whether a record contributes to depth at all.
    fn passes_filters(&self, record: &bam::Record) -> bool {
        record.flags() & self.exclude_flags == 0 && record.mapq() >= self.mapq_cut
    }

    /// Record a half-open aligned block `[start, end)`, clamped to the contig.
    fn add_block(&mut self, start: i64, end: i64) {
        let start = start.max(0) as usize;
        let end = (end.max(0) as usize).min(self.len);
        if start >= end {
            return;
        }
        self.deltas[start] += 1;
        self.deltas[end] -= 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};

    /// Build a minimal mapped record at `pos` with the given CIGAR, MAPQ and flags.
    ///
    /// `seq` and `qual` must both be as long as the query-consuming part of the
    /// CIGAR, otherwise the record is malformed and every assertion made against
    /// it is meaningless, so the helper asserts that itself.
    fn rec(pos: i64, cigar: Vec<Cigar>, mapq: u8, flags: u16) -> Record {
        let query_len: usize = cigar
            .iter()
            .map(|op| match op {
                Cigar::Match(n) | Cigar::Ins(n) | Cigar::SoftClip(n) => *n as usize,
                Cigar::Equal(n) | Cigar::Diff(n) => *n as usize,
                _ => 0,
            })
            .sum();
        let seq = vec![b'A'; query_len];
        let qual = vec![30u8; query_len];
        assert_eq!(seq.len(), qual.len(), "malformed test record");

        let mut r = Record::new();
        r.set(b"q", Some(&CigarString(cigar)), &seq, &qual);
        r.set_tid(0);
        r.set_pos(pos);
        r.set_mapq(mapq);
        r.set_flags(flags);
        r
    }

    #[test]
    fn match_block_covers_exactly_its_span() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(5, vec![Cigar::Match(4)], 60, 0));
        assert_eq!(&d.into_depths()[4..10], &[0, 1, 1, 1, 1, 0]);
    }

    #[test]
    fn deletion_and_skip_advance_without_covering() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(
            0,
            vec![Cigar::Match(2), Cigar::Del(3), Cigar::Match(2)],
            60,
            0,
        ));
        assert_eq!(&d.into_depths()[0..8], &[1, 1, 0, 0, 0, 1, 1, 0]);
    }

    #[test]
    fn ref_skip_advances_without_covering() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(
            0,
            vec![Cigar::Match(2), Cigar::RefSkip(3), Cigar::Match(2)],
            60,
            0,
        ));
        assert_eq!(&d.into_depths()[0..8], &[1, 1, 0, 0, 0, 1, 1, 0]);
    }

    #[test]
    fn insertion_and_soft_clip_do_not_advance_the_reference() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(
            0,
            vec![
                Cigar::SoftClip(3),
                Cigar::Match(2),
                Cigar::Ins(4),
                Cigar::Match(2),
            ],
            60,
            0,
        ));
        assert_eq!(&d.into_depths()[0..6], &[1, 1, 1, 1, 0, 0]);
    }

    #[test]
    fn duplicate_flagged_reads_are_excluded_by_default() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(0, vec![Cigar::Match(4)], 60, BAM_FDUP));
        assert_eq!(d.into_depths().iter().sum::<u32>(), 0);
    }

    #[test]
    fn secondary_qcfail_and_unmapped_reads_are_excluded_by_default() {
        for flag in [BAM_FSECONDARY, BAM_FQCFAIL, BAM_FUNMAP] {
            let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
            d.process_read(&rec(0, vec![Cigar::Match(4)], 60, flag));
            assert_eq!(
                d.into_depths().iter().sum::<u32>(),
                0,
                "flag {flag:#x} should be excluded"
            );
        }
    }

    #[test]
    fn reads_below_the_mapq_cutoff_are_excluded() {
        let mut d = DepthAccum::new(20, 30, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(0, vec![Cigar::Match(4)], 29, 0));
        assert_eq!(d.into_depths().iter().sum::<u32>(), 0);

        let mut d = DepthAccum::new(20, 30, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(0, vec![Cigar::Match(4)], 30, 0));
        assert_eq!(d.into_depths().iter().sum::<u32>(), 4);
    }

    #[test]
    fn a_read_running_past_the_contig_end_is_clipped_not_panicking() {
        let mut d = DepthAccum::new(6, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(4, vec![Cigar::Match(10)], 60, 0));
        assert_eq!(d.into_depths(), vec![0, 0, 0, 0, 1, 1]);
    }
}
