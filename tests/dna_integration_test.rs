//! Parity tests for the DNA pipeline against the committed reference outputs.
//!
//! The fixtures under `tests/expected/dna/` are the output of the upstream
//! tools themselves, at the versions pinned in `VERSIONS.txt`. A failure here
//! is a defect in RustQC, not a reason to regenerate the fixture.
//!
//! Compressed outputs are compared on their decompressed bytes. Two bgzf
//! writers at the same compression level need not emit identical compressed
//! bytes, so comparing the `.gz` files directly would test the compressor
//! rather than this code.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use rust_htslib::bam::Read as BamRead;
use rust_htslib::{bam, bgzf};

use rustqc::dna::depth::{DepthAccum, MOSDEPTH_DEFAULT_EXCLUDE};
use rustqc::dna::mosdepth::{output, ContigDepth, MosdepthResult};

/// Window size and thresholds the fixtures were generated with.
const WINDOW_SIZE: u32 = 500;
const THRESHOLDS: [u32; 7] = [1, 5, 10, 15, 20, 30, 50];

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/expected/dna")
        .join(name)
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join("rustqc-dna-parity");
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

/// Run the depth engine over the committed test BAM and summarise it exactly
/// as the fixtures were generated.
fn compute() -> MosdepthResult {
    let bam_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/dna/test.dna.bam");
    let mut reader = bam::Reader::from_path(&bam_path).unwrap();
    let header = reader.header().to_owned();

    let mut contigs = Vec::new();
    for tid in 0..header.target_count() {
        let name = String::from_utf8(header.tid2name(tid).to_vec()).unwrap();
        let length = header.target_len(tid).unwrap();
        let mut accum = DepthAccum::new(length, 0, MOSDEPTH_DEFAULT_EXCLUDE);

        let mut record = bam::Record::new();
        let mut per_contig = bam::Reader::from_path(&bam_path).unwrap();
        while let Some(result) = per_contig.read(&mut record) {
            result.unwrap();
            if record.tid() == tid as i32 {
                accum.process_read(&record);
            }
        }
        let depths = accum.into_depths();
        contigs.push(ContigDepth::from_depths(
            &name,
            &depths,
            Some(WINDOW_SIZE),
            &THRESHOLDS,
        ));
    }

    MosdepthResult {
        contigs,
        window_size: Some(WINDOW_SIZE),
        thresholds: THRESHOLDS.to_vec(),
    }
}

fn read_bgzf(path: &Path) -> String {
    let mut reader = bgzf::Reader::from_path(path).unwrap();
    let mut buf = Vec::new();
    reader.read_to_end(&mut buf).unwrap();
    String::from_utf8(buf).unwrap()
}

/// Compare line by line so a failure names the offending row.
fn assert_same_lines(actual: &str, expected: &str, what: &str) {
    let a: Vec<&str> = actual.lines().collect();
    let e: Vec<&str> = expected.lines().collect();
    for (i, (got, want)) in a.iter().zip(e.iter()).enumerate() {
        assert_eq!(got, want, "{what}: line {} differs", i + 1);
    }
    assert_eq!(a.len(), e.len(), "{what}: line count differs");
}

#[test]
fn fixture_tool_versions_are_the_pinned_ones() {
    let versions = std::fs::read_to_string(fixture("VERSIONS.txt")).unwrap();
    assert!(
        versions.contains("mosdepth\t0.3.14"),
        "unexpected mosdepth fixture version: {versions}"
    );
    assert!(
        versions.contains("samtools\t1.24"),
        "unexpected samtools fixture version: {versions}"
    );
}

#[test]
fn summary_matches_mosdepth() {
    let path = scratch("test.mosdepth.summary.txt");
    output::write_summary(&compute(), &path).unwrap();
    assert_same_lines(
        &std::fs::read_to_string(&path).unwrap(),
        &std::fs::read_to_string(fixture("test.mosdepth.summary.txt")).unwrap(),
        "summary",
    );
}

#[test]
fn global_dist_matches_mosdepth() {
    let path = scratch("test.mosdepth.global.dist.txt");
    output::write_global_dist(&compute(), &path).unwrap();
    assert_same_lines(
        &std::fs::read_to_string(&path).unwrap(),
        &std::fs::read_to_string(fixture("test.mosdepth.global.dist.txt")).unwrap(),
        "global dist",
    );
}

#[test]
fn region_dist_matches_mosdepth() {
    let path = scratch("test.mosdepth.region.dist.txt");
    output::write_region_dist(&compute(), &path).unwrap();
    assert_same_lines(
        &std::fs::read_to_string(&path).unwrap(),
        &std::fs::read_to_string(fixture("test.mosdepth.region.dist.txt")).unwrap(),
        "region dist",
    );
}

#[test]
fn per_base_matches_mosdepth() {
    let path = scratch("test.per-base.bed.gz");
    output::write_per_base(&compute(), &path).unwrap();
    assert_same_lines(
        &read_bgzf(&path),
        &read_bgzf(&fixture("test.per-base.bed.gz")),
        "per-base",
    );
}

#[test]
fn regions_match_mosdepth() {
    let path = scratch("test.regions.bed.gz");
    output::write_regions(&compute(), &path).unwrap();
    assert_same_lines(
        &read_bgzf(&path),
        &read_bgzf(&fixture("test.regions.bed.gz")),
        "regions",
    );
}

#[test]
fn thresholds_match_mosdepth() {
    let path = scratch("test.thresholds.bed.gz");
    output::write_thresholds(&compute(), &path).unwrap();
    assert_same_lines(
        &read_bgzf(&path),
        &read_bgzf(&fixture("test.thresholds.bed.gz")),
        "thresholds",
    );
}

/// The depth histogram is the input to both distribution files, so pinning it
/// separately makes a distribution failure easy to attribute.
#[test]
fn depth_histogram_matches_the_per_base_fixture() {
    let result = compute();
    let mut expected: BTreeMap<u32, u64> = BTreeMap::new();
    for line in read_bgzf(&fixture("test.per-base.bed.gz")).lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let start: u64 = fields[1].parse().unwrap();
        let end: u64 = fields[2].parse().unwrap();
        let depth: u32 = fields[3].parse().unwrap();
        *expected.entry(depth).or_insert(0) += end - start;
    }
    assert_eq!(result.contigs[0].histogram, expected);
}
