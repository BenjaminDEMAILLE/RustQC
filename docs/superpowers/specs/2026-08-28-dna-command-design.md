# Design: `rustqc dna` subcommand

Date: 2026-08-28
Status: approved for planning
Scope: new DNA (WGS and targeted) QC subcommand, delivered as a stack of five PRs

## 1. Goal

Add a `dna` subcommand that is to DNA sequencing what `rustqc rna` is to RNA-Seq:
a single-pass, single-binary reimplementation of the QC tools people otherwise
chain together, producing byte-compatible output files.

Tools in scope:

| Upstream tool | Module | Notes |
| --- | --- | --- |
| samtools `stats`, `flagstat`, `idxstats` | `common/samtools/` | already implemented for `rna`, moved to a shared module |
| mosdepth | `dna/mosdepth/` | per-base, per-window, per-region depth, distributions, thresholds |
| Picard `CollectWgsMetrics` | `dna/wgs_metrics.rs` | needs reference FASTA (genome territory excludes N bases) |
| Picard `CollectInsertSizeMetrics` | `dna/insert_size.rs` | per pair orientation, with histogram and plot |
| Picard `CollectGcBiasMetrics` | `dna/gc_bias.rs` | needs reference FASTA |
| Picard `CollectHsMetrics` | `dna/hs_metrics.rs` | targeted mode only, activated by `--targets` |
| Qualimap `bamqc` | `dna/qualimap/` | `genome_results.txt`, raw data tables, HTML report |
| preseq `lc_extrap` | `common/preseq.rs` | already implemented for `rna`, moved to a shared module |

Out of scope for this stack: contamination estimation (VerifyBamID), sex or ploidy
inference, variant-aware metrics, CRAM-specific optimisations beyond what
`rust-htslib` already provides.

## 2. Command-line surface

```
rustqc dna <INPUT>... [OPTIONS]
```

No `--gtf`, no `--stranded`. Everything else that `rna` exposes keeps the same
long name, the same short flag and the same `RUSTQC_*` environment variable, so
that muscle memory and existing wrapper scripts carry over:

`-o/--outdir`, `--sample-name`, `--flat-output`, `-c/--config`,
`-j/--json-summary`, `-t/--threads`, `-Q/--mapq`, `-p/--paired`, `-q/--quiet`,
`-v/--verbose`, `--skip-dup-check`, `--skip-preseq`, `--preseq-seed`,
`--preseq-max-extrap`, `--preseq-step-size`, `--preseq-n-bootstraps`,
`--preseq-seg-len`.

One shared flag takes a different default: **`-Q/--mapq` defaults to 0 for
`dna`**, not to 30 as it does for `rna`. That is mosdepth's `-Q` default, and
matching it is a precondition for exact parity. Changing it silently would make
every depth figure disagree with the reference tool.

New flags:

| Flag | Env | Default | Purpose |
| --- | --- | --- | --- |
| `-r/--reference <FASTA>` | `RUSTQC_REFERENCE` | none | required for CRAM input, for `CollectWgsMetrics` and for `CollectGcBiasMetrics` |
| `--targets <BED>` | `RUSTQC_TARGETS` | none | presence switches the run to targeted mode and enables `CollectHsMetrics` |
| `--baits <BED>` | `RUSTQC_BAITS` | value of `--targets` | capture bait intervals for `CollectHsMetrics` |
| `--depth-thresholds <N,...>` | `RUSTQC_DEPTH_THRESHOLDS` | `1,5,10,15,20,30,50` | mosdepth `--thresholds` and the `PCT_xX` columns of `CollectWgsMetrics` |
| `--window-size <N>` | `RUSTQC_WINDOW_SIZE` | none | mosdepth `--by <N>` fixed-width windows |
| `--coverage-cap <N>` | `RUSTQC_COVERAGE_CAP` | `250` | Picard `COVERAGE_CAP` |
| `--min-base-quality <N>` | `RUSTQC_MIN_BASE_QUALITY` | `20` | Picard `MINIMUM_BASE_QUALITY` |
| `--skip-per-base` | `RUSTQC_SKIP_PER_BASE` | `false` | suppress `*.per-base.bed.gz`, by far the largest output |
| `--skip-gc-bias` | `RUSTQC_SKIP_GC_BIAS` | `false` | skip `CollectGcBiasMetrics` |
| `--max-depth-workers <N>` | `RUSTQC_MAX_DEPTH_WORKERS` | derived | RAM guard for the depth engine, see section 4 |

Behaviour when `--reference` is absent: `CollectWgsMetrics` and
`CollectGcBiasMetrics` are skipped with a warning through `Ui`, the rest of the
pipeline runs. CRAM input without `--reference` is a hard error, matching `rna`.

`--baits` without `--targets` is a hard error.

### Configuration file

`config::Config` gains a `dna: DnaConfig` field alongside the existing `rna`.
`DnaConfig` mirrors `RnaConfig`: shared settings at the top level
(`chromosome_prefix`, `chromosome_mapping`, `flat_output`, each declared on
`DnaConfig` itself, since `RnaConfig` declares its own rather than inheriting
from the root `Config`) and one nested struct per tool, each with an `enabled`
toggle and its tool-specific parameters:

```yaml
dna:
  mosdepth:
    enabled: true
    window_size: 500
    thresholds: [1, 10, 30]
  wgs_metrics:
    enabled: true
    coverage_cap: 250
  insert_size: { enabled: true }
  gc_bias: { enabled: true }
  hs_metrics: { enabled: true }
  qualimap: { enabled: true }
  preseq: { enabled: false }
```

CLI flags override config values, config values override defaults, matching the
precedence `rna` already implements.

## 3. Module layout

Shared, non-RNA-specific code moves out of `src/rna/` into a new `src/common/`:

```
src/
  common/
    mod.rs
    bam_flags.rs                      moved from rna/bam_flags.rs
    bam_stat.rs                       moved from rna/rseqc/bam_stat.rs
    bam_stat_accum.rs                 BamStatAccum, lifted out of rna/rseqc/accumulators.rs
    cpp_rng.rs                        moved from rna/cpp_rng.rs
    preseq.rs                         moved from rna/preseq.rs
    samtools/
      mod.rs
      stats.rs                        moved from rna/rseqc/stats.rs
      flagstat.rs                     moved from rna/rseqc/flagstat.rs
      idxstats.rs                     moved from rna/rseqc/idxstats.rs
  rna/                                unchanged behaviour; re-exports the moved
                                      items so the 0.2.x library surface and all
                                      `crate::rna::...` paths keep working
  dna/                                built in PR2 unless marked
    mod.rs
    depth.rs                          full-contig i32 delta array, CIGAR walk,
                                      mate-overlap correction, prefix sum
    mosdepth/
      mod.rs                          per-contig summarisation, distribution rules
      output.rs                       per-base, regions, thresholds, summary, distributions
    wgs_metrics.rs                    PR3
    insert_size.rs                    PR3
    gc_bias.rs                        PR4
    hs_metrics.rs                     PR4
    qualimap/                         PR5
      mod.rs
      report.rs
      plots.rs
    plots.rs                          PR3 onwards: coverage, insert size and GC curves
```

The `accumulators.rs` sketched here was not built. Per-record dispatch turned
out to be three lines inside the contig worker in `run_dna`, feeding a
`DepthAccum`, a `BamStatAccum` and a `PreseqAccum` directly, so a module whose
only job was to forward one call to three others would have been indirection
for its own sake. Revisit that when PR3 and PR4 add accumulators with their own
filter regimes and the forwarding stops being trivial.

The samtools writers consume `bam_stat`'s result type, and the counters that
build it live in `BamStatAccum`, whose `process_read` takes only a record and a
MAPQ cutoff. Both are therefore assay-agnostic and move into `common` as well,
where the `dna` pipeline drives the same accumulator. This was found while
planning PR1 and is the one departure from the layout first sketched here.

`lib.rs` gains `pub mod common;` and `pub mod dna;`. The `rna` module keeps
`pub use crate::common::{bam_flags, cpp_rng, preseq};` and
`pub use crate::common::samtools::{flagstat, idxstats, stats};` inside
`rna::rseqc` so that no published path breaks. The re-export shims are documented
as deprecated-in-spirit and can be dropped at 1.0.

`main.rs` gains `run_dna(args, ui)` next to `run_rna`, and the `Verbosity` match
in `main()` is generalised over both subcommand variants rather than matching
`Commands::Rna` alone.

## 4. Single-pass engine

The DNA pipeline reuses the concurrency shape `rna` already uses: one rayon
worker per contig, each opening its own `bam::IndexedReader` and calling
`fetch(tid)`, with per-worker accumulators merged after the parallel section. A
final `fetch(FetchDefinition::Unmapped)` pass feeds the unmapped counters that
`flagstat` and `idxstats` need, exactly as `rna` does today.

### Depth representation

Each worker allocates a `Vec<i32>` the length of its contig. For every accepted
record the CIGAR is walked and reference-consuming aligned blocks contribute
`+1` at the block start and `-1` one past the block end. After the contig is
consumed a prefix sum turns the delta array into per-base depth, and a single
sweep over that array feeds, in one go: the depth histogram (mosdepth global and
region distributions, plus the `CollectWgsMetrics` histogram), the threshold
counters, the fixed-width windows, the BED regions and the `CollectHsMetrics`
per-target statistics.

This is the memory profile of mosdepth itself: roughly 1 GB for GRCh38 chr1
(248,956,422 bases at 4 bytes). With unbounded rayon parallelism that becomes
`threads x 1 GB`, so contigs are scheduled longest-first and the number of
concurrently live depth arrays is capped by `--max-depth-workers`, defaulting to
`min(threads, budget / (largest_contig_length * 4 bytes))` with a floor of
1. As built, `budget` is a fixed 4 GB constant: there is no portable way to ask
the operating system how much memory is free from std, and pulling in a
dependency for it was not worth doing when `--max-depth-workers` already gives
the user an exact override. The cap applies only to the depth stage; read-level accumulators are cheap and
stay fully parallel.

### One depth accumulator is not enough

mosdepth and Picard do not agree on which reads and which bases count. mosdepth
applies a flag filter and a mapping-quality floor and, outside fast mode, walks
internal CIGAR operations and reconciles overlapping mates.
`CollectWgsMetrics` additionally drops bases below `MINIMUM_BASE_QUALITY`,
applies its own mapping-quality floor and caps depth at `COVERAGE_CAP`, and
reports the excluded fractions as separate `PCT_EXC_*` columns.
`CollectHsMetrics` differs again.

Byte-exact parity is the stated acceptance criterion, so the design keeps
**separate depth accumulators per filter regime** rather than one shared array
with post-hoc corrections. Concretely: one accumulator for mosdepth semantics,
one for the Picard `WgsMetrics`/`HsMetrics` semantics. They are built in the same
contig worker from the same record stream, so the BAM is still read once.

The exact filter and CIGAR semantics of each upstream tool are **not** restated
here from memory. The first implementation task of each tool PR is to read the
pinned upstream source, write the semantics down as a table in the module doc
comment, and encode it in code. Guessing these rules is the single most likely
cause of parity failures.

## 5. Outputs

Nested subdirectories by default, all flattened into `--outdir` when
`--flat-output` is set, matching `rna`:

```
<outdir>/
  mosdepth/    <sample>.mosdepth.global.dist.txt, .mosdepth.summary.txt,
               .per-base.bed.gz(+.csi), .regions.bed.gz, .mosdepth.region.dist.txt,
               .thresholds.bed.gz
  picard/
    wgs_metrics/    <sample>.wgs_metrics.txt
    insert_size/    <sample>.insert_size_metrics.txt, .insert_size_histogram.svg
    gc_bias/        <sample>.gc_bias.detail_metrics.txt, .gc_bias.summary_metrics.txt, .gc_bias.svg
    hs_metrics/     <sample>.hs_metrics.txt, .per_target_coverage.txt
  samtools/    <sample>.stats.txt, .flagstat.txt, .idxstats.txt
  qualimap/    genome_results.txt, qualimapReport.html, raw_data_qualimapReport/
  preseq/      <sample>.lc_extrap.txt
  rustqc_summary.json
  CITATIONS.md
```

`summary.rs` gains a `DnaSummary` hung off `InputSummary`, holding genome
length, covered bases, mean, median and maximum coverage, the threshold
percentages, and the duplicate rate. Insert size figures join it in PR3, and
fold enrichment and fold-80 base penalty in PR4. The threshold percentages are
a list of `{threshold, pct_bases}` objects rather than a map, so the requested
order survives serialisation; a map keyed by the threshold would sort `"10"`
before `"5"`.
`CountingSummary` stays RNA-specific; the two are mutually exclusive per input.

`CITATIONS.md` gains the DNA tool citations with the pinned upstream versions.

## 6. Test data and parity

Reference data is a real public human alignment: the chr22 slice published in
**nf-core/test-datasets** (`data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam`,
40001 bp, 5644 records), together with the matching `genome.fasta` and its
index. The upstream BAM is not duplicate-marked, so the generation script marks
duplicates locally with `samtools markdup`, yielding 1656 duplicate-flagged
records.

This replaces the 1000 Genomes NA12878 chr20 slice originally planned here.
`ftp.1000genomes.ebi.ac.uk` does not resolve from the build sandbox, whereas
GitHub-hosted raw content does. The nf-core dataset is still a real public
human alignment, it is idiomatic for this project, and at roughly 240 kB of
inputs plus 140 kB of fixtures it sits far inside the size budget. A targets
BED for the eventual `CollectHsMetrics` work is derived from the same slice
when PR4 needs it.

Decision on size: committed fixtures stay under 10 MB in total. No git-lfs. If
a metric turns out to need more reads to be meaningful (preseq extrapolation is
the likely candidate), the dataset is widened only as far as the 10 MB budget
allows and the affected test is marked `#[ignore]` with an explanation rather
than the budget being raised.

Two properties of this dataset shape the tests, and both were measured rather
than assumed:

- **Mate pairs overlap almost completely.** Correcting for that takes total
  covered bases from 469875 to 247878, exactly the gap between
  `mosdepth --fast-mode` and its default. The dataset therefore exercises
  overlap correction hard, which is desirable: it is the single easiest thing
  to get wrong in the depth engine.
- **The MAPQ filter is not exercised.** `mosdepth -Q 30` returns the same
  totals as the default on this data, so that code path is covered by unit
  tests over synthetic records instead of by the fixtures.

`tests/create_dna_test_data.sh` downloads, duplicate-marks and indexes the
inputs, and then runs every upstream tool to regenerate `tests/expected/dna/`.
The script pins tool versions explicitly and refuses to run against a different
one, recording them in a `tests/expected/dna/VERSIONS.txt` that the test suite
asserts against, so a fixture regenerated with a different upstream version
fails loudly instead of silently changing the baseline.

`tests/dna_integration_test.rs` compares RustQC output against those fixtures
field by field: exact equality for integer fields and for text formatting,
relative tolerance of 1e-6 for floating point fields. Comparisons are structured
per tool so a failure names the offending column rather than diffing whole files.

Unit tests cover the pieces that fixtures cannot pin down cleanly: CIGAR walking
into the delta array, prefix-sum correctness, threshold counting at boundaries,
BED interval parsing and merging, and the depth cap.

## 7. PR stack

All branches are based on `main` and stacked in order.

1. **`feat/dna-common-extract`** (delivered, PR #152) — moved `bam_flags`,
   `cpp_rng`, `preseq`, `bam_stat`, `BamStatAccum` and the samtools trio into
   `src/common/`, with re-export shims in `rna`. No behaviour change: all 61
   output files of a `rustqc rna` run were compared byte for byte before and
   after.
2. **`feat/dna-skeleton`** (delivered, PR #153) — `DnaArgs`, `DnaConfig`,
   `run_dna`, the contig-worker depth engine with mate-overlap correction, the
   six mosdepth outputs, the samtools trio and preseq wired up for DNA, the
   test dataset and its generation script, the DNA JSON summary block and
   citations, and parity tests against mosdepth 0.3.14 and samtools 1.24. Every
   mosdepth output matches exactly; `samtools stats` matches on all 1889 data
   lines, its header differing by design.
3. **`feat/dna-picard-core`** — `CollectWgsMetrics` and
   `CollectInsertSizeMetrics`, their plots, and their parity fixtures.
4. **`feat/dna-targeted`** — `--targets` / `--baits`, `CollectHsMetrics`,
   `CollectGcBiasMetrics`, and targeted-mode parity fixtures.
5. **`feat/dna-qualimap-docs`** — Qualimap `bamqc` output and HTML report,
   `docs/src/content/docs/dna/*` pages, README, CHANGELOG and AGENTS.md updates.

Each PR is independently buildable, `cargo fmt --check` and
`cargo clippy -- -D warnings` clean, and ships its own tests.

## 8. Risks

- **Parity drift from upstream version skew.** Mitigated by pinning versions in
  the generation script and asserting them from the test suite.
- **Memory on real genomes.** The chosen engine holds a full-contig array per
  worker by explicit decision. `--max-depth-workers` bounds it, but a user who
  raises the cap on a 32-thread machine can still exhaust RAM. The default must
  be conservative and the flag documented with the arithmetic.
- **Picard interval semantics.** Picard consumes `.interval_list`, RustQC accepts
  BED. The one-based-inclusive versus zero-based-half-open conversion is a
  classic off-by-one; it gets its own unit tests and its own fixture.
- **Qualimap HTML.** The report is the least specified output and the least
  valuable to match byte for byte. Parity is asserted on `genome_results.txt` and
  the raw data tables only; the HTML is checked for structural presence, not
  equality.
- **Test fixture size budget.** Ten megabytes is tight for coverage-based
  metrics. Section 6 states what gives way if it binds.

## 9. Open items

The three points raised during design are resolved above: fixture size budget
(section 6), pinned upstream versions (section 6), and the RAM guard
(section 4). Three things surfaced during PR1 and PR2 and are still open:

- **`.csi` companion indexes** for the bgzf outputs are not written. mosdepth
  writes them, and `tabix` needs them. They are buildable through the
  `tbx_index_build` FFI that hts-sys already exposes.
- **The samtools citation constant says v1.22.1**, the version the `rna`
  pipeline was validated against, while the DNA fixtures were generated with
  1.24. Bumping it would imply the RNA pipeline had been revalidated, which it
  has not, so the constant was left alone and the discrepancy flagged instead.
- **`--targets` is accepted but inert** in PR2, warning that targeted metrics
  are not implemented. It becomes real in PR4.
