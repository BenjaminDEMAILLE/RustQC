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
  dna/
    mod.rs
    accumulators.rs                   per-read dispatch, one set per contig worker, merged at the end
    depth.rs                          full-contig i32 delta array plus prefix sum
    mosdepth/
      mod.rs
      output.rs                       per-base, regions, thresholds, summary, distributions
    wgs_metrics.rs
    insert_size.rs
    gc_bias.rs
    hs_metrics.rs
    qualimap/
      mod.rs
      report.rs
      plots.rs
    plots.rs                          coverage, insert size and GC curves via plotters
```

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
`min(threads, available_ram / (largest_contig_length * 4 bytes))` with a floor of
1. The cap applies only to the depth stage; read-level accumulators are cheap and
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

`summary.rs` gains a `DnaSummary` variant hung off `InputSummary`, holding mean
and median coverage, the threshold percentages, duplicate rate, insert size
median and MAD, and, in targeted mode, fold enrichment and fold-80 base penalty.
`CountingSummary` stays RNA-specific; the two are mutually exclusive per input.

`CITATIONS.md` gains the DNA tool citations with the pinned upstream versions.

## 6. Test data and parity

Reference data is a slice of a public high-coverage NA12878 alignment (1000
Genomes) on chr20, together with the matching GRCh38 chr20 FASTA slice and a
targets BED derived from an exome capture kit restricted to that interval.

Decision on size: the slice is trimmed to roughly 200 kb so the committed
fixtures stay under 10 MB in total. No git-lfs. If a metric turns out to need
more reads to be meaningful (preseq extrapolation is the likely candidate), the
interval is widened only as far as the 10 MB budget allows and the affected test
is marked `#[ignore]` with an explanation rather than the budget being raised.

`tests/create_dna_test_data.sh` downloads, slices, duplicate-marks and indexes
the inputs, converts the targets BED to a Picard `.interval_list`, and then runs
every upstream tool to regenerate `tests/expected/dna/`. The script pins tool
versions explicitly (mosdepth, samtools, Picard, Qualimap, preseq) and records
them in a `tests/expected/dna/VERSIONS.txt` that the test suite asserts against,
so a fixture regenerated with a different upstream version fails loudly instead
of silently changing the baseline.

`tests/dna_integration_test.rs` compares RustQC output against those fixtures
field by field: exact equality for integer fields and for text formatting,
relative tolerance of 1e-6 for floating point fields. Comparisons are structured
per tool so a failure names the offending column rather than diffing whole files.

Unit tests cover the pieces that fixtures cannot pin down cleanly: CIGAR walking
into the delta array, prefix-sum correctness, threshold counting at boundaries,
BED interval parsing and merging, and the depth cap.

## 7. PR stack

All branches are based on `main` and stacked in order.

1. **`feat/dna-common-extract`** — move `bam_flags`, `cpp_rng`, `preseq` and the
   samtools trio into `src/common/`, add re-export shims in `rna`, update
   `lib.rs` and `AGENTS.md`. No behaviour change; the existing test suite must
   pass untouched. This PR is deliberately mechanical so the later diffs are
   readable.
2. **`feat/dna-skeleton`** — `DnaArgs`, `DnaConfig`, `run_dna`, the contig-worker
   depth engine, samtools trio and preseq wired up for DNA, the test dataset and
   the generation script, plus parity tests against mosdepth and samtools.
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

## 9. Open questions

None blocking. The three points raised during design are resolved above: fixture
size budget (section 6), pinned upstream versions (section 6), and the RAM guard
(section 4).
