# DNA PR2: the `dna` subcommand skeleton and depth engine

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship `rustqc dna <BAM>...` running the depth engine, the mosdepth-compatible outputs, the samtools trio and preseq in a single pass, with exact-parity tests against mosdepth 0.3.14 and samtools 1.24.

**Architecture:** One rayon worker per contig, each holding a full-contig `i32` delta array (spec option B), merged after the parallel section. Read-level counters reuse `BamStatAccum` from PR1 unchanged. Picard and Qualimap outputs are explicitly out of scope; they land in PR3 to PR5.

**Tech Stack:** Rust 2021, rust-htslib (including its bgzf writer and the `tbx_index_build` FFI), rayon, anyhow, flate2 for reading fixtures in tests.

**Spec:** `docs/superpowers/specs/2026-08-28-dna-command-design.md`

## Global Constraints

- Rust edition 2021, `rust-version = "1.87"`. Do not raise either.
- Default `rustfmt`; `cargo clippy -- -D warnings` must pass; `anyhow::Result<T>` everywhere; `unwrap()` in tests only.
- Every file opens with `//!`; every public item gets `///`.
- No new third-party dependency without saying so in the PR description. `flate2` and `rust-htslib` are already present.
- Committed test fixtures stay under 10 MB in total.
- Base branch is `feat/dna-common-extract` (PR #152). Do not restart from `main`.

## Ground truth measured on the actual tools

These were measured, not recalled. Do not substitute remembered behaviour.

**mosdepth 0.3.14**, from `mosdepth --help`:
- `-F --flag` defaults to **1796** = `UNMAP(0x4) | SECONDARY(0x100) | QCFAIL(0x200) | DUP(0x400)`. A record is skipped when it has **any** of those bits.
- `-Q --mapq` defaults to **0**. Records with MAPQ strictly less than the threshold are ignored.
- `-x --fast-mode` is documented as "dont look at internal cigar operations or correct mate overlaps". Default mode therefore **does** both: it walks internal CIGAR operations, and it counts a base once when two mates of the same pair overlap it.

**Measured on `tests/data/dna/test.dna.bam`** (chr22 slice, 40001 bp, 5644 records, 1656 duplicate-flagged), total covered bases:

| invocation | bases | mean |
| --- | --- | --- |
| default | 247878 | 6.20 |
| `-x` (fast mode) | 469875 | 11.75 |
| `-F 0` | 354522 | 8.86 |
| `-Q 30` | 247878 | 6.20 |
| `-a` (fragment mode) | 247589 | 6.19 |

Read this table carefully before writing the engine:
- Default versus fast mode nearly doubles the count. In this dataset mate pairs overlap almost completely, so **mate-overlap correction is the dominant behaviour, not an edge case**. An engine that ignores it will be wrong by a factor of two, not by a rounding error.
- `-F 0` versus default shows duplicate exclusion is exercised.
- `-Q 30` matches the default exactly, so **this dataset does not exercise the MAPQ filter**. Cover that path with a unit test on a synthetic record instead of relying on the fixtures.

## File structure

| File | Responsibility |
| --- | --- |
| `src/dna/mod.rs` | module re-exports |
| `src/dna/depth.rs` | `DepthAccum`: contig delta array, CIGAR walk, mate-overlap correction, prefix sum |
| `src/dna/accumulators.rs` | `DnaAccumulators`: one `DepthAccum`, one `BamStatAccum` and one `PreseqAccum` per contig worker |
| `src/dna/mosdepth/mod.rs` | result types produced from the merged accumulators |
| `src/dna/mosdepth/output.rs` | the six mosdepth output writers plus the bgzf and CSI plumbing |
| `src/cli.rs` | `Commands::Dna(DnaArgs)` and `DnaArgs` |
| `src/config.rs` | `DnaConfig` and its per-tool nested structs |
| `src/main.rs` | `run_dna` orchestration, verbosity match generalised over both subcommands |
| `src/summary.rs` | `DnaSummary` |
| `tests/create_dna_test_data.sh` | regenerates inputs and fixtures from nf-core/test-datasets |
| `tests/dna_integration_test.rs` | parity assertions against the committed fixtures |

## Deviation from the spec, section 6

The spec named a 1000 Genomes NA12878 chr20 slice. `ftp.1000genomes.ebi.ac.uk` does not resolve from the build sandbox, so the dataset is instead the **nf-core/test-datasets** chr22 slice (`data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam`, GitHub-hosted and therefore reachable), duplicate-marked locally with `samtools markdup`. It is still a real public human alignment, it is idiomatic for this project, and at roughly 230 kB of inputs plus 120 kB of fixtures it sits far inside the 10 MB budget. Update spec section 6 as part of Task 1.

---

### Task 1: Commit the test dataset, its generation script and the fixtures

**Files:**
- Create: `tests/create_dna_test_data.sh`
- Create: `tests/data/dna/{test.dna.bam,test.dna.bam.bai,genome.fasta,genome.fasta.fai}`
- Create: `tests/expected/dna/` (mosdepth and samtools outputs, plus `VERSIONS.txt`)
- Modify: `docs/superpowers/specs/2026-08-28-dna-command-design.md` (section 6)

**Interfaces:**
- Produces: the fixture paths that every later task asserts against.

- [ ] **Step 1: Write the generation script**

```bash
#!/usr/bin/env bash
# Regenerate the DNA test inputs and the reference outputs they are compared against.
#
# Inputs come from nf-core/test-datasets (a real human chr22 slice, 40 kb).
# The upstream BAM is not duplicate-marked, so this script marks duplicates
# with samtools; RustQC requires duplicate-marked input.
#
# The reference outputs are produced by the pinned tool versions recorded in
# tests/expected/dna/VERSIONS.txt. Regenerating with a different version will
# make the parity tests fail, which is the intended behaviour: fixtures and
# tool versions travel together.
set -euo pipefail

MOSDEPTH_VERSION="0.3.14"
SAMTOOLS_VERSION="1.24"

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data="$here/data/dna"
expected="$here/expected/dna"
base="https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens"

have() { command -v "$1" >/dev/null || { echo "missing tool: $1" >&2; exit 1; }; }
have samtools; have mosdepth; have curl

check_version() {
  local tool="$1" want="$2" got
  got="$($tool --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)"
  if [[ "$got" != "$want" ]]; then
    echo "$tool version $got does not match the pinned $want" >&2
    echo "Install the pinned version, or update VERSIONS.txt and the fixtures together." >&2
    exit 1
  fi
}
check_version samtools "$SAMTOOLS_VERSION"
check_version mosdepth "$MOSDEPTH_VERSION"

mkdir -p "$data" "$expected"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

curl -sSfL -o "$tmp/upstream.bam"      "$base/illumina/bam/test.paired_end.sorted.bam"
curl -sSfL -o "$data/genome.fasta"     "$base/genome/genome.fasta"
curl -sSfL -o "$data/genome.fasta.fai" "$base/genome/genome.fasta.fai"

# Mark duplicates: name-sort, add mate tags, coordinate-sort, then markdup.
samtools sort -n -o "$tmp/ns.bam" "$tmp/upstream.bam"
samtools fixmate -m "$tmp/ns.bam" "$tmp/fm.bam"
samtools sort -o "$tmp/cs.bam" "$tmp/fm.bam"
samtools markdup -S "$tmp/cs.bam" "$data/test.dna.bam"
samtools index "$data/test.dna.bam"

mosdepth --by 500 --thresholds 1,5,10,15,20,30,50 "$expected/test" "$data/test.dna.bam"
samtools stats    "$data/test.dna.bam" > "$expected/test.stats.txt"
samtools flagstat "$data/test.dna.bam" > "$expected/test.flagstat.txt"
samtools idxstats "$data/test.dna.bam" > "$expected/test.idxstats.txt"

printf 'mosdepth\t%s\nsamtools\t%s\n' "$MOSDEPTH_VERSION" "$SAMTOOLS_VERSION" > "$expected/VERSIONS.txt"

echo "Regenerated $(find "$data" "$expected" -type f | wc -l | tr -d ' ') files."
```

Write that to `tests/create_dna_test_data.sh` and `chmod +x` it.

- [ ] **Step 2: Run it**

```bash
./tests/create_dna_test_data.sh
```

Expected: it reports the file count and exits 0. If a version check fails, install the pinned tool rather than loosening the check.

- [ ] **Step 3: Verify the dataset matches the ground-truth table**

```bash
awk '$1=="total"' tests/expected/dna/test.mosdepth.summary.txt
samtools flagstat tests/data/dna/test.dna.bam | sed -n '1,6p'
du -sh tests/data/dna tests/expected/dna
```

Expected: `total 40001 247878 6.20 0 867`, 5644 records with 1656 duplicates, and well under 10 MB. If the numbers differ, the upstream dataset moved; stop and report rather than adjusting the table.

- [ ] **Step 4: Update spec section 6**

Replace the NA12878 paragraph with the nf-core/test-datasets description given under "Deviation from the spec" above, keeping the 10 MB budget rule and the pinned-versions rule as written.

- [ ] **Step 5: Commit**

```bash
git add tests/create_dna_test_data.sh tests/data/dna tests/expected/dna
git commit -m "test: add the DNA test dataset and its reference outputs

A real public human chr22 slice from nf-core/test-datasets, duplicate-marked
locally, with mosdepth 0.3.14 and samtools 1.24 reference outputs. The
generation script pins both tool versions and refuses to run against others,
so fixtures and tool versions cannot drift apart.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: The `dna` CLI surface

**Files:**
- Modify: `src/cli.rs`
- Test: `src/cli.rs` (the existing `mod tests`)

**Interfaces:**
- Produces: `cli::Commands::Dna(cli::DnaArgs)`; `DnaArgs` with the fields listed in spec section 2. Task 7 reads them.

- [ ] **Step 1: Write the failing tests**

Add to `mod tests` in `src/cli.rs`:

```rust
#[test]
fn test_dna_default_args() {
    let cli = Cli::parse_from(["rustqc", "dna", "test.bam"]);
    match cli.command {
        Commands::Dna(args) => {
            assert_eq!(args.input, vec!["test.bam"]);
            assert_eq!(args.outdir, ".");
            assert_eq!(args.threads, 1);
            assert_eq!(args.mapq_cut, 0);
            assert_eq!(args.coverage_cap, 250);
            assert_eq!(args.min_base_quality, 20);
            assert_eq!(args.depth_thresholds, vec![1, 5, 10, 15, 20, 30, 50]);
            assert_eq!(args.window_size, None);
            assert!(args.targets.is_none());
            assert!(args.baits.is_none());
            assert!(!args.skip_per_base);
            assert!(!args.skip_gc_bias);
        }
        _ => panic!("Expected Dna subcommand"),
    }
}

#[test]
fn test_dna_no_gtf_required() {
    assert!(Cli::try_parse_from(["rustqc", "dna", "test.bam"]).is_ok());
}

#[test]
fn test_dna_targeted_args() {
    let cli = Cli::parse_from([
        "rustqc", "dna", "a.bam", "b.bam",
        "--targets", "t.bed",
        "--baits", "b.bed",
        "--depth-thresholds", "1,10,100",
        "--window-size", "500",
        "--reference", "genome.fa",
        "-Q", "20",
        "--threads", "4",
    ]);
    match cli.command {
        Commands::Dna(args) => {
            assert_eq!(args.input, vec!["a.bam", "b.bam"]);
            assert_eq!(args.targets, Some("t.bed".to_string()));
            assert_eq!(args.baits, Some("b.bed".to_string()));
            assert_eq!(args.depth_thresholds, vec![1, 10, 100]);
            assert_eq!(args.window_size, Some(500));
            assert_eq!(args.reference, Some("genome.fa".to_string()));
            assert_eq!(args.mapq_cut, 20);
            assert_eq!(args.threads, 4);
        }
        _ => panic!("Expected Dna subcommand"),
    }
}

#[test]
fn test_dna_baits_without_targets_is_rejected() {
    let result = Cli::try_parse_from(["rustqc", "dna", "test.bam", "--baits", "b.bed"]);
    assert!(result.is_err(), "--baits without --targets must be rejected");
}
```

The existing `rna` tests match on `Commands::Rna(args)` with an `#[allow(unreachable_patterns)] _ =>` arm. Once a second variant exists that attribute becomes unnecessary; remove it from those tests so clippy stays quiet.

- [ ] **Step 2: Run the tests to verify they fail**

```bash
cargo test --lib test_dna 2>&1 | tail -20
```

Expected: compile error, no variant named `Dna`.

- [ ] **Step 3: Add the subcommand variant**

In `enum Commands`:

```rust
    /// DNA QC — single-pass analysis of BAM/SAM/CRAM files.
    ///
    /// Runs depth of coverage, samtools stats and library complexity
    /// estimation in one pass. Needs no annotation. Pass `--targets` to
    /// switch to targeted (exome or panel) mode.
    Dna(DnaArgs),
```

- [ ] **Step 4: Add the `DnaArgs` struct**

Mirror `RnaArgs`: the same `help_template`, the same help headings, the same short flags and `RUSTQC_*` environment variables for the shared options. Defaults come from the table in spec section 2. `mapq_cut` defaults to **0**, matching mosdepth's `-Q` default, not to 30 like `rna`.

```rust
/// Arguments for the `dna` subcommand.
#[derive(Parser, Debug)]
#[command(
    next_line_help = false,
    term_width = 120,
    help_template = "\
{about-with-newline}
{usage-heading} {usage}

{all-args}"
)]
pub struct DnaArgs {
    // ── Input / Output ──────────────────────────────────────────────────
    /// Duplicate-marked alignment file(s)
    #[arg(value_name = "INPUT", num_args = 1.., required = true, help_heading = "Input / Output")]
    pub input: Vec<String>,

    /// Reference FASTA (required for CRAM, CollectWgsMetrics and GC bias)
    #[arg(short, long, value_name = "FASTA", env = "RUSTQC_REFERENCE", help_heading = "Input / Output")]
    pub reference: Option<String>,

    /// Target intervals BED; switches on targeted (exome or panel) mode
    #[arg(long, value_name = "BED", env = "RUSTQC_TARGETS", help_heading = "Input / Output")]
    pub targets: Option<String>,

    /// Capture bait intervals BED [default: same as --targets]
    #[arg(long, value_name = "BED", env = "RUSTQC_BAITS", requires = "targets", help_heading = "Input / Output")]
    pub baits: Option<String>,

    // ... the remaining shared fields, copied field for field from RnaArgs:
    // outdir, sample_name, flat_output, config, json_summary, threads, quiet,
    // verbose, skip_dup_check, paired, and the preseq_* family.

    // ── Tool parameters ─────────────────────────────────────────────────
    /// MAPQ cutoff; reads below it are ignored [default: 0]
    #[arg(short = 'Q', long = "mapq", default_value_t = 0, hide_default_value = true, env = "RUSTQC_MAPQ", help_heading = "General")]
    pub mapq_cut: u8,

    /// Coverage thresholds to report [default: 1,5,10,15,20,30,50]
    #[arg(
        long = "depth-thresholds",
        value_name = "N,...",
        value_delimiter = ',',
        default_values_t = vec![1u32, 5, 10, 15, 20, 30, 50],
        hide_default_value = true,
        env = "RUSTQC_DEPTH_THRESHOLDS",
        help_heading = "Tool parameters"
    )]
    pub depth_thresholds: Vec<u32>,

    /// Fixed-width window size for per-window depth
    #[arg(long = "window-size", value_name = "N", env = "RUSTQC_WINDOW_SIZE", help_heading = "Tool parameters")]
    pub window_size: Option<u32>,

    /// Picard COVERAGE_CAP [default: 250]
    #[arg(long = "coverage-cap", default_value_t = 250, hide_default_value = true, env = "RUSTQC_COVERAGE_CAP", help_heading = "Tool parameters")]
    pub coverage_cap: u32,

    /// Picard MINIMUM_BASE_QUALITY [default: 20]
    #[arg(long = "min-base-quality", default_value_t = 20, hide_default_value = true, env = "RUSTQC_MIN_BASE_QUALITY", help_heading = "Tool parameters")]
    pub min_base_quality: u8,

    /// Skip the per-base depth output, by far the largest file
    #[arg(long, default_value_t = false, env = "RUSTQC_SKIP_PER_BASE", help_heading = "Tool parameters")]
    pub skip_per_base: bool,

    /// Skip GC bias metrics
    #[arg(long, default_value_t = false, env = "RUSTQC_SKIP_GC_BIAS", help_heading = "Tool parameters")]
    pub skip_gc_bias: bool,

    /// Cap on concurrently live per-contig depth arrays [default: derived from RAM]
    #[arg(long = "max-depth-workers", value_name = "N", env = "RUSTQC_MAX_DEPTH_WORKERS", help_heading = "Tool parameters")]
    pub max_depth_workers: Option<usize>,
}
```

`requires = "targets"` on `baits` is what makes the fourth test pass; do not hand-roll that check.

- [ ] **Step 5: Run the tests**

```bash
cargo test --lib test_dna 2>&1 | tail -10
```

Expected: 4 passed. `cargo build` will still fail because `main()` does not handle the new variant. To keep this task independently green, add the match arm now as a stub returning `anyhow::bail!("the dna subcommand is not implemented yet")`, and generalise the `Verbosity` match in `main()` to cover both variants.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "feat(cli): add the dna subcommand surface

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: `DnaConfig`

**Files:**
- Modify: `src/config.rs`

**Interfaces:**
- Produces: `config::Config::dna: DnaConfig`, with `flat_output`, `chromosome_prefix` and `chromosome_mapping` declared on `DnaConfig` itself (as `RnaConfig` does), plus `mosdepth: MosdepthConfig`, `samtools: SamtoolsConfig` and `preseq: PreseqConfig` (the existing type, reused).

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn test_dna_config_defaults() {
    let config = Config::default();
    assert!(config.dna.mosdepth.enabled);
    assert!(config.dna.samtools.enabled);
    assert!(config.dna.preseq.enabled);
    assert!(!config.dna.flat_output);
    assert_eq!(config.dna.mosdepth.thresholds, vec![1, 5, 10, 15, 20, 30, 50]);
    assert_eq!(config.dna.mosdepth.window_size, None);
}

#[test]
fn test_dna_config_from_yaml() {
    let yaml = "dna:\n  flat_output: true\n  mosdepth:\n    window_size: 500\n    thresholds: [1, 30]\n  preseq:\n    enabled: false\n";
    let config: Config = serde_yaml_ng::from_str(yaml).unwrap();
    assert!(config.dna.flat_output);
    assert_eq!(config.dna.mosdepth.window_size, Some(500));
    assert_eq!(config.dna.mosdepth.thresholds, vec![1, 30]);
    assert!(!config.dna.preseq.enabled);
    // A dna-only config leaves the rna side untouched.
    assert!(config.rna.preseq.enabled);
}

#[test]
fn test_dna_config_deep_merge() {
    let mut merged: Value = serde_yaml_ng::from_str("dna:\n  mosdepth:\n    window_size: 100\n    thresholds: [1]\n").unwrap();
    let overlay: Value = serde_yaml_ng::from_str("dna:\n  mosdepth:\n    window_size: 500\n").unwrap();
    deep_merge(&mut merged, overlay);
    let config: Config = serde_yaml_ng::from_value(merged).unwrap();
    assert_eq!(config.dna.mosdepth.window_size, Some(500));
    assert_eq!(config.dna.mosdepth.thresholds, vec![1]);
}
```

- [ ] **Step 2: Run them to verify they fail**

```bash
cargo test --lib test_dna_config 2>&1 | tail -10
```

Expected: no field `dna` on `Config`.

- [ ] **Step 3: Implement**

Add `pub dna: DnaConfig` to `Config`, then `DnaConfig`, `MosdepthConfig` and `SamtoolsConfig` following the exact shape of `RnaConfig` and its nested structs: `#[derive(Debug, Deserialize, Default)]`, `#[serde(default)]`, a hand-written `Default` impl wherever a field's default is not the type's default, and a `///` doc comment carrying a YAML example. Reuse the existing `PreseqConfig`; do not fork it.

`MosdepthConfig::default()` sets `enabled: true`, `thresholds: vec![1, 5, 10, 15, 20, 30, 50]`, `window_size: None`, `skip_per_base: false`.

- [ ] **Step 4: Run the tests**

```bash
cargo test --lib test_dna_config 2>&1 | tail -10
```

Expected: 3 passed, and the whole `--lib` suite still green.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(config): add the dna configuration section

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: The depth engine

The core of the PR. Write the unit tests first; they are the only place the MAPQ path and the CIGAR edge cases get covered, because the fixture dataset does not exercise them.

**Files:**
- Create: `src/dna/mod.rs`, `src/dna/depth.rs`
- Modify: `src/lib.rs` (`pub mod dna;`)

**Interfaces:**
- Consumes: `crate::common::bam_flags::*`.
- Produces:

```rust
pub struct DepthAccum { /* private */ }

impl DepthAccum {
    /// Allocate for one contig of `length` bases.
    pub fn new(length: u64, mapq_cut: u8, exclude_flags: u16) -> Self;
    /// Add one record's aligned blocks. Records failing the filters are ignored.
    pub fn process_read(&mut self, record: &bam::Record);
    /// Consume the delta array and return per-base depth for the contig.
    pub fn into_depths(self) -> Vec<u32>;
}

/// Bit mask matching mosdepth's `-F` default: UNMAP | SECONDARY | QCFAIL | DUP.
pub const MOSDEPTH_DEFAULT_EXCLUDE: u16 = 1796;
```

- [ ] **Step 1: Write the failing tests**

Build a `rec(pos, cigar, mapq, flags)` helper first. `Record::set` takes `(qname, cigar, seq, qual)` and requires `seq` and `qual` to have the same length, equal to the number of query-consuming CIGAR bases (`Match`, `Ins`, `SoftClip`, `Equal`, `Diff`). Get the helper right before writing assertions against it: a helper that silently builds an invalid record makes every assertion meaningless. Assert in the helper itself that the two lengths agree.

```rust
#[test]
fn match_block_covers_exactly_its_span() {
    let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
    d.process_read(&rec(5, vec![Cigar::Match(4)], 60, 0));
    assert_eq!(&d.into_depths()[4..10], &[0, 1, 1, 1, 1, 0]);
}

#[test]
fn deletion_and_skip_advance_without_covering() {
    let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
    d.process_read(&rec(0, vec![Cigar::Match(2), Cigar::Del(3), Cigar::Match(2)], 60, 0));
    assert_eq!(&d.into_depths()[0..8], &[1, 1, 0, 0, 0, 1, 1, 0]);
}

#[test]
fn insertion_and_soft_clip_do_not_advance_the_reference() {
    let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
    d.process_read(&rec(0, vec![Cigar::SoftClip(3), Cigar::Match(2), Cigar::Ins(4), Cigar::Match(2)], 60, 0));
    assert_eq!(&d.into_depths()[0..6], &[1, 1, 1, 1, 0, 0]);
}

#[test]
fn duplicate_flagged_reads_are_excluded_by_default() {
    let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
    d.process_read(&rec(0, vec![Cigar::Match(4)], 60, BAM_FDUP));
    assert_eq!(d.into_depths().iter().sum::<u32>(), 0);
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
```

- [ ] **Step 2: Run them to verify they fail**

```bash
cargo test --lib dna::depth 2>&1 | tail -20
```

Expected: module not found.

- [ ] **Step 3: Implement `DepthAccum`**

`new` allocates `vec![0i32; length + 1]` so the `-1` at a block's end never needs a bounds branch. `process_read` returns early when `record.flags() & exclude_flags != 0` or `record.mapq() < mapq_cut`, then walks `record.cigar()` tracking the reference position: `Match | Equal | Diff` emit `+1` at the block start and `-1` at the block end and advance; `Del | RefSkip` advance without emitting; `Ins | SoftClip | HardClip | Pad` do not advance. Clamp both endpoints to the contig length. `into_depths` runs the prefix sum and drops the sentinel.

`i32` is the accumulator type because the array holds signed increments; the returned depths are `u32`.

- [ ] **Step 4: Run the tests**

```bash
cargo test --lib dna::depth 2>&1 | tail -10
```

Expected: 6 passed.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(dna): add the per-contig depth accumulator

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Mate-overlap correction

Kept separate from Task 4 because it is the single most likely thing to get wrong, and because the ground-truth table shows it changes the answer by a factor of two on the fixture data. A reviewer should be able to accept Task 4 and reject this one.

**Files:**
- Modify: `src/dna/depth.rs`

**Interfaces:**
- Produces: `DepthAccum::process_read` gains overlap correction; the public signature does not change.

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn overlapping_mates_cover_a_base_once() {
    let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
    let flags = BAM_FPAIRED | BAM_FPROPER_PAIR;
    let mut r1 = rec(0, vec![Cigar::Match(4)], 60, flags | BAM_FREAD1);
    let mut r2 = rec(0, vec![Cigar::Match(4)], 60, flags | BAM_FREAD2);
    r1.set_qname(b"pair1");
    r2.set_qname(b"pair1");
    d.process_read(&r1);
    d.process_read(&r2);
    assert_eq!(&d.into_depths()[0..5], &[1, 1, 1, 1, 0], "an overlapped base counts once");
}

#[test]
fn non_overlapping_mates_each_contribute() {
    let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
    let flags = BAM_FPAIRED | BAM_FPROPER_PAIR;
    let mut r1 = rec(0, vec![Cigar::Match(4)], 60, flags | BAM_FREAD1);
    let mut r2 = rec(8, vec![Cigar::Match(4)], 60, flags | BAM_FREAD2);
    r1.set_qname(b"pair2");
    r2.set_qname(b"pair2");
    d.process_read(&r1);
    d.process_read(&r2);
    assert_eq!(&d.into_depths()[0..13], &[1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0]);
}

#[test]
fn reads_from_different_pairs_at_the_same_locus_both_count() {
    let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
    let flags = BAM_FPAIRED | BAM_FPROPER_PAIR;
    let mut a = rec(0, vec![Cigar::Match(4)], 60, flags | BAM_FREAD1);
    let mut b = rec(0, vec![Cigar::Match(4)], 60, flags | BAM_FREAD1);
    a.set_qname(b"pairA");
    b.set_qname(b"pairB");
    d.process_read(&a);
    d.process_read(&b);
    assert_eq!(d.into_depths()[0], 2);
}

#[test]
fn the_pending_mate_map_does_not_leak() {
    let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
    let flags = BAM_FPAIRED | BAM_FPROPER_PAIR;
    let mut r1 = rec(0, vec![Cigar::Match(4)], 60, flags | BAM_FREAD1);
    let mut r2 = rec(0, vec![Cigar::Match(4)], 60, flags | BAM_FREAD2);
    r1.set_qname(b"pair1");
    r2.set_qname(b"pair1");
    d.process_read(&r1);
    d.process_read(&r2);
    assert_eq!(d.pending_mates_len(), 0, "both mates seen, the entry must be dropped");
}
```

`pending_mates_len` is a `#[cfg(test)]` accessor; do not widen the public API for it.

- [ ] **Step 2: Run them to verify they fail**

```bash
cargo test --lib dna::depth 2>&1 | tail -10
```

Expected: `overlapping_mates_cover_a_base_once` fails with `[2, 2, 2, 2, 0]`, proving the uncorrected engine double-counts.

- [ ] **Step 3: Implement**

A record can only overlap its own mate, and only when both are on the same contig. The worker walks one contig in coordinate order, so hold a map from read name to the mate's already-counted reference intervals, populated only for records whose mate is on the same contig at a position at or after the current one. When the second mate arrives, subtract the intersection of its aligned blocks with the stored intervals before adding, then drop the entry.

Bound the map: an entry whose stored mate position lies further behind the current position than the largest observed template length can never be claimed and must be evicted, otherwise a large contig leaks memory.

- [ ] **Step 4: Run the tests**

```bash
cargo test --lib dna::depth 2>&1 | tail -10
```

Expected: all green, Task 4's tests included.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(dna): correct mate overlaps in the depth engine

mosdepth counts a base once when both mates of a pair cover it, unless
--fast-mode is given. On the test dataset this halves total covered bases
(469875 uncorrected versus 247878 corrected), so it is the dominant
behaviour rather than an edge case.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: mosdepth-compatible outputs

**Files:**
- Create: `src/dna/mosdepth/mod.rs`, `src/dna/mosdepth/output.rs`

**Interfaces:**
- Consumes: `Vec<u32>` per-contig depths from `DepthAccum::into_depths`.
- Produces: `MosdepthResult` and writers for `{prefix}.mosdepth.summary.txt`, `{prefix}.mosdepth.global.dist.txt`, `{prefix}.mosdepth.region.dist.txt`, `{prefix}.per-base.bed.gz`, `{prefix}.regions.bed.gz` and `{prefix}.thresholds.bed.gz`.

- [ ] **Step 1: Read the fixture formats before writing any code**

```bash
head -3 tests/expected/dna/test.mosdepth.summary.txt
head -3 tests/expected/dna/test.mosdepth.global.dist.txt
gunzip -c tests/expected/dna/test.per-base.bed.gz | head -5
gunzip -c tests/expected/dna/test.regions.bed.gz | head -3
gunzip -c tests/expected/dna/test.thresholds.bed.gz | head -3
```

Record, in the module doc comment, the exact column set, the separator, the numeric formatting (the summary's `6.20` is two decimals) and whether each file carries a header. The fixtures, not this plan, are the specification.

Two formats were already reverse-engineered from the fixture and verified
against every row, so Task 6 does not need to rediscover them.

**`{prefix}.mosdepth.summary.txt`** has the header
`chrom	length	bases	mean	min	max`, then one row per contig, then one
`{contig}_region` row per contig when `--by` was given, then `total`, then
`total_region`. `mean` is `bases / length` to two decimals.

**`{prefix}.mosdepth.global.dist.txt`** and `.region.dist.txt` hold
`chrom	depth	proportion` rows in descending depth order, where `proportion`
is the fraction of that contig's bases at depth **at or above** `depth`,
formatted `%.2f`, ending at depth 0 with `1.00`. The per-contig block is
followed by the same table for `total`. Which depths get a row is the part
worth writing down:

- depths **0 through 300 always get a row**, even when no base sits at that
  exact depth (301 is mosdepth's internal fixed depth-array size);
- **above 300**, only depths where at least one base has exactly that depth;
- the **maximum observed depth never gets a row**. On the fixture the maximum
  is 867 and the first row is 866.

That rule was checked against all 547 chr22 rows of the fixture, and the
cumulative-proportion formula reproduces every value with no mismatch.

- [ ] **Step 2: Write the failing unit tests**

One test per writer: build a small hand-made depth vector, write to a scratch path under `std::env::temp_dir()`, assert the exact bytes. Include a test asserting `per-base.bed.gz` collapses runs of equal depth into one interval, since that is what keeps the file small.

- [ ] **Step 3: Run them to verify they fail, then implement**

Compression: use `rust_htslib::bgzf::Writer`, which is what mosdepth writes and what both `tabix` and `gunzip` read. Parity for the compressed outputs is asserted on the **decompressed** bytes: two bgzf writers at the same level need not emit identical compressed bytes, so comparing the `.gz` byte for byte is not a valid criterion. Say so in the test module doc comment.

The `.csi` companion indexes are built with `tbx_index_build` from the htslib FFI. Assert that the index file exists and that htslib can open it; do not compare it byte for byte.

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "feat(dna): add the mosdepth-compatible output writers

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: `run_dna` orchestration

**Files:**
- Create: `src/dna/accumulators.rs`
- Modify: `src/main.rs`, `src/summary.rs`, `src/citations.rs`

**Interfaces:**
- Consumes: everything above.
- Produces: a working `rustqc dna` writing the output tree from spec section 5, minus the Picard and Qualimap directories.

- [ ] **Step 1: Write `DnaAccumulators`**

One per contig worker, owning a `DepthAccum`, a `BamStatAccum` (from `crate::common::bam_stat_accum`, unchanged) and an optional `PreseqAccum`. `process_read` feeds all three from the one record; `merge` folds a finished worker into the running total. `BamStatAccum` and `PreseqAccum` already have `merge`; the depth results concatenate rather than merge, because each worker owns a distinct contig.

- [ ] **Step 2: Write `run_dna` in `main.rs`**

Model it on `run_rna`: resolve config, apply CLI overrides, create the output directory, then for each input open an `IndexedReader`, schedule contigs longest-first across a rayon pool bounded by `max_depth_workers`, run the unmapped-records pass for flagstat and idxstats, write outputs, and collect the JSON summary. Reject CRAM without `--reference`, exactly as `run_rna` does.

Default for `max_depth_workers`: `min(threads, budget / (largest_contig_len * 4))` with a floor of 1. There is no portable "available RAM" call in std, so use a conservative constant budget of 4 GB unless a dependency already exposes the real figure, and say which it is in the doc comment.

- [ ] **Step 3: Add `DnaSummary` to `summary.rs`**

`InputSummary` gains `#[serde(skip_serializing_if = "Option::is_none")] pub dna: Option<DnaSummary>`, holding mean and median coverage, the per-threshold percentages, the duplicate rate and the covered-bases total. Leave `counting` and `dupradar` untouched: an input carries either the RNA fields or the DNA one.

- [ ] **Step 4: Add the DNA citations**

`citations.rs` gains mosdepth (`https://github.com/brentp/mosdepth`) alongside the existing samtools and preseq entries, emitted when the corresponding `dna` tool is enabled.

- [ ] **Step 5: Manual smoke test**

```bash
cargo build --release
./target/release/rustqc dna tests/data/dna/test.dna.bam --outdir /tmp/dna-smoke --window-size 500
find /tmp/dna-smoke -type f | sort
awk '$1=="total"' /tmp/dna-smoke/mosdepth/*.mosdepth.summary.txt
```

Expected: the output tree exists and the total line reads `total 40001 247878 6.20 0 867`. If it differs, Task 8 will name the offending file and column; do not adjust the fixtures.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "feat(dna): wire up the run_dna pipeline

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 8: Parity tests

**Files:**
- Create: `tests/dna_integration_test.rs`

**Interfaces:**
- Consumes: the fixtures from Task 1 and the binary from Task 7.

- [ ] **Step 1: Write the version guard**

```rust
/// The fixtures were generated by the tool versions recorded alongside them.
/// Regenerating with a different version changes the expected output, so the
/// suite refuses to run against a mismatched VERSIONS.txt rather than
/// reporting a parity failure that is really a version skew.
#[test]
fn fixture_tool_versions_are_the_pinned_ones() {
    let versions = std::fs::read_to_string("tests/expected/dna/VERSIONS.txt").unwrap();
    assert!(versions.contains("mosdepth\t0.3.14"), "unexpected mosdepth fixture version: {versions}");
    assert!(versions.contains("samtools\t1.24"), "unexpected samtools fixture version: {versions}");
}
```

- [ ] **Step 2: Write the parity tests**

One test per output file. Each runs the binary into a scratch directory once, shared through a `std::sync::OnceLock`, then compares against the fixture:
- `test.mosdepth.summary.txt`, `.global.dist.txt` and `.region.dist.txt`: parse into `(chrom, key, value)` rows and compare row by row, integers exactly and floats within a relative 1e-6, so a failure names the offending row rather than dumping a whole-file diff.
- `per-base.bed.gz`, `regions.bed.gz` and `thresholds.bed.gz`: decompress both sides and compare interval by interval.
- `test.stats.txt`: compare the `SN` section key by key, then each histogram section row by row. The `rna` suite already proves this writer; the point here is that the DNA pipeline feeds it the same counters.
- `test.flagstat.txt` and `test.idxstats.txt`: exact string equality.

- [ ] **Step 3: Run**

```bash
cargo test --release --test dna_integration_test 2>&1 | tail -30
```

Expected: all green. A failure here is a real defect in the engine; the fixtures come from the upstream tools themselves.

- [ ] **Step 4: Full gate, then commit and open the PR**

```bash
cargo fmt --check && cargo clippy -- -D warnings && cargo test --release
git add -A
git commit -m "test(dna): add parity tests against mosdepth and samtools

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
git push -u fork feat/dna-skeleton
```

Open the PR against `seqeralabs/RustQC:main` from the fork, naming PR #152 as its base in the description since the stack is sequential.

---

## Notes for the executor

- Tasks 4 and 5 are where this PR succeeds or fails. Do not merge them into one commit; the overlap correction must be independently reviewable.
- If a parity test fails in Task 8, resist adjusting the fixture. The fixtures are upstream output. Find the defect.
- Picard and Qualimap outputs are out of scope. If `--targets` is passed in this PR, accept the flag and store it, but emit nothing that depends on it; `CollectHsMetrics` lands in PR4.
