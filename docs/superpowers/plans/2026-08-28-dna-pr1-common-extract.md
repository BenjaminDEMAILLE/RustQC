# DNA PR1: extract shared modules into `src/common/`

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the non-RNA-specific analysis code out of `src/rna/` into a new `src/common/` module so the forthcoming `dna` subcommand can consume it, with zero change in observable behaviour.

**Architecture:** Pure refactor. Six existing files move wholesale, one struct (`BamStatAccum`) is lifted out of `src/rna/rseqc/accumulators.rs` into its own module, and `src/rna/` gains re-export shims so every existing `crate::rna::...` path and the published 0.2.x library surface keep resolving. Correctness is proven by a golden-output comparison of the `rustqc rna` binary before and after, not by new unit tests.

**Tech Stack:** Rust 2021, rust-htslib, anyhow, rayon, cargo test / clippy / fmt.

**Spec:** `docs/superpowers/specs/2026-08-28-dna-command-design.md` (section 3, Module layout)

## Global Constraints

- Rust edition 2021, `rust-version = "1.87"`. Do not raise either.
- Default `rustfmt`. No `rustfmt.toml` exists and none may be created.
- `cargo clippy -- -D warnings` must pass. Targeted `#[allow(...)]` needs a justifying comment.
- `anyhow::Result<T>` for all fallible functions. No custom error types.
- `unwrap()` / `expect()` in test code only.
- Every source file opens with a `//!` module doc comment. Every public item gets a `///` doc comment.
- This PR changes **no behaviour**. Any output byte that differs is a bug in the refactor.
- Do not reformat, rename, or "improve" moved code. Moves are verbatim except for import paths.

---

### Task 1: Capture the golden-output regression net

A pure refactor has no new tests to write. The safety net is a byte-for-byte comparison of the binary's output before and after the move.

**Files:**
- Create: `/tmp/rustqc-golden/before/` (scratch, not committed)

**Interfaces:**
- Produces: a `before/` output tree that Task 6 compares against.

- [ ] **Step 1: Confirm the tree is clean and on the right branch**

```bash
git status --porcelain   # expect empty
git rev-parse --abbrev-ref HEAD
```

Expected: no output from the first command. If the branch is not `feat/dna-common-extract`, create it from `feat/dna-command-spec`:

```bash
git switch -c feat/dna-common-extract
```

- [ ] **Step 2: Build the release binary**

```bash
cargo build --release
```

Expected: success. This is the slow step (static htslib plus the C++ preseq shim); allow several minutes.

- [ ] **Step 3: Run the full rna pipeline on the committed test data**

```bash
rm -rf /tmp/rustqc-golden && mkdir -p /tmp/rustqc-golden/before
./target/release/rustqc rna tests/data/test.bam \
  --gtf tests/data/test.gtf \
  --outdir /tmp/rustqc-golden/before \
  --json-summary /tmp/rustqc-golden/before/summary.json \
  --quiet
```

Expected: exit 0 and a populated output tree.

- [ ] **Step 4: Record a stable manifest**

The JSON summary and `CITATIONS.md` embed a timestamp, a runtime and absolute paths, so they are excluded from the checksum set.

```bash
cd /tmp/rustqc-golden/before && \
  find . -type f ! -name 'summary.json' ! -name 'CITATIONS.md' -print0 \
  | sort -z | xargs -0 shasum -a 256 > /tmp/rustqc-golden/before.sha256
wc -l /tmp/rustqc-golden/before.sha256
```

Expected: a non-empty manifest listing every output file. Note the line count; Task 6 asserts the same count.

- [ ] **Step 5: Record the baseline test result**

```bash
cargo test --release 2>&1 | tail -20 > /tmp/rustqc-golden/before-tests.txt
cat /tmp/rustqc-golden/before-tests.txt
```

Expected: all tests pass. Any test that is already failing on `main` must be noted now so it is not blamed on this refactor.

No commit for this task; nothing in the repository changed.

---

### Task 2: Create `src/common/` and move the three leaf modules

`bam_flags.rs`, `cpp_rng.rs` and `preseq.rs` sit directly under `src/rna/` and have no RNA-specific dependency. `preseq.rs` depends on `cpp_rng`, so the three move together.

**Files:**
- Create: `src/common/mod.rs`
- Move: `src/rna/bam_flags.rs` to `src/common/bam_flags.rs`
- Move: `src/rna/cpp_rng.rs` to `src/common/cpp_rng.rs`
- Move: `src/rna/preseq.rs` to `src/common/preseq.rs`
- Modify: `src/rna/mod.rs`, `src/lib.rs`

**Interfaces:**
- Produces: `crate::common::bam_flags`, `crate::common::cpp_rng`, `crate::common::preseq`, all with their existing public items unchanged.
- Produces: shims `crate::rna::bam_flags`, `crate::rna::cpp_rng`, `crate::rna::preseq` that re-export the above.

- [ ] **Step 1: Move the files with git so history follows**

```bash
mkdir -p src/common
git mv src/rna/bam_flags.rs src/common/bam_flags.rs
git mv src/rna/cpp_rng.rs   src/common/cpp_rng.rs
git mv src/rna/preseq.rs    src/common/preseq.rs
```

- [ ] **Step 2: Write `src/common/mod.rs`**

```rust
//! Analysis modules shared between the `rna` and `dna` pipelines.
//!
//! Nothing in this module is specific to a library preparation or an assay:
//! BAM flag helpers, the C++ RNG shim used for preseq bootstrap
//! reproducibility, the preseq `lc_extrap` implementation, read-level
//! alignment statistics, and the samtools-compatible output writers.

pub mod bam_flags;
pub mod cpp_rng;
pub mod preseq;
```

- [ ] **Step 3: Declare the module in `src/lib.rs`**

Add `pub mod common;` to the module list, keeping the existing alphabetical order (it goes first, before `pub mod config;`).

- [ ] **Step 4: Replace the moved declarations in `src/rna/mod.rs` with shims**

Delete the `pub mod bam_flags;`, `pub mod cpp_rng;` and `pub mod preseq;` lines and put this in their place. A `pub mod` and a `pub use` of the same name cannot coexist, so the old lines must go.

```rust
// These analyses are not RNA-specific and now live in `crate::common`.
// Re-exported here so existing `crate::rna::...` paths and the published
// 0.2.x library surface keep resolving. Drop the shims at 1.0.
pub use crate::common::{bam_flags, cpp_rng, preseq};
```

- [ ] **Step 5: Fix the one internal path inside the moved preseq module**

`src/common/preseq.rs` line 14 reads `use super::cpp_rng::CppMt19937;`. Under `src/common/` that `super::` now resolves to `crate::common`, which is where `cpp_rng` lives, so the line is already correct. Verify it rather than editing it:

```bash
grep -n "cpp_rng" src/common/preseq.rs
```

Expected: `14:use super::cpp_rng::CppMt19937;` and nothing else.

- [ ] **Step 6: Compile and let the compiler find the remaining call sites**

```bash
cargo check 2>&1 | grep -E "^(error|warning: unused)" | head -30
```

Expected: errors only from `crate::rna::bam_flags::*` glob imports if the shim is wrong. A `pub use` re-export supports glob imports, so `use crate::rna::bam_flags::*;` in `dupradar/counting.rs`, `qualimap/accumulator.rs` and `rseqc/accumulators.rs` should keep working untouched. If any error remains, fix the import at the call site to `crate::common::bam_flags::*` rather than weakening the shim.

- [ ] **Step 7: Full build and test**

```bash
cargo build --release && cargo test --release 2>&1 | tail -20
```

Expected: same result as `/tmp/rustqc-golden/before-tests.txt`.

- [ ] **Step 8: Commit**

```bash
git add -A
git commit -m "refactor: move bam_flags, cpp_rng and preseq into src/common/

These three modules carry no RNA-specific logic and are needed by the
forthcoming dna subcommand. src/rna re-exports them so every existing
crate::rna::... path and the published 0.2.x library surface keep working.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Move `bam_stat` and the samtools trio into `src/common/`

`stats.rs`, `flagstat.rs` and `idxstats.rs` each `use super::bam_stat::BamStatResult`, so `bam_stat.rs` moves with them. RSeQC's `bam_stat` is read-level only and needs no annotation, so it belongs in `common`.

**Files:**
- Move: `src/rna/rseqc/bam_stat.rs` to `src/common/bam_stat.rs`
- Move: `src/rna/rseqc/stats.rs` to `src/common/samtools/stats.rs`
- Move: `src/rna/rseqc/flagstat.rs` to `src/common/samtools/flagstat.rs`
- Move: `src/rna/rseqc/idxstats.rs` to `src/common/samtools/idxstats.rs`
- Create: `src/common/samtools/mod.rs`
- Modify: `src/common/mod.rs`, `src/rna/rseqc/mod.rs`, `src/main.rs:1565,1575,1585`

**Interfaces:**
- Consumes: `crate::common` from Task 2.
- Produces: `crate::common::bam_stat::{BamStatResult, GcDepthBin}` and `crate::common::samtools::{stats, flagstat, idxstats}` with `write_stats`, `write_flagstat`, `write_idxstats` unchanged.
- Produces: shims `crate::rna::rseqc::{bam_stat, stats, flagstat, idxstats}`.

- [ ] **Step 1: Move the files**

```bash
mkdir -p src/common/samtools
git mv src/rna/rseqc/bam_stat.rs  src/common/bam_stat.rs
git mv src/rna/rseqc/stats.rs     src/common/samtools/stats.rs
git mv src/rna/rseqc/flagstat.rs  src/common/samtools/flagstat.rs
git mv src/rna/rseqc/idxstats.rs  src/common/samtools/idxstats.rs
```

- [ ] **Step 2: Write `src/common/samtools/mod.rs`**

```rust
//! samtools-compatible output writers.
//!
//! Reproduce the exact output formats of `samtools stats`, `samtools flagstat`
//! and `samtools idxstats` from the counters gathered in
//! [`crate::common::bam_stat::BamStatResult`], so that MultiQC and
//! `plot-bamstats` parse RustQC output as if samtools had produced it.

pub mod flagstat;
pub mod idxstats;
pub mod stats;
```

- [ ] **Step 3: Add the new modules to `src/common/mod.rs`**

The module list becomes, in alphabetical order:

```rust
pub mod bam_flags;
pub mod bam_stat;
pub mod cpp_rng;
pub mod preseq;
pub mod samtools;
```

- [ ] **Step 4: Fix the `super::bam_stat` imports in the three moved writers**

Each of the three files imports the result type from its old sibling. `src/common/samtools/stats.rs` line 13 and `src/common/samtools/flagstat.rs` line 11 read `use super::bam_stat::...`; under `src/common/samtools/` that `super::` now resolves to `crate::common::samtools`, which has no `bam_stat`. Point them at the parent module instead:

```bash
sed -i '' 's|^use super::bam_stat::|use crate::common::bam_stat::|' \
  src/common/samtools/stats.rs \
  src/common/samtools/flagstat.rs \
  src/common/samtools/idxstats.rs
grep -n "bam_stat" src/common/samtools/*.rs
```

Expected: every hit now reads `use crate::common::bam_stat::...`.

- [ ] **Step 5: Replace the moved declarations in `src/rna/rseqc/mod.rs` with shims**

Delete the `pub mod bam_stat;`, `pub mod flagstat;`, `pub mod idxstats;` and `pub mod stats;` lines, and add:

```rust
// bam_stat and the samtools writers are read-level and assay-agnostic; they
// now live in `crate::common`. Re-exported so existing
// `crate::rna::rseqc::...` paths and the published 0.2.x library surface
// keep resolving. Drop the shims at 1.0.
pub use crate::common::bam_stat;
pub use crate::common::samtools::{flagstat, idxstats, stats};
```

- [ ] **Step 6: Update the three call sites in `src/main.rs`**

Lines 1565, 1575 and 1585 currently call through the `rna::rseqc` path. Point them at the real location so the binary does not depend on a compatibility shim:

```bash
sed -i '' \
  -e 's|rna::rseqc::flagstat::write_flagstat|common::samtools::flagstat::write_flagstat|' \
  -e 's|rna::rseqc::idxstats::write_idxstats|common::samtools::idxstats::write_idxstats|' \
  -e 's|rna::rseqc::stats::write_stats|common::samtools::stats::write_stats|' \
  src/main.rs
grep -n "write_flagstat\|write_idxstats\|write_stats" src/main.rs
```

Then add `common` to the `use rustqc::{...}` list near the top of `src/main.rs` (it currently reads `use rustqc::{config, cpu, gtf, rna, summary};`).

- [ ] **Step 7: Compile and resolve whatever the compiler reports**

```bash
cargo check 2>&1 | grep -E "^error" -A3 | head -40
```

Expected: clean. The likely residue is `use super::bam_stat::{BamStatResult, GcDepthBin};` inside `src/rna/rseqc/accumulators.rs`, which still resolves through the shim and needs no change yet; Task 4 rewrites it.

- [ ] **Step 8: Build and test**

```bash
cargo build --release && cargo test --release 2>&1 | tail -20
```

Expected: same result as the baseline.

- [ ] **Step 9: Commit**

```bash
git add -A
git commit -m "refactor: move bam_stat and the samtools writers into src/common/

bam_stat is read-level and needs no annotation, and the samtools stats,
flagstat and idxstats writers consume its result type, so all four move
together into src/common/. src/rna/rseqc re-exports them.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Lift `BamStatAccum` out of the RSeQC accumulators

`BamStatAccum` is the struct that gathers every read-level counter feeding bam_stat and the samtools trio. It is already self-contained: `process_read(&mut self, record, mapq_cut)` takes no annotation, and `merge` and `into_result` touch nothing RNA-specific. It is merely *stored* in `src/rna/rseqc/accumulators.rs` alongside the RNA accumulators. The `dna` pipeline needs it, so it moves to its own module.

**Files:**
- Create: `src/common/bam_stat_accum.rs`
- Modify: `src/rna/rseqc/accumulators.rs` (remove lines 118-1315 and 2346-2434, add an import)
- Modify: `src/common/mod.rs`

**Interfaces:**
- Consumes: `crate::common::bam_stat::{BamStatResult, GcDepthBin}` from Task 3, `crate::common::bam_flags::*` from Task 2.
- Produces: `crate::common::bam_stat_accum::BamStatAccum` with the unchanged public API `process_read(&mut self, record: &bam::Record, mapq_cut: u8)`, `flush_cov_buf_all(&mut self)`, `merge(&mut self, other: BamStatAccum)`, `into_result(self) -> BamStatResult`, and `Default`.
- Produces: shim `crate::rna::rseqc::accumulators::BamStatAccum`.

- [ ] **Step 1: Identify the exact blocks to move**

Two contiguous regions of `src/rna/rseqc/accumulators.rs`, verified before cutting:

```bash
sed -n '118,124p'   src/rna/rseqc/accumulators.rs   # doc comment + struct opening
sed -n '1313,1318p' src/rna/rseqc/accumulators.rs   # end of impl, start of InferExpAccum
sed -n '2346,2352p' src/rna/rseqc/accumulators.rs   # section divider + impl BamStatAccum
sed -n '2431,2437p' src/rna/rseqc/accumulators.rs   # end of into_result, start of InferExpAccum
```

Expected: block A is lines 118 through 1315 (doc comment, `pub struct BamStatAccum`, and the `impl BamStatAccum` holding `process_read`, `flush_cov_buf_all` and `merge`). Block B is the `impl BamStatAccum` holding `into_result`, ending just before `impl InferExpAccum`. If the printed boundaries do not match this description, adjust the line numbers to the real ones before proceeding; the file may have shifted.

- [ ] **Step 2: Extract both blocks into the new file**

```bash
{
  cat <<'EOF'
//! Read-level alignment statistics accumulator.
//!
//! [`BamStatAccum`] gathers, in a single pass over the records, every counter
//! consumed by RSeQC `bam_stat` and by the samtools-compatible `stats`,
//! `flagstat` and `idxstats` writers. It needs no annotation and no library
//! protocol, so both the `rna` and `dna` pipelines drive the same struct: each
//! parallel worker owns one, and they are merged before conversion.

use std::collections::HashMap;

use rust_htslib::bam;

use crate::common::bam_flags::*;
use crate::common::bam_stat::{BamStatResult, GcDepthBin};

/// Default GC-depth bin size in base pairs (matches upstream samtools default).
const GCD_BIN_SIZE: u64 = 20_000;

EOF
  sed -n '118,1315p' src/rna/rseqc/accumulators.rs
  echo
  sed -n '2346,2434p' src/rna/rseqc/accumulators.rs
} > src/common/bam_stat_accum.rs
wc -l src/common/bam_stat_accum.rs
```

Expected: roughly 1300 lines.

- [ ] **Step 3: Delete the moved blocks from the original file**

Delete the higher range first so the lower line numbers stay valid.

```bash
sed -i '' '2346,2434d' src/rna/rseqc/accumulators.rs
sed -i '' '118,1315d'  src/rna/rseqc/accumulators.rs
```

- [ ] **Step 4: Re-export the type from its old location**

`RseqcAccumulators` holds a `BamStatAccum` field and `main.rs` names the type, so the old path must keep working. Add to the import block at the top of `src/rna/rseqc/accumulators.rs`:

```rust
// BamStatAccum is read-level and assay-agnostic; it lives in `crate::common`
// and is shared with the dna pipeline. Re-exported so existing paths resolve.
pub use crate::common::bam_stat_accum::BamStatAccum;
```

Also delete the now-orphaned `const GCD_BIN_SIZE` near line 17 and the `use super::bam_stat::{BamStatResult, GcDepthBin};` import if the compiler reports them unused. Leave every other import alone until the compiler asks.

- [ ] **Step 5: Register the module**

Add `pub mod bam_stat_accum;` to `src/common/mod.rs`, after `pub mod bam_stat;`.

- [ ] **Step 6: Compile and iterate on the errors**

```bash
cargo check 2>&1 | grep -E "^error" -A5 | head -60
```

Expected failure modes and their fixes:
- *unresolved import / unused import in `accumulators.rs`*: delete the import; it belonged to the moved code.
- *private item accessed from `bam_stat_accum.rs`*: a helper function that `BamStatAccum` calls stayed behind. Move that helper into `src/common/bam_stat_accum.rs` too, and delete it from `accumulators.rs`. The array-merging helper just above the old line 2346 is the likely candidate.
- *`GCD_BIN_SIZE` unused in `accumulators.rs`*: delete it there; the new file declares its own.

Repeat until `cargo check` is clean.

- [ ] **Step 7: Format, lint, test**

```bash
cargo fmt
cargo fmt --check && cargo clippy -- -D warnings && cargo test --release 2>&1 | tail -20
```

Expected: all clean, tests matching the baseline. Note that `cargo fmt` may reindent the moved blocks; that is the one formatting change this PR is allowed to make.

- [ ] **Step 8: Commit**

```bash
git add -A
git commit -m "refactor: lift BamStatAccum into src/common/bam_stat_accum.rs

BamStatAccum gathers the read-level counters behind bam_stat and the
samtools writers. Its process_read takes only a record and a MAPQ cutoff,
so it is assay-agnostic and the dna pipeline will drive the same struct.
rna::rseqc::accumulators re-exports it.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Update the documentation that describes the module tree

**Files:**
- Modify: `src/lib.rs` (module list in the crate doc comment, around lines 20-30)
- Modify: `AGENTS.md` (the Project Structure tree, and the sentence stating there is no `lib.rs`)
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: the final layout from Tasks 2 through 4.
- Produces: nothing code-facing.

- [ ] **Step 1: Update the crate doc module list in `src/lib.rs`**

Add a `common` bullet above the `rna` bullet:

```rust
//! - [`common`] — analyses shared by every pipeline: BAM flag helpers,
//!   read-level statistics ([`common::bam_stat`], [`common::bam_stat_accum`]),
//!   the samtools-compatible writers ([`common::samtools`]), and preseq
//!   library complexity extrapolation ([`common::preseq`]).
```

And amend the `rna` bullet so it no longer claims `rna::preseq` as its own; it now reads:

```rust
//! - [`rna`] — the RNA-Seq analysis modules:
//!   - [`rna::dupradar`], [`rna::featurecounts`], [`rna::qualimap`], [`rna::rseqc`].
```

- [ ] **Step 2: Update the tree in `AGENTS.md`**

Insert the `common/` block before `rna/` and delete the moved entries from under `rna/`:

```
  common/
    mod.rs            — Re-exports the shared modules
    bam_flags.rs      — BAM flag constants and aux-tag helpers
    bam_stat.rs       — bam_stat.py reimplementation, result types
    bam_stat_accum.rs — Read-level counter accumulator feeding bam_stat and samtools
    cpp_rng.rs        — C++ RNG FFI shim for preseq bootstrap reproducibility
    preseq.rs         — preseq lc_extrap library complexity extrapolation
    samtools/
      mod.rs          — Re-exports the samtools writers
      stats.rs        — samtools stats full output (SN + all histogram sections)
      flagstat.rs     — samtools flagstat-compatible output
      idxstats.rs     — samtools idxstats-compatible output
```

- [ ] **Step 3: Correct the stale sentence in `AGENTS.md`**

It currently claims "top-level modules (`cli`, `config`, `io`, `gtf`, `rna`) declared in `main.rs`, no `lib.rs`". A `lib.rs` has existed since #101. Replace with:

```
Nested module structure. The library crate root is `src/lib.rs`, which declares
`common`, `config`, `cpu`, `gtf`, `io`, `rna` and `summary`; the binary
(`src/main.rs`) additionally declares `cli`, `citations` and `ui`.
Inter-module access uses `crate::` paths (e.g. `use crate::common::bam_stat_accum::BamStatAccum;`).
Assay-agnostic analyses belong in `common`; put new code under `rna` only if it
needs a gene annotation or a library strand protocol.
```

- [ ] **Step 4: Add the CHANGELOG entry**

Under the unreleased heading, in the existing style of the file:

```markdown
### Changed

- Internal: assay-agnostic analyses (BAM flag helpers, read-level statistics,
  the samtools stats/flagstat/idxstats writers, preseq) moved from `rna` to a
  new `common` module. The old `rustqc::rna::...` paths still resolve through
  re-exports, so this is not a breaking change for library users.
```

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs: describe the src/common module split

Also corrects the AGENTS.md claim that the crate has no lib.rs, which has
been untrue since #101.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: Prove the refactor changed nothing, then open the PR

**Files:** none modified.

**Interfaces:**
- Consumes: the golden manifest from Task 1.

- [ ] **Step 1: Re-run the pipeline on the same input**

```bash
mkdir -p /tmp/rustqc-golden/after
cargo build --release
./target/release/rustqc rna tests/data/test.bam \
  --gtf tests/data/test.gtf \
  --outdir /tmp/rustqc-golden/after \
  --json-summary /tmp/rustqc-golden/after/summary.json \
  --quiet
```

- [ ] **Step 2: Compare against the golden manifest**

```bash
cd /tmp/rustqc-golden/after && \
  find . -type f ! -name 'summary.json' ! -name 'CITATIONS.md' -print0 \
  | sort -z | xargs -0 shasum -a 256 > /tmp/rustqc-golden/after.sha256
diff /tmp/rustqc-golden/before.sha256 /tmp/rustqc-golden/after.sha256 && echo "IDENTICAL"
```

Expected: `IDENTICAL`. Any difference is a refactor bug and must be fixed before the PR is opened, not explained away. Plots are the one legitimate source of noise if they embed a timestamp; check that hypothesis by diffing the offending file, and only then exclude it, naming it explicitly in the PR description.

- [ ] **Step 3: Compare the JSON summary, ignoring the volatile fields**

```bash
for f in before after; do
  python3 -c "import json,sys; d=json.load(open('/tmp/rustqc-golden/$f/summary.json')); \
    [d.pop(k,None) for k in ('timestamp_start','timestamp_end','runtime_seconds')]; \
    [i.pop('runtime_seconds',None) for i in d.get('inputs',[])]; \
    [o.update(path=o['path'].split('/')[-1]) for i in d.get('inputs',[]) for o in i.get('outputs',[])]; \
    print(json.dumps(d,sort_keys=True,indent=1))" > /tmp/rustqc-golden/$f.json
done
diff /tmp/rustqc-golden/before.json /tmp/rustqc-golden/after.json && echo "SUMMARY IDENTICAL"
```

Expected: `SUMMARY IDENTICAL`.

- [ ] **Step 4: Full gate**

```bash
cargo fmt --check && cargo clippy -- -D warnings && cargo test --release
```

Expected: all pass, with the same test count as the baseline. A *lower* count means a test was lost in the move; find it.

- [ ] **Step 5: Confirm the shims actually shim**

The re-exports are load-bearing for library users, and nothing in the binary exercises them any more after Task 3 step 6. Add a compile-time guard as the last unit test in `src/rna/rseqc/mod.rs`:

```rust
#[cfg(test)]
mod compat_tests {
    //! Guards the re-export shims that keep the published 0.2.x paths alive.
    //! These are compile-time assertions; there is nothing to observe at runtime.

    #[test]
    fn moved_modules_are_still_reachable_from_their_old_paths() {
        let _: fn(&crate::rna::rseqc::bam_stat::BamStatResult, &std::path::Path) -> anyhow::Result<()> =
            crate::rna::rseqc::flagstat::write_flagstat;
        let _ = crate::rna::rseqc::accumulators::BamStatAccum::default();
        let _ = crate::rna::bam_flags::BAM_FDUP;
    }
}
```

Both signatures were verified against the tree at `bea5571`: `write_flagstat(&BamStatResult, &Path) -> Result<()>` and `impl Default for BamStatAccum` (declared inside the block that Task 4 moves). If either fails to compile, the move dropped something.

```bash
cargo test --release compat_tests
```

Expected: PASS.

- [ ] **Step 6: Commit and push**

```bash
git add -A
git commit -m "test: guard the rna re-export shims against silent breakage

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
git push -u origin feat/dna-common-extract
```

- [ ] **Step 7: Open the PR**

```bash
gh pr create --base main --title "refactor: extract assay-agnostic analyses into src/common/" --body "$(cat <<'BODY'
Groundwork for the forthcoming `rustqc dna` subcommand. Pure refactor: no
behaviour change, no new dependency, no public API break.

## What moves

| From | To |
| --- | --- |
| `src/rna/bam_flags.rs` | `src/common/bam_flags.rs` |
| `src/rna/cpp_rng.rs` | `src/common/cpp_rng.rs` |
| `src/rna/preseq.rs` | `src/common/preseq.rs` |
| `src/rna/rseqc/bam_stat.rs` | `src/common/bam_stat.rs` |
| `src/rna/rseqc/{stats,flagstat,idxstats}.rs` | `src/common/samtools/` |
| `BamStatAccum` (from `src/rna/rseqc/accumulators.rs`) | `src/common/bam_stat_accum.rs` |

None of these need a gene annotation or a library strand protocol, and the
`dna` pipeline needs all of them. `src/rna` re-exports every moved item, so
`rustqc::rna::preseq`, `rustqc::rna::rseqc::stats` and friends keep resolving
for library users on 0.2.x. A `compat_tests` module pins those paths.

## How it is verified

Beyond the existing suite: `rustqc rna` was run on `tests/data/` before and
after the refactor and every output file compared by SHA-256, plus the JSON
summary compared field by field with timestamps, runtimes and absolute paths
normalised away. Output is identical.

## Design

`docs/superpowers/specs/2026-08-28-dna-command-design.md`, section 3.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
BODY
)"
```

---

## Notes for the executor

- Tasks 2, 3 and 4 each end in a green build and their own commit. If a task cannot be made green, stop and report rather than pressing on; a half-moved module compiles in confusing ways.
- The line numbers in Task 4 are read off the file as committed at `bea5571`. Verify them with the `sed -n` probes in step 1 before cutting. If they have drifted, find the real boundaries by the landmarks named in the step (the `BamStatAccum` doc comment, the `InferExpAccum` doc comment, the "Converter methods" divider) instead of trusting the numbers.
- Do not add `dna` code in this PR. The `src/dna/` tree lands in PR2.
