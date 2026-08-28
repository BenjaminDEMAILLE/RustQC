# Design: `rustqc protein` subcommand

Date: 2026-08-29
Status: approved for planning
Scope: new protein QC subcommand covering three disjoint domains, delivered as a stack of five PRs

## 1. Goal

Add a `protein` subcommand covering the three things people mean by "protein
QC", each validated against a real upstream tool the way the `rna` and `dna`
pipelines are.

The three are genuinely disjoint. They take different inputs, produce different
outputs, and share no analysis code:

| Mode | Input | Question it answers |
| --- | --- | --- |
| `sequence` | protein FASTA | is this proteome well-formed and complete? |
| `coding` | BAM plus GTF | where do reads fall relative to coding sequence? |
| `spectra` | mzML | is this mass spectrometry run healthy? |

## 2. Command surface

```
rustqc protein sequence <FASTA>... [OPTIONS]
rustqc protein coding   <BAM>... --gtf <GTF> [OPTIONS]
rustqc protein spectra  <MZML>... [OPTIONS]
```

One subcommand with three modes, rather than three subcommands or one flat
command with mutually exclusive flags. The inputs share nothing, so a flat
command would need a validation matrix explaining which flags apply when.
Modes make the incompatibility structural instead.

Shared options keep the names, short flags and `RUSTQC_*` environment variables
they have in `rna` and `dna`: `-o/--outdir`, `--sample-name`, `--flat-output`,
`-c/--config`, `-j/--json-summary`, `-t/--threads`, `-q/--quiet`,
`-v/--verbose`.

Mode-specific options:

| Mode | Option | Default | Purpose |
| --- | --- | --- | --- |
| `sequence` | `--min-length <N>` | 0 | ignore sequences shorter than this |
| `sequence` | `--expect-stop` | off | treat a missing terminal `*` as a defect |
| `coding` | `-g/--gtf <GTF>` | required | annotation, as for `rna` |
| `coding` | `-Q/--mapq <N>` | 30 | mapping quality floor, as for `rna` |
| `coding` | `--riboseq` | off | also compute periodicity and P-site metrics |
| `spectra` | `--ms-levels <N,...>` | `1,2` | levels to report on |
| `spectra` | `--mzqc` | off | also write an mzQC document |

## 3. Module layout

```
src/
  protein/
    mod.rs
    sequence/
      mod.rs          FASTA parsing, per-sequence records
      stats.rs        length distribution, N50, composition
      defects.rs      internal stops, non-standard residues, duplicates
      output.rs       the seqkit-compatible table and the RustQC report
    coding/
      mod.rs
      regions.rs      CDS, UTR, intron and intergenic interval sets from the GTF
      assignment.rs   per-base assignment of aligned bases to a region class
      riboseq.rs      three-nucleotide periodicity, P-site offsets
      output.rs       CollectRnaSeqMetrics-compatible output
    spectra/
      mod.rs
      ingest.rs       mzML reading via mzdata
      metrics.rs      per-run and per-level spectrum metrics
      mzqc.rs         mzQC document writing
      output.rs
```

`sequence` and `coding` add no dependency. `coding` reuses `crate::gtf`, the
`crate::common` modules extracted for the DNA stack, and the per-contig worker
pattern. `spectra` adds `mzdata`.

## 4. The mzdata dependency

`mzdata 0.66.5` reads mzML and builds on Rust 1.87, the project's MSRV. It was
checked before this spec was written: it parses the 3.4 MB `small.mzML` from
its own test suite into 14 MS1 and 34 MS2 spectra holding 305213 peaks.

It pulls 68 transitive packages with `--no-default-features --features mzml`,
which is a real cost for a project that ships a single static binary. It
therefore goes behind a cargo feature:

```toml
[features]
default = ["snmalloc", "proteomics"]
proteomics = ["dep:mzdata"]
```

Building with `--no-default-features` drops the `spectra` mode, which then
reports that it was compiled out rather than failing obscurely.

One caution for whoever updates it: `cargo add mzdata@0.66.5` failed once
against a stale index because the version requires `quick-xml ^0.41`, published
only shortly before. `cargo update` first if resolution fails.

## 5. Parity targets

Every mode is validated against a real tool, as the `rna` and `dna` pipelines
are. This was checked before committing to the design; all three are
obtainable.

**`sequence` against seqkit 2.13.0.** `seqkit stats -a` reports sequence count,
total, minimum, average and maximum length, the length quartiles, N50 and its
index, and correctly detects `type` as `Protein`. RustQC reproduces that table
and adds what seqkit does not report: amino acid composition, internal stop
codons, non-standard residues, and duplicate sequences.

**`coding` against Picard `CollectRnaSeqMetrics`.** It reports `CODING_BASES`,
`UTR_BASES`, `INTRONIC_BASES`, `INTERGENIC_BASES` and their percentages, which
is exactly the per-base region assignment this mode computes. Picard 3.4.0 is
already pinned in `tests/create_dna_test_data.sh`.

**`spectra` against pyteomics.** The reference figures are derived with
`pyteomics.mzml`, installed from PyPI by the generation script. Where the mzQC
document is concerned the reference is the HUPO-PSI schema rather than another
implementation.

Ribo-seq periodicity has no settled reference implementation that installs
cleanly, so it is validated by unit tests over synthetic reads with a known
frame, not by fixtures. That is stated rather than glossed over.

## 6. Test data

All fixtures are small, real and public.

| Mode | File | Source | Size |
| --- | --- | --- | --- |
| `sequence` | `yeast_UPS_mini.fasta` | nf-core/test-datasets | 4 kB |
| `sequence` | `protein_mini_with_cazymes.faa` | nf-core/test-datasets | 7 kB |
| `coding` | `test.bam`, `test.gtf` | already in the repository | 5 kB |
| `spectra` | `peakpicker_tutorial_1.mzML` | nf-core/test-datasets | 1.9 MB |

Roughly 2 MB in total, inside the 10 MB budget the DNA stack set.

`tests/create_protein_test_data.sh` downloads the inputs and regenerates the
reference outputs, pinning seqkit 2.13.0, Picard 3.4.0 and the pyteomics
version, and refusing to run against others, exactly as the DNA script does.

## 7. PR stack

Based on `feat/dna-qualimap-docs`, the tip of the DNA stack, because `coding`
reuses `src/common` from PR #152 and the per-contig worker pattern from
PR #153. If the DNA stack merges first this rebases onto `main` cleanly.

1. **`feat/protein-skeleton`** — `ProteinArgs` with the three modes,
   `ProteinConfig`, `run_protein`, and the `sequence` mode complete with
   seqkit parity.
2. **`feat/protein-coding`** — region interval sets from the GTF, per-base
   assignment, `CollectRnaSeqMetrics` parity.
3. **`feat/protein-riboseq`** — `--riboseq`: three-nucleotide periodicity,
   P-site offset estimation, per-frame read counts.
4. **`feat/protein-spectra`** — the `mzdata` dependency behind the `proteomics`
   feature, mzML ingestion, per-run and per-level metrics, pyteomics parity.
5. **`feat/protein-mzqc-docs`** — mzQC document output, the documentation
   pages, CHANGELOG and a tracking issue.

## 8. Risks

- **Scope.** Three domains under one subcommand is more surface than the DNA
  stack carried. The mode split keeps each PR reviewable, but the series is
  long and any one mode could be dropped without hurting the others.
- **The mzdata dependency.** 68 packages is a large addition for one mode. The
  feature flag contains it, but the default build grows. Worth a maintainer's
  opinion before PR4 lands.
- **Ribo-seq has no clean reference.** Unit tests over synthetic reads are
  weaker evidence than fixture parity. The metrics are standard enough that
  this is acceptable, but it is the weakest link in the series.
- **`coding` overlaps `rna`.** `CollectRnaSeqMetrics` is an RNA-seq tool, and
  someone will reasonably ask why this is not in `rustqc rna`. The answer is
  that the question being asked is about coding sequence rather than about
  transcripts, but the boundary is genuinely arguable and should be settled in
  review rather than after.

## 9. Open items

None blocking. The three technical unknowns were resolved before writing:
mzdata builds on the MSRV and reads mzML, seqkit is installed and reports
protein FASTA statistics, and small public fixtures exist for every mode.
