//! RSeQC reimplementations.
//!
//! Standalone reimplementations of selected RSeQC Python quality control scripts
//! for RNA-Seq data analysis.

pub mod accumulators;
pub mod common;
pub mod plots;

pub mod infer_experiment;
pub mod inner_distance;
pub mod junction_annotation;
pub mod junction_saturation;
pub mod read_distribution;
pub mod read_duplication;
pub mod tin;

// bam_stat and the samtools writers are read-level and assay-agnostic; they
// now live in `crate::common`. Re-exported so existing
// `crate::rna::rseqc::...` paths and the published 0.2.x library surface
// keep resolving. Drop the shims at 1.0.
pub use crate::common::bam_stat;
pub use crate::common::samtools::{flagstat, idxstats, stats};
