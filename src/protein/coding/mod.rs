//! Coding-region quality control from an alignment and an annotation.
//!
//! Answers where reads fall relative to protein-coding sequence, which is the
//! question `rna` does not: that pipeline asks about transcripts, this one
//! about the coding part of them.
//!
//! Reproduces Picard `CollectRnaSeqMetrics`'s base assignment.

pub mod output;
pub mod regions;

pub use regions::{Region, RegionSets};
