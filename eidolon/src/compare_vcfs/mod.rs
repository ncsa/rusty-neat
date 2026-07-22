//! `compare-vcfs` subcommand. Compares a downstream variant-caller VCF (called)
//! against a NEAT-simulated truth VCF (golden) and reports true positives,
//! false negatives, and false positives.
//!
//! Pipeline: exact-match classification by `(position, ref, alt)` key →
//! optional ±N bp equivalence sweep (catches denotation-different alternates
//! like left- vs. right-aligned indels) → NEAT-aware FN attribution. Outputs
//! `comparison_summary.{json,txt}`, `FN_with_reasons.vcf`, and (optionally)
//! `FP.vcf`.
pub mod errors;
pub mod utils;

use crate::compare_vcfs::{
    errors::CompareVcfsError, utils::config::RunConfiguration, utils::runner::runner,
};
use std::path::PathBuf;

pub fn main(config_file: &PathBuf) -> Result<(), CompareVcfsError> {
    let run_config = RunConfiguration::from(config_file)?;
    runner(&run_config)?;
    Ok(())
}
