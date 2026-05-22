//! `compare-vcfs` subcommand. Compares a downstream variant-caller VCF (called)
//! against a NEAT-simulated truth VCF (golden) and reports true positives,
//! false negatives, and false positives.
//!
//! Phase 1 implements exact-match classification by `(position, ref, alt)` key
//! plus a JSON / TXT rollup. Equivalence detection (denotation-different
//! indels) and NEAT-aware FN attribution are planned for later phases — see
//! issue #127.
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
