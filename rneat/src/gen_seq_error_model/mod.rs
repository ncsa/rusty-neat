pub mod errors;
pub mod utils;

use std::path::PathBuf;
use crate::gen_seq_error_model::{
    errors::GenSeqErrorModelError,
    utils::config::RunConfiguration,
};

pub fn main(config_file: &PathBuf) -> Result<(), GenSeqErrorModelError> {
    let config = RunConfiguration::from(config_file)?;
    utils::runner::runner(&config)
}
