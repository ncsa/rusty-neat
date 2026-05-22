pub mod errors;
pub mod utils;

use crate::gen_bam_models::errors::GenBamModelsError;
use std::path::PathBuf;

pub fn main(config_file: &PathBuf) -> Result<(), GenBamModelsError> {
    utils::runner::runner(config_file)
}
