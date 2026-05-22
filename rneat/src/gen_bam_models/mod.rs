pub mod errors;
pub mod utils;

use std::path::PathBuf;
use crate::gen_bam_models::errors::GenBamModelsError;

pub fn main(config_file: &PathBuf) -> Result<(), GenBamModelsError> {
    utils::runner::runner(config_file)
}
