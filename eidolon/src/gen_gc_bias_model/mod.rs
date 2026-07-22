pub mod errors;
pub mod utils;

use crate::gen_gc_bias_model::{errors::GenGcBiasModelError, utils::runner::runner};
use std::path::PathBuf;

pub fn main(config_file: &PathBuf) -> Result<(), GenGcBiasModelError> {
    runner(config_file)
}
