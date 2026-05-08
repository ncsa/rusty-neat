pub mod errors;
pub mod utils;

use std::path::PathBuf;
use crate::gen_gc_bias_model::{
    errors::GenGcBiasModelError,
    utils::runner::runner,
};

pub fn main(config_file: &PathBuf) -> Result<(), GenGcBiasModelError> {
    runner(config_file)
}