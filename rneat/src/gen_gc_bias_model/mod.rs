
/// This file will create a mutation model based on input data.
///
pub mod errors;
pub mod utils;

use log::*;
use std::{
    collections::HashMap, 
    path::PathBuf
};
use common::{
    file_tools::fasta_reader::read_fasta, 
    structs::variants::Variant,
};
use crate::gen_gc_bias_model::{
    errors::GenGcBiasModelError,
    utils::config::RunConfiguration,
    utils::runner::runner,
};

pub fn main(config_file: &PathBuf) -> Result<(), GenGcBiasModelError> {
    Ok(())
}