/// This file will create a mutation model based on input data. 
///
pub mod errors;
pub mod utils;

use std::path::PathBuf;
use common::models::mutation_model;
use crate::create_mutation_model::errors::CreateMutationModelError;

pub fn main(config_file: &PathBuf) -> Result<(), CreateMutationModelError> {
    
    Ok(())
}