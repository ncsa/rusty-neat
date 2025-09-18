pub mod utils;
use thiserror::Error;
use std::num::ParseIntError;
use std::path::PathBuf;
use log::*;

#[derive(Error, Debug)]
pub enum FilterReadsError {
    #[error("Error while filtering the vcf file {0}")]
    VcfFilterError(String),
    #[error("Error while filtering the fastq file {0}")]
    FastqFilterError(String),
    #[error("File does not exist: {0}")]
    FileNotFound(String),
    #[error("Unable to parse file name for extension {0}")]
    MalformedFileName(String),
    #[error("I/O error during filtering of files")]
    IoError(#[from] std::io::Error),
    #[error("Error parsing coordinates from read name: {0}")]
    CoordParseError(#[from] ParseIntError)
}

fn main(configuration: &PathBuf) -> Result<(), FilterReadsError> {
    info!("Configuration file receieved: {:?}", configuration);
    Ok(())
}