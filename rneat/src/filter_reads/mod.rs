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

pub fn main(
    bed_file: PathBuf,
    fastq1_file: Option<PathBuf>,
    fastq2_file: Option<PathBuf>,
    vcf_file: Option<PathBuf>,
    output_dir: PathBuf,
) -> Result<(), FilterReadsError> {
    info!("values receieved: {:?}, {:?}, {:?}, {:?}, {:?}", bed_file, fastq1_file, fastq2_file, vcf_file, output_dir);
    
    Ok(())
}