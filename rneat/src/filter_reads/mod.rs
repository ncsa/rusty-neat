pub mod utils;
use thiserror::Error;
use std::num::ParseIntError;
use std::path::PathBuf;
use log::*;
use common::{file_tools::bed_reader::{self, BedReaderError}};
use utils::filter_lib::{filter_fastq, filter_vcf};

use crate::filter_reads::utils::config::RunConfiguration;

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
    CoordParseError(#[from] ParseIntError),
    #[error("The bed reader returned an error: {0}")]
    BedReaderErr(#[from] BedReaderError),
    #[error("Unknown char while parsing suspected read name")]
    UnknownChar,
    #[error("Cannot overwrite existing file: {0}")]
    OverwriteFileError(String),
    #[error("Invalid CLI inputs: {0}")]
    CliError(String),
    #[error("Invalid Configuration inputs: {0}")]
    ConfigurationError(String),
}

pub fn main(config: &PathBuf) -> Result<(), FilterReadsError> {
    info!("////////////// Welcome to rusty-neat filter reads! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\");
    info!("Using Configuration file input: {:?}", &config);
    let run_config = RunConfiguration::from(config);
    // bed_file: path to bed file to use for filtering
    // file_to_filter: fastq or vcf ready for filtering, can be gzipped or not
    // is_gzip: whether the input file is gzipped or not
    // is_fastq: wether the
    // output_file: will be gzipped.
    let bed_table = bed_reader::read_bed(&run_config.bed_file)?;
    info!("Running filter-reads on input file.");
    for input_file in run_config.file_map.keys() {
        let (output_file, is_gzip, is_fastq) = run_config.file_map[input_file].clone();
        if is_fastq {
            info!("Filtering input fastq file: {:?}", &input_file);
            filter_fastq(&bed_table, &input_file, is_gzip, &output_file)?;
        } else {
            info!("Filtering input vcf file: {:?}", &input_file);
            filter_vcf(&bed_table, &input_file, is_gzip, &output_file)?;
        }
        info!("Successfully filtered input file {:?}, written to {:?}", input_file, output_file);
    }
    Ok(())
}