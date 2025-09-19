pub mod utils;
use thiserror::Error;
use std::num::ParseIntError;
use std::path::PathBuf;
use log::*;
use common::{file_tools::bed_reader::{self, BedReaderError}};
use utils::filter_lib::{filter_fastq, filter_vcf};

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
}

pub fn main(
    bed_file: PathBuf,
    file_to_filter: PathBuf,
    is_gzip: bool,
    is_fastq: bool,
    output_file: PathBuf,
) -> Result<(), FilterReadsError> {
    // bed_file: path to bed file to use for filtering
    // file_to_filter: fastq or vcf ready for filtering, can be gzipped or not
    // is_gzip: whether the input file is gzipped or not
    // is_fastq: wether the
    // output_file: will be gzipped.
    info!("filter-reads values receieved: {:?}, {:?}, {:?}", bed_file, file_to_filter, output_file);
    let bed_table = bed_reader::read_bed(&bed_file)?;
    info!("Running filter-reads on input file.");
    if is_fastq {
        info!("Filtering input fastq file: {:?}", file_to_filter);
        filter_fastq(&bed_table, &file_to_filter, is_gzip, &output_file)?;
    } else {
        info!("Filtering input vcf file: {:?}", file_to_filter);
        filter_vcf(&bed_table, &file_to_filter, is_gzip, &output_file)?;
    }
    info!("Successfully filtered input file {:?}, written to {:?}, with bed {:?}", file_to_filter, output_file, bed_file);
    
    Ok(())
}