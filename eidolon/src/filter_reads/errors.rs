use eidolon_core::file_tools::bed_reader::BedReaderError;
use std::num::ParseIntError;
use thiserror::Error;

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
