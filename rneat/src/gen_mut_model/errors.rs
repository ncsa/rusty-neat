use thiserror::Error;
use std::num::ParseIntError;
use common::{file_tools::{bed_reader::BedReaderError, fasta_reader::FastaReaderError, vcf_tools::VcfToolsError}, models::{mutation_model::MutationModelError, snp_trinuc_model::SnpTrinucError}, structs::fasta_map::FastaMapError};

#[derive(Error, Debug)]
pub enum GenMutationModelError {
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
    #[error("Error reading the vcf file: {0}")]
    VcfReaderError(#[from] VcfToolsError),
    #[error("Error generating trinucleotide counts: {0}")]
    TrinucCountError(String),
    #[error("Error reading data from FastaMap: {0}")]
    FastaMapError(#[from] FastaMapError),
    #[error("Error while reading in fasta file: {0}")]
    FastaReaderError(#[from] FastaReaderError),
    #[error("Variant location did not match reference: {0}")]
    BaseMismatch(String),
    #[error("Error building mutation model: {0}")]
    MutModelError(#[from] MutationModelError),
}