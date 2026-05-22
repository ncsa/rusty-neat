use thiserror::Error;
use common::file_tools::bam_reader::BamReaderError;
use common::file_tools::bed_reader::BedReaderError;
use common::file_tools::fasta_stream::FastaStreamError;
use common::models::gc_bias_model::GcBiasModelError;

#[derive(Error, Debug)]
pub enum GenGcBiasModelError {
    #[error("Error with GC Bias configuration: {0}")]
    ConfigError(String),
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Serde error: {0}")]
    SerdeError(#[from] serde_yml::Error),
    #[error("Error reading bed: {0}")]
    BedReadError(#[from] BedReaderError),
    #[error("Error building GC bias model: {0}")]
    GcBiasModelError(#[from] GcBiasModelError),
    #[error("Error parsing coverage file: {0}")]
    CoverageParseError(String),
    #[error("Error reading FASTA: {0}")]
    FastaStreamError(#[from] FastaStreamError),
    #[error("Error reading BAM: {0}")]
    BamReaderError(#[from] BamReaderError),
}