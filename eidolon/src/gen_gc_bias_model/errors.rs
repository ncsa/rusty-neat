use eidolon_core::file_tools::bam_reader::BamReaderError;
use eidolon_core::file_tools::bed_reader::BedReaderError;
use eidolon_core::file_tools::fasta_stream::FastaStreamError;
use eidolon_core::models::gc_bias_model::GcBiasModelError;
use thiserror::Error;

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
    #[error("Error reading FASTA: {0}")]
    FastaStreamError(#[from] FastaStreamError),
    #[error("Error reading BAM: {0}")]
    BamReaderError(#[from] BamReaderError),
}
