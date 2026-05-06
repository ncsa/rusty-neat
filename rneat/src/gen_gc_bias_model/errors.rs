use thiserror::Error;
use common::file_tools::bed_reader::BedReaderError;

#[derive(Error, Debug)]
pub enum GenGcBiasModelError {
    #[error("Error with GC Bias configuration: {0}")]
    ConfigError(String),
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Serde error: {0}")]
    SerdeError(#[from] serde_yml::Error),
    #[error("Error reading bed: {0}")]
    BedReadError(#[from] BedReaderError)
}