use thiserror::Error;
use common::models::{
    quality_scores::QualityModelError,
    sequencing_error_model::SeqModelError,
};

#[derive(Error, Debug)]
pub enum GenSeqErrorModelError {
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Configuration error: {0}")]
    ConfigurationError(String),
    #[error("CLI error: {0}")]
    CliError(String),
    #[error("Malformed FASTQ: {0}")]
    MalformedFastq(String),
    #[error("Quality score model error: {0}")]
    QualModelError(#[from] QualityModelError),
    #[error("Sequencing error model error: {0}")]
    SeqModelError(#[from] SeqModelError),
}
