use common::file_tools::{bed_reader::BedReaderError, vcf_tools::VcfToolsError};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum CompareVcfsError {
    #[error("I/O error during compare-vcfs: {0}")]
    IoError(#[from] std::io::Error),
    #[error("YAML parse error: {0}")]
    YamlError(#[from] serde_yml::Error),
    #[error("JSON serialization error: {0}")]
    JsonError(#[from] serde_json::Error),
    #[error("Invalid Configuration inputs: {0}")]
    ConfigurationError(String),
    #[error("Cannot overwrite existing file: {0}")]
    OverwriteFileError(String),
    #[error("Error reading the VCF file: {0}")]
    VcfReaderError(#[from] VcfToolsError),
    #[error("The bed reader returned an error: {0}")]
    BedReaderErr(#[from] BedReaderError),
}
