use thiserror::Error;

use crate::gen_frag_length_model::errors::GenFragLengthModelError;
use crate::gen_gc_bias_model::errors::GenGcBiasModelError;
use eidolon_core::file_tools::bam_reader::BamReaderError;

#[derive(Error, Debug)]
pub enum GenBamModelsError {
    #[error("Error with gen-bam-models configuration: {0}")]
    ConfigError(String),
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Serde error: {0}")]
    SerdeError(#[from] serde_yml::Error),
    #[error("Error reading BAM: {0}")]
    BamReaderError(#[from] BamReaderError),
    #[error("Fragment length model error: {0}")]
    FragLengthError(#[from] GenFragLengthModelError),
    #[error("GC bias model error: {0}")]
    GcBiasError(#[from] GenGcBiasModelError),
}
