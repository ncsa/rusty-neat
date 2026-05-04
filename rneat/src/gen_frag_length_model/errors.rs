use thiserror::Error;
use common::{
    file_tools::bam_reader::BamReaderError,
    models::fragment_length::FragmentModelError,
};

#[derive(Error, Debug)]
pub enum GenFragLengthModelError {
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Configuration error: {0}")]
    ConfigurationError(String),
    #[error("BAM/SAM reader error: {0}")]
    BamReaderError(#[from] BamReaderError),
    #[error("Fragment model error: {0}")]
    FragModelError(#[from] FragmentModelError),
    #[error("No valid fragment lengths found in input file; check that the file contains paired reads")]
    EmptyData,
    #[error("No data passed the min_reads filter; try lowering min_reads or set it to 0 to disable")]
    FilteredToEmpty,
}