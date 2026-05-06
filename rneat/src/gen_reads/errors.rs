use std::{num::{ParseFloatError, ParseIntError}, str::ParseBoolError};

use common::{
    file_tools::{
        bam_writer::BamWriterError,
        bed_reader::BedReaderError,
        fasta_reader::FastaReaderError,
        fastq_tools::FastqToolsError,
        vcf_tools::VcfToolsError,
    },
    
    models::{
        fragment_length::FragmentModelError, 
        mutation_model::MutationModelError, 
        quality_scores::QualityModelError, 
        sequencing_error_model::SeqModelError,
        gc_bias_model::GcBiasModelError,
    }, 
    
    structs::{
        distributions::DistributionErrors, 
        fasta_map::FastaMapError,
        mutated_map::MutatedMapError,
    }
};
use common::rng::NeatRngError;
use thiserror::Error;
use serde_yml;

#[derive(Debug, Error)]
pub enum GenerateReadsError {
    #[error("Invalid CLI inputs: {0}")]
    CliError(String),
    #[error("NEAT generate-reads requires a reference to proceed.")]
    MissingReferenceError,
    #[error("Error processing Configuration for the run")]
    MainConfigurationError,
    #[error("Error while applying variants!")]
    ApplyVariantsError,
    #[error("Error reading config key: {0} accepts values of {1}")]
    ConfigReadError(String, String),
    #[error("Error generating configuration!")]
    ConfigError,
    #[error("Error generating fragments!")]
    GenerateFragmentsError,
    #[error("Error generating variants!")]
    GenerateVariantsError,
    #[error("Error while mutating fasta!")]
    MutateFastaError,
    #[error("Generate Reads Runner reported an error!")]
    RunnerError,
    #[error("File not found: {0}")]
    FileNotFound(String),
    #[error("IO Error called by NEAT: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Error generating fragment length model: {0}")]
    FragModelError(#[from] FragmentModelError),
    #[error("Error creating quality score model: {0}")]
    QualModelError(#[from] QualityModelError),
    #[error("Error creating mutation model: {0}")]
    MutModelError(#[from] MutationModelError),
    #[error("Error creating sequencing error model: {0}")]
    SeqModelError(#[from] SeqModelError),
    #[error("Error accessing FastaMap: {0}")]
    FastaMapError(#[from] FastaMapError),
    #[error("Error reading fasta file: {0}")]
    FastaReaderError(#[from] FastaReaderError),
    #[error("Error creating distributions: {0}")]
    GenReadsDistroError(#[from] DistributionErrors),
    #[error("Error sampling distro: {0}")]
    GenReadsRngError(#[from] NeatRngError),
    #[error("Error creating mutated map for block: {0}")]
    MutMapError(#[from] MutatedMapError),
    #[error("Error during configuration creation: {0}")]
    ConfigYamlError(#[from] serde_yml::Error),
    #[error("Configuration reader reported an error parsing an int: {0}")]
    ConfigParseIntError(#[from] ParseIntError),
    #[error("Configuration reader reported an error parsing a float: {0}")]
    ConfigParseFloatError(#[from] ParseFloatError),
    #[error("Configuration reader reported an error parsing a bool: {0}")]
    ConfigParseBoolError(#[from] ParseBoolError),
    #[error("Main threw an error: {0}")]
    MainError(#[from] Box<dyn std::error::Error>),
    #[error("FastqTools error: {0}")]
    FqToolsError(#[from] FastqToolsError),
    #[error("Bed reader error: {0}")]
    BedError(#[from] BedReaderError),
    #[error("Input VCF error: {0}")]
    InputVcfError(#[from] VcfToolsError),
    #[error("BAM writer error: {0}")]
    BamWriterError(#[from] BamWriterError),
    #[error("Sequence too short to process")]
    ShortSequence,
    #[error("GC Bias model threw an error: {0}")]
    BiasModelError(#[from] GcBiasModelError),
}