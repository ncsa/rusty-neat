use common::{models::{
    fragment_length::FragmentModelError, mutated_map::MutatedMapError, mutation_model::MutationModelError, quality_scores::QualityModelError, sequencing_error_model::SeqModelError
}, structs::{distributions::DistributionErrors, fasta_map::FastaMapError}};
use simple_rng::NeatRngError;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum GenerateReadsErrors {
    #[error("Error while applying variants!")]
    ApplyVariantsError,
    #[error("Error generating configuration!")]
    ConfigError,
    #[error("Error generating reads!")]
    GenerateReadsError,
    #[error("Error generating variants!")]
    GenerateVariantsError,
    #[error("Error while mutating fasta!")]
    MutateFastaError,
    #[error("Generate Reads Runner reported an error!")]
    RunnerError,
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
    #[error("Error creating distributions: {0}")]
    GenReadsDistroError(#[from] DistributionErrors),
    #[error("Error sampling distro: {0}")]
    GenReadsRngError(#[from] NeatRngError),
    #[error("Error creating mutated map for block: {0}")]
    MutMapError(#[from] MutatedMapError)
}