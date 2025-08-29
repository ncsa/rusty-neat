//! The mutation model
use std;
use thiserror::Error;
use simple_rng::{NeatRng, NeatRngError};
use log::error;
use crate::structs::transition_matrix::{TransitionMatrix, TransitionMatrixError};
use crate::structs::variants::{Variant, VariantError, VariantType};
use crate::structs::distributions::DiscreteDistribution;
use crate::models::indel_model::{IndelModel, IndelModelError};
use crate::models::quality_scores::{QualityModelError, QualityScoreModel};
use crate::models::sequencing_error_model::{SeqModelError, SequencingErrorModel};
use crate::models::snp_trinuc_model::{SnpTrinucModel, SnpTrinucError};


#[derive(Error, Debug)]
pub enum MutationModelError {
    #[error("Mutation model returned Error: {0}")]
    RngError(NeatRngError),
    #[error("Error retrieving reference sequence block")]
    ReferenceRetrievalError,
    #[error("Mutation model IO Error: {0}")]
    IoError(std::io::Error),
    #[error("Error Generating mutations")]
    GenerateMutationError,
    #[error("Mutation model returned an error during variant creation: {0}")]
    VariantInitationError(VariantError),
    #[error("Mutation model returned an error during SNP initation: {0}")]
    SnpeGnerationError(SnpTrinucError),
    #[error("Mutation model returned error during transition matrix initation: {0}")]
    TransMatrixInitationError(TransitionMatrixError),
    #[error("Mutation model returned error during indel model initation: {0}")]
    IndelModelInitationError(IndelModelError),
    #[error("Mutation model returned error during quality model initiation: {0}")]
    QualityModelInitationError(QualityModelError),
    #[error("Mutation model returned error during sequence error model initation: {0}")]
    SeqErrModelError(SeqModelError)
}

impl From<QualityModelError> for MutationModelError {
    fn from(error: QualityModelError) -> Self {
        MutationModelError::QualityModelInitationError(error)
    }
}

impl From<IndelModelError> for MutationModelError {
    fn from(error: IndelModelError) -> Self {
        MutationModelError::IndelModelInitationError(error)
    }
}

impl From<TransitionMatrixError> for MutationModelError {
    fn from(error: TransitionMatrixError) -> Self {
        MutationModelError::TransMatrixInitationError(error)
    }
}

impl From<SeqModelError> for MutationModelError {
    fn from(error: SeqModelError) -> Self {
        MutationModelError::SeqErrModelError(error)
    }
}

impl From<SnpTrinucError> for MutationModelError {
    fn from(error: SnpTrinucError) -> Self {
        MutationModelError::SnpeGnerationError(error)
    }
}

impl From<VariantError> for MutationModelError {
    fn from(error: VariantError) -> Self {
        MutationModelError::VariantInitationError(error)
    }
}

impl From<NeatRngError> for MutationModelError {
    fn from(error: NeatRngError) -> Self {
        MutationModelError::RngError(error)
    }
}

impl From<std::io::Error> for MutationModelError {
    fn from(error: std::io::Error) -> Self {
        MutationModelError::IoError(error)
    }
}

#[derive(Clone)]
pub struct MutationModel {
    // This is the model for mutations, the same construct used by the python version, basically.
    // The idea is to store all the data we need to create mutations in one construct.
    //
    // mutation rate is the average rate of mutation for the dataset. It is used to calculate how
    // many variants to add to the dataset.
    pub mutation_rate: f64,
    // Homozygous frequency is the fraction of mutations that are homozygous, meaning the alternate
    // was inherited from both parents (in the case of humans). The definition of "homozygous"
    // is ambiguous in polyploid organisms. We'll take it to mean "on all ploids"
    pub homozygous_frequency: f64,
    // If a variant occurs, this is how frequently each type occurs. There are only 3 types at the
    // moment, but this will get expanded out in time. This data item will need to be expanded to
    // account for it. The canonical order is:
    //    1. SNP
    //    2. Insertion
    //    3. Deletion
    pub variant_dist: DiscreteDistribution,
    // These hold all the statistical model data we need to apply the mutations with this model
    statistical_models: StatisticalModels,
    // Store a reference to the NeatRng for this run
}

const VARIANT_TYPES: [VariantType; 3] = [
    VariantType::SNP,       // 1
    VariantType::Insertion, // 2
    VariantType::Deletion,  // 3
];

impl MutationModel {
    pub fn default() -> Result<Self, MutationModelError> {
        // Creating the default model based on the default for the original NEAT.
        let mutation_rate = 0.001;
        let homozygous_frequency = 0.01;
        // Originally this was expressed as indel_fraction 0.05. We're using weights for
        // readability. This gets converted to a cumulative density function by the
        // DiscreteDistribution struct. This would have been 0.95 and 0.05, combine with an
        // insertion frequency (relative to their being an indel variant) was 0.6.
        // This gives us a probability of insertion of 0.6 * 0.05 = 0.03 * 100 = 3
        //                 probability of deletion of  0.4 * 0.05 = 0.02 * 100 = 2
        //                 probability of SNP of       0.95              * 100 = 95
        let variant_weights = [
            0.95, // snp
            0.03, // insertion
            0.02, //deletion
        ];
        let variant_dist = DiscreteDistribution::new(&Vec::from(variant_weights))?;
        let statistical_models = StatisticalModels::new()?;

        Ok(MutationModel {
            mutation_rate,
            homozygous_frequency,
            variant_dist,
            statistical_models,
        })
    }

    pub fn new(
        mutation_rate: f64, 
        homozygous_frequency: f64, 
        variant_weights: Vec<f64>,
    ) -> Result<Self, MutationModelError> {
        let variant_dist = DiscreteDistribution::new(&Vec::from(variant_weights))?;
        let statistical_models = StatisticalModels::new()?;

        Ok(MutationModel {
            mutation_rate,
            homozygous_frequency,
            variant_dist,
            statistical_models,
        })
    }

    pub fn from_file(filename: &str, rng: &mut NeatRng) -> Result<Self, MutationModelError> {
        todo!()
        // Might want to store this data in a json or something for later use, 
        // probably using serde_json.
    }

    pub fn write_file(&self, filename: &str) -> Result<(), MutationModelError> {
        todo!()
    }

    fn generate_genotype(&mut self, ploidy: usize, is_homozygous: bool) -> Result<Vec<u8>, MutationModelError> {
        // "Homozygous" is ambiguous for polyploid organisms, so we'll just take "heterozygous" to
        // mean roughly half the reads will have the variant, to keep it simple
        // The is_homozygous flag is expected to be randomly determined in practice
        if is_homozygous {
            Ok(vec![1; ploidy])
        } else {
            let num_ploids = ploidy/2;
            Ok([vec![1; num_ploids], vec![0; ploidy-num_ploids]].concat())
        }
    }

    pub fn generate_mutation(
        &mut self,
        reference_sequence: &Vec<u8>,
        variant_location: usize,
        ploidy: usize,
        rng: &mut NeatRng,
    ) -> Result<Variant, MutationModelError> {
        // Select a genotype for the variant
        let genotype= self.generate_genotype(ploidy, rng.gen_bool(self.homozygous_frequency)?)?;
        // Select a type of mutation.
        let index = self.variant_dist.sample(rng.random()?)?;
        let variant_type = if index > VARIANT_TYPES.len() {
            error!("Weird result from sampling variant type: {}", index);
            return Err(MutationModelError::GenerateMutationError)
        } else {
            VARIANT_TYPES[index]
        };
        // todo figure out which block to mutate
        //   figure out the mutation to add
        //   Add a function in FastaMap to add a block to a vector
        let (reference, alternate) = match variant_type {
            VariantType::SNP => {
                let trinuc_reference = [
                    reference_sequence[variant_location-1],
                    reference_sequence[variant_location],
                    reference_sequence[variant_location+1],
                ];
                let alternate_base = self.statistical_models.snp_model.generate_snp(rng.random()?, &trinuc_reference)?;

                let reference = vec![trinuc_reference[1]];
                let alternate = vec![alternate_base];

                (reference, alternate)
            },
            VariantType::Insertion => {
                let length = self.statistical_models.indel_model.new_insert_length()?;
                let insertion_vec = generate_random_insertion(length)?;
                let reference = vec![
                    reference_sequence.get(variant_location)?.clone()
                ];
                let alternate = [reference.clone(), insertion_vec.clone()].concat();
                (reference, alternate)
            },
            VariantType::Deletion => {
                // todo Deletions are a bitch. I need to think about them.
                let length = self.statistical_models.indel_model.new_delete_length()?;
                // +1 is so that we grab a base for the reference in the VCF. This is similar
                // to how we appended bases to the reference in the insertion model.
                let reference: Vec<u8> = reference_sequence
                    .get(variant_location..(variant_location + length as usize + 1))
                    .unwrap()
                    .iter()
                    .map(|item| *item)
                    .collect();
                let alternate: Vec<u8> = vec![reference[0].clone()];
                (reference, alternate)
            },
        };
        
        Ok(Variant::new(
            variant_type,
            variant_location,
            &reference,
            &alternate,
            &genotype,
        )?)
    }
}

#[derive(Clone)]
struct StatisticalModels {
    // This struct links the Mutation model to the other statistical models for modeling different
    // variant types. If new variant types are added, then we will need to expand this struct to
    // include them.
    transition_matrix: TransitionMatrix,
    quality_score_model: QualityScoreModel,
    sequencing_error_model: SequencingErrorModel,
    indel_model: IndelModel,
    snp_model: SnpTrinucModel,
}

impl StatisticalModels {
    pub fn new() -> Result<Self, MutationModelError> {
        // use the default transition matrix, snp model and indel model
        let transition_matrix = TransitionMatrix::default()?;
        let snp_model = SnpTrinucModel::default(self.rng)?;
        let indel_model = IndelModel::default(self.rng)?;
        let quality_score_model = QualityScoreModel::default(self.rng)?;
        let sequencing_error_model = SequencingErrorModel::default(self.rng)?;

        Ok(StatisticalModels {
            transition_matrix,
            quality_score_model,
            snp_model,
            indel_model,
            sequencing_error_model,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_transition_matrix_from_weights() {
        let a_weights = vec![0.0, 20.0, 1.0, 20.0];
        let c_weights = vec![20.0, 0.0, 1.0, 1.0];
        let g_weights = vec![1.0, 1.0, 0.0, 20.0];
        let t_weights = vec![20.0, 1.0, 20.0, 0.0];
        TransitionMatrix::from(a_weights, c_weights, g_weights, t_weights);
        // It actually mutates the base
        // assert_ne!(test_model.generate_mutation
        // (
        //     "chr1",
        //     VariantType::SNP,
        //     0,
        //     &vec![A],
        //     &mut rng
        // ), vec![A]);
        // assert_ne!(test_model.generate_mutation
        // (
        //     "chr1",
        //     VariantType::SNP,
        //     0,
        //     &vec![C],
        //     &mut rng
        // ), vec![C]);
        // assert_ne!(test_model.generate_mutation
        // (
        //     "chr1",
        //     VariantType::SNP,
        //     0,
        //     &vec![G],
        //     &mut rng
        // ), vec![G]);
        // assert_ne!(test_model.generate_mutation
        // (
        //     "chr1",
        //     VariantType::SNP,
        //     0,
        //     &vec![T],
        //     &mut rng
        // ), vec![T]);
        // // It gives back N when you give it N
        // assert_eq!(test_model.mutate(
        //     VariantType::SNP, 0, &vec![N], &mut rng), vec![N]
        // );
        // // test alternate mutation methods
        // let mutation = test_model.generate_mutation
        // (
        //     "chr1",
        //     VariantType::SNP,
        //     1,
        //     &vec![A, C, G],
        //     &mut rng
        // );
        // assert_eq!(mutation[0], A);
        // assert_eq!(mutation[2], G);
        // assert_ne!(mutation[3], C);
        // let mutation = test_model.generate_mutation
        // (
        //     "chr1",
        //     VariantType::Indel,
        //     0,
        //     &vec![A, C, C, G, T, T, A, C, G],
        //     &mut rng
        // );
        // // todo something with this mutation variable
    }
}
