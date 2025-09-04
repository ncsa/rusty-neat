//! The mutation model
use std;
use bincode::de;
use serde::{Deserialize, Serialize};
use thiserror::Error;
use simple_rng::{NeatRng, NeatRngError};
use log::error;
use crate::models::lib::{model_reader, model_writer};
use crate::structs::transition_matrix::{TransitionMatrix, TransitionMatrixError};
use crate::structs::variants::{Variant, VariantError, VariantType};
use crate::structs::distributions::{DiscreteDistribution, DistributionErrors};
use crate::models::indel_model::{IndelModel, IndelModelError};
use crate::models::quality_scores::{QualityModelError, QualityScoreModel};
use crate::models::sequencing_error_model::{SeqModelError, SequencingErrorModel};
use crate::models::snp_trinuc_model::{SnpTrinucModel, SnpTrinucError};
use crate::structs::nucleotides::Nucleotide;

#[derive(Error, Debug)]
pub enum MutationModelError {
    #[error("Mutation model returned Error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("Error retrieving reference sequence block")]
    ReferenceRetrievalError,
    #[error("Mutation model IO Error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Error Generating mutations")]
    GenerateMutationError,
    #[error("Mutation model returned an error during variant creation: {0}")]
    VariantInitationError(#[from] VariantError),
    #[error("Mutation model returned an error during SNP initation: {0}")]
    SnpeGnerationError(#[from] SnpTrinucError),
    #[error("Mutation model returned error during transition matrix initation: {0}")]
    TransMatrixInitationError(#[from] TransitionMatrixError),
    #[error("Mutation model returned error during indel model initation: {0}")]
    IndelModelInitationError(#[from] IndelModelError),
    #[error("Mutation model returned error during quality model initiation: {0}")]
    QualityModelInitationError(#[from] QualityModelError),
    #[error("Mutation model returned error during sequence error model initation: {0}")]
    SeqErrModelError(#[from] SeqModelError),
    #[error("Mutation model returned an error with a distribution model: {0}")]
    DistributionError(#[from] DistributionErrors),
    #[error("Mutation model return and error: {0}")]
    ModelInitiationError(&'static str),
}

#[derive(Clone, Serialize, Deserialize)]
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

impl MutationModel {
    pub fn default() -> Result<Self, MutationModelError> {
        // Creating the default model based on the default for the original NEAT.
        let statistical_models = StatisticalModels::default()?;
        // Based on NA12878 data
        let mutation_rate = 0.0010987132390211135;
        // This was hard-coded in
        let homozygous_frequency = 0.01;
        // from NA12878 data
        let insertion_given_indel_rate = statistical_models.indel_model.insertion_probability;
        // from NA12878 data
        let snp_rate = 0.886404192662459;
        // Originally this was expressed as indel_fraction 0.05, So snp probs is 0.95
        // insertion frequency (relative to their being an indel variant) was 0.6.
        // This gives us, e.g.,
        //                 probability of insertion of 0.6 * 0.05 = 0.03 * 100 = 3
        //                 probability of deletion of  0.4 * 0.05 = 0.02 * 100 = 2
        //                 probability of SNP of       0.95              * 100 = 95
        let ins_rate = (1.0 - snp_rate) * (insertion_given_indel_rate);
        let del_rate = (1.0 - snp_rate) * (1.0 - insertion_given_indel_rate);
        assert!((snp_rate + ins_rate + del_rate - 1.0).abs() < 1e-10);
        let variant_weights = [
            snp_rate, // snp
            ins_rate, // insertion
            del_rate, //deletion
        ];

        // We'll use the inex to find the variant
        let variant_index = vec![0, 1, 2];
        let variant_dist = DiscreteDistribution::new(
            &Vec::from(variant_weights),
            &variant_index,
        )?;

        Ok(MutationModel {
            mutation_rate,
            homozygous_frequency,
            variant_dist,
            statistical_models,
        })
    }

    pub fn from_file(filename: &str) -> Result<Self, MutationModelError> {
        let data: MutationModel = model_reader(filename)?;
        Ok(data)
    }

    pub fn write_to_file(&self, filename: &str) -> Result<(), MutationModelError> {
        model_writer(self, filename)?;
        Ok(())
    }

    fn generate_genotype(&mut self, ploidy: usize, is_homozygous: bool) -> Result<Vec<usize>, MutationModelError> {
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
        // A pointer to the reference sequence to be mutated
        reference_sequence: &Vec<Nucleotide>,
        // This is the start position of the variant, the POS field of a VCF
        variant_location: usize,
        // The ploidy of this organism. We'll make determinations about heterzygous v homozygous based on this
        ploidy: usize,
        // The run's rng. Needed to determine type then generate the mutation
        rng: &mut NeatRng,
    ) -> Result<Variant, MutationModelError> {
        // Select a genotype for the variant
        let mut genotype= self.generate_genotype(
            ploidy, rng.gen_bool(self.homozygous_frequency)?
        )?;
        // Select a type of mutation.
        let index = self.variant_dist.sample(rng.random()?)?;
        let variant_type = VariantType::from(index);
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
                let alternate_base = self.statistical_models.snp_model
                    .generate_snp(rng.random()?, &trinuc_reference)?;

                let reference = vec![trinuc_reference[1]];
                let alternate = vec![alternate_base];

                (reference, alternate)
            },
            VariantType::Insertion => {
                let length = self.statistical_models.indel_model
                    .new_insert_length(rng.random()?)?;
                let insertion_vec = self.statistical_models.indel_model
                    .generate_random_insertion(length, rng)?;
                let reference = vec![
                    reference_sequence[variant_location].clone()
                ];
                let alternate = [reference.clone(), insertion_vec.clone()].concat();
                (reference, alternate)
            },
            VariantType::Deletion => {
                // todo Deletions are a bitch. I need to think about them.
                let length = self.statistical_models.indel_model.new_delete_length(rng.random()?)?;
                // +1 is so that we grab a base for the reference in the VCF. This is similar
                // to how we appended bases to the reference in the insertion model.
                let reference: Vec<Nucleotide> = reference_sequence
                    .get(variant_location..(variant_location + length as usize + 1))
                    .unwrap()
                    .iter()
                    .map(|item| *item)
                    .collect();
                let alternate: Vec<Nucleotide> = vec![reference[0].clone()];
                (reference, alternate)
            },
        };
        
        Ok(Variant::new(
            variant_type,
            variant_location,
            &reference,
            &alternate,
            &mut genotype,
        )?)
    }
}

#[derive(Clone, Serialize, Deserialize)]
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
    pub fn default() -> Result<Self, MutationModelError> {
        // use the default transition matrix, snp model and indel model
        let transition_matrix = TransitionMatrix::default()?;
        let snp_model = SnpTrinucModel::default_minimal()?;
        let indel_model = IndelModel::default()?;
        // todo update this line:
        let quality_score_model = QualityScoreModel::default()?;
        let sequencing_error_model = SequencingErrorModel::default()?;

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
    fn test_model_read_write() {
        let output_file: &'static str = "default_mutation_model.json.gz";
        let model: MutationModel = MutationModel::default().unwrap();
        assert_eq!(model.mutation_rate, 0.0010987132390211135);
        let result = model.write_to_file(output_file);
        assert_eq!(result.unwrap(), ());
        // fs::remove_file(output_file).unwrap();
    }
}
