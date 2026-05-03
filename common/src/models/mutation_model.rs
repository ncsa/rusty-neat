//! The mutation model
use std::{self, collections::HashMap};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use flate2::read::GzDecoder;
use thiserror::Error;
use simple_rng::{NeatRng, NeatRngError};
use log::error;

use crate::{models::snp_trinuc_model::TrinucFrame, structs::{
    distributions::{
        DiscreteDistribution, 
        DistributionErrors
    }, 
    nucleotides::{
        allowed_vec, 
        Nucleotide
    }, 
    transition_matrix::{
        TransitionMatrix, 
        TransitionMatrixError
    }, 
    variants::{
        Variant, 
        VariantError, 
        VariantType
    }
}};
use crate::models::{
    indel_model::{IndelModel, IndelModelError},
    quality_scores::{QualityModelError},
    sequencing_error_model::{SeqModelError},
    snp_trinuc_model::{SnpTrinucModel, SnpTrinucError},
    lib::{model_reader, model_writer},
};

#[derive(Error, Debug)]
pub enum MutationModelError {
    #[error("Error in inputs to MutationModel")]
    InputError,
    #[error("Mutation model returned RNG error: {0}")]
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
    ModelInitiationError(String),
    #[error("Error reading Mutation Model from file {0}")]
    SerdeError(#[from] serde_json::Error),
}

fn allowed_variant_types() -> Vec<VariantType> { 
    vec![
        VariantType::SNP,
        VariantType::Insertion,
        VariantType::Deletion,
    ]
}

#[derive(Clone, Debug, Serialize, Deserialize)]
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
    pub variant_dist: DiscreteDistribution<VariantType>,
    // These hold all the statistical model data we need to apply the mutations with this model
    statistical_models: StatisticalModels,
    // Store a reference to the NeatRng for this run
}

static DATA_FILE: &'static [u8] = include_bytes!("model_data/default_mutation_model.json.gz");

impl MutationModel {
    pub fn from_raw_data(
        average_mutation_rate: f64,
        homozygous_frequency: f64,
        variant_probs: Vec<f64>,
        snp_transition_frequency: HashMap<(Nucleotide, Nucleotide), f64>,
        trinuc_frequency: HashMap<TrinucFrame, f64>,
        trinuc_transition_frequency: HashMap<(TrinucFrame, TrinucFrame), f64>,
        ins_lengths: Vec<usize>,
        ins_weights: Vec<f64>,
        del_lengths: Vec<usize>,
        del_weights: Vec<f64>,
    ) -> Result<Self, MutationModelError> {
        // Inputs:
        // average_mutation_rate: Average rate by total reference (or bed track) length
        //    of a base mutating
        // homozygous_frequency: Average rate that, if a mutation occurs, it is homozygous
        // variant_probs: Relative probabilities of each type of allowed variant (SNP, INS, DEL)
        // snp_transition_frequency: Probability that one base mutates to another, 
        //    independent of context
        // trinuc_frequency: The relative frequency of each trinucleotide in the dataset
        //    mutating to anything else (to help detect if some trinucs are more volatile 
        //    than others)
        // trinuc_transition_frequency: The probability of trinuc A mutating into trinuc B,
        //    as seen in data. Some may not have been observed and will be assigned a very
        //    low probability.
        //
        // custom tranisition matrix for snps
        let mut temp_trans_matrix: [[f64; 4]; 4] = [
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ];
        for (key, value) in snp_transition_frequency {
            temp_trans_matrix[key.0 as usize][key.1 as usize] = value;
        }
        let transition_matrix = TransitionMatrix::from(
            temp_trans_matrix[Nucleotide::A as usize],
            temp_trans_matrix[Nucleotide::C as usize],
            temp_trans_matrix[Nucleotide::G as usize],
            temp_trans_matrix[Nucleotide::T as usize],
        )?;
        // build transition matrices from data for snps and trinucs
        let snp_trinuc_model = SnpTrinucModel::from_raw_data(
            trinuc_frequency,
            trinuc_transition_frequency,
        )?;
        // Probability, given that it is an indel, that it is an insertion.
        // Fall back to the default indel model when no indel data was observed.
        let indel_denom = variant_probs[1] + variant_probs[2];
        let indel_model = if indel_denom == 0.0 || ins_lengths.is_empty() || del_lengths.is_empty() {
            IndelModel::default()?
        } else {
            let insertion_probability = variant_probs[1] / indel_denom;
            IndelModel::from_raw_data(
                insertion_probability,
                ins_lengths,
                ins_weights,
                del_lengths,
                del_weights,
            )?
        };
        let statistical_models = StatisticalModels {
            transition_matrix,
            indel_model,
            snp_trinuc_model,
        };
        // SNP, Insertion, Deletion, respectively.
        if variant_probs.len() != 3 {
            error!("Input to mutation model from raw incorrect. Variant probs should have a length of 3: SNP, INS, DEL");
            return Err(MutationModelError::InputError)
        }
        let variant_distribution = DiscreteDistribution::new(
            &variant_probs,
            &(allowed_variant_types()),
        )?;
        Ok(MutationModel {
            mutation_rate: average_mutation_rate,
            variant_dist: variant_distribution,
            homozygous_frequency,
            statistical_models,
        })
    }

    pub fn default() -> Result<Self, MutationModelError> {
        // Creating the default model based on the default for the original NEAT.
        let reader = GzDecoder::new(DATA_FILE);
        let data: MutationModel = serde_json::from_reader(reader)
            .map_err(MutationModelError::SerdeError)?;
        Ok(data)
    }

    pub fn from_file(filename: &PathBuf) -> Result<Self, MutationModelError> {
        let data: MutationModel = model_reader(filename)?;
        Ok(data)
    }

    pub fn write_to_file(&self, filename: &PathBuf) -> Result<(), MutationModelError> {
        model_writer(self, filename)?;
        Ok(())
    }

    fn generate_genotype(&self, ploidy: usize, is_homozygous: bool) -> Result<Vec<usize>, MutationModelError> {
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
        &self,
        // A pointer to the reference sequence to be mutated
        reference_sequence: &Vec<Nucleotide>,
        // This is the start position of the variant, the POS field of a VCF
        variant_location: usize,
        // The ploidy of this organism. We'll make determinations about heterzygous
        // v homozygous based on this
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
        let mut variant_type = VariantType::from(index);
        let ref_base = reference_sequence[variant_location];
        // For insertions and SNPs, the reference will be the ref_base
        let mut reference: Vec<Nucleotide> = vec![ref_base];
        // For deletions, the alternato will be the ref_base
        let mut alternate: Vec<Nucleotide> = vec![ref_base];
        match variant_type {
            VariantType::SNP => {
                if (variant_location < 1) || ((variant_location + 1) >= reference_sequence.len()) {
                    // We don't want to accidentally go out of bounds, so we'll 
                    // just pick a random base.
                    alternate = pick_random_snp(ref_base, rng)?;
                } else {
                    // In this case we use the trinculeotide model
                    let trinuc_reference = [
                        check_base(reference_sequence[variant_location-1]),
                        check_base(reference_sequence[variant_location]),
                        check_base(reference_sequence[variant_location+1]),
                    ];
                    let alternate_base = self.statistical_models.snp_trinuc_model
                        .generate_snp(rng.random()?, &trinuc_reference)?;

                    alternate = vec![alternate_base];
                }
            },
            VariantType::Insertion => {
                let length = self.statistical_models.indel_model
                    .new_insert_length(rng.random()?)?;
                let insertion_vec = self.statistical_models.indel_model
                    .generate_random_insertion(length, rng)?;
                alternate = [reference.clone(), insertion_vec.clone()].concat();
            },
            VariantType::Deletion => {
                // todo Deletions are a bitch. I need to think about them.
                let length = self.statistical_models.indel_model.new_delete_length(rng.random()?)?;
                // +1 is so that we grab a base for the reference in the VCF. This is similar
                // to how we appended bases to the reference in the insertion model.
                if (variant_location + length as usize + 1) > reference_sequence.len() {
                    // Too close to the end, so let's skip this
                    let ref_base = reference_sequence[variant_location];
                    alternate = pick_random_snp(ref_base, rng)?;
                    variant_type = VariantType::SNP;
                } else {
                    reference = reference_sequence
                        .get(variant_location..(variant_location + length as usize + 1))
                        .unwrap()
                        .to_vec();
                }
            },
            _ => {
                panic!("Unsupported option for generating mutation")
            }
        };
        
        if reference.is_empty() || alternate.is_empty() || reference == alternate {
            error!("Malformed mutation: {:?} {:?} {:?}", variant_type, reference, alternate);
            return Err(MutationModelError::GenerateMutationError)
        }
        Ok(Variant::new(
            variant_type,
            variant_location,
            &reference,
            &alternate,
            &mut genotype,
        )?)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct StatisticalModels {
    // This struct links the Mutation model to the other statistical models for modeling different
    // variant types. If new variant types are added, then we will need to expand this struct to
    // include them.
    transition_matrix: TransitionMatrix,
    indel_model: IndelModel,
    snp_trinuc_model: SnpTrinucModel,
}

fn pick_random_snp(ref_base: Nucleotide, rng: &mut NeatRng)-> Result<Vec<Nucleotide>, MutationModelError> {
    // Grab allowed nucs
    let mut allowed: Vec<Nucleotide> = allowed_vec();
    // Throw out the reference base
    allowed.retain(|&x| x != ref_base);
    // Pick one of the remaining
    let alt = rng.choose(&allowed)?;
    // sanity check
    assert_ne!(alt, ref_base);
    Ok(vec![alt])
}

fn check_base(nuc: Nucleotide) -> Nucleotide {
    // For now we'll replace occasional N's with A's.
    if nuc.is_masked() {
        return nuc.get_unmasked_base()
    }
    match nuc {
        Nucleotide::N => Nucleotide::A,
        _ => nuc,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use crate::models::snp_trinuc_model::TrinucFrame;

    #[test]
    fn test_model_read_write() {
        let output_file: PathBuf = PathBuf::from("test.json.gz");
        let model: MutationModel = MutationModel::default().unwrap();
        assert_eq!(model.mutation_rate, 0.0010987132390211135);
        let result = model.write_to_file(&output_file);
        assert_eq!(result.unwrap(), ());
        fs::remove_file(output_file).unwrap();
    }

    #[test]
    fn test_from_raw_data_no_indels() {
        // Exercises the indel_denom == 0 || empty vecs fallback to IndelModel::default()
        let snp_trans: HashMap<(Nucleotide, Nucleotide), f64> = HashMap::from([
            ((Nucleotide::A, Nucleotide::C), 0.25),
            ((Nucleotide::A, Nucleotide::G), 0.5),
            ((Nucleotide::A, Nucleotide::T), 0.25),
        ]);
        let trinuc_freq: HashMap<TrinucFrame, f64> = HashMap::new();
        let trinuc_trans: HashMap<(TrinucFrame, TrinucFrame), f64> = HashMap::new();
        let result = MutationModel::from_raw_data(
            0.001,
            0.5,
            vec![1.0, 0.0, 0.0], // SNP only, no indels
            snp_trans,
            trinuc_freq,
            trinuc_trans,
            vec![], // no insertion data
            vec![],
            vec![], // no deletion data
            vec![],
        );
        assert!(result.is_ok(), "Expected Ok with no indels, got {:?}", result);
        let model = result.unwrap();
        assert_eq!(model.mutation_rate, 0.001);
    }
}
