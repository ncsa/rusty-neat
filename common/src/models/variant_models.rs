use crate::utils;
use crate::models;
use crate::variant_classes;

use utils::neat_rng::NeatRng;
use models::nucleotides::Nuc;
use models::transition_matrix::TransitionMatrix;
use variant_classes::variants::{Variant, VariantType};

// Statistical models are based on the models used in the original NEAT. NEAT made no attempt to
// distinguish between insertions or deletions, nor of different types of either, in terms of
// statistics. NOTE: any new statistical models would need to be created in the structs below and
// added to this list in order to be fully implemented.
pub struct StatisticalModels {
    // the model governing the single nucleotide polymorphisms for this run.
    pub(crate) snp_model: SnpModel,
    // The model governing indels for this run.
    pub(crate) indel_model: IndelModel,
    // The 4 x 4 matrix that shows the probability of one nucleotide transitioning to another.
    pub(crate) transition_matrix: TransitionMatrix,
}

impl StatisticalModels {
    pub fn new() -> Self {
        // use the default transition matrix, snp model and indel model
        let transition_matrix = TransitionMatrix::new();
        let snp_model = SnpModel::new();
        let indel_model = IndelModel::new();

        StatisticalModels {
            transition_matrix,
            snp_model,
            indel_model
        }
    }

    pub fn from(
        transition_matrix: TransitionMatrix,
        snp_model: SnpModel,
        indel_model: IndelModel
    ) -> Self {
        StatisticalModels {
            transition_matrix,
            snp_model,
            indel_model
        }
    }

    pub fn get_variant(
        &self,
        variant_type: VariantType,
        variant_location: usize,
        input_sequence: &Vec<Nuc>,
        rng: &mut NeatRng
    ) -> Vec<Nuc> {
        match variant_type {
            VariantType::SNP => { self.snp_model.generate_variant(
                input_sequence, variant_location, &self.transition_matrix, rng
            ) },
            VariantType::Indel => { self.indel_model.generate_variant(
                input_sequence, variant_location, &self.transition_matrix, rng
            ) },
        }
    }
}

// The following section are the models for each type of variant. In order to create the variant,
// we need to model its statistical property. NEAT included two types of variants: SNPs and Indels.
// Note that indels are actually two types, insertions and deletions, but they are usually classed
// together in the literature and are usually short. Insertions are often slips that cause sections
// to be duplicated, but NEAT made no attempt to distinguish between the types of insertions or
// deletions, since they are treated similarly by variant calling software.

// All new variants must use the Mutate implementation to work.
pub trait Mutate {
    fn generate_variant(
        &self,
        chromosome: &str,
        input_sequence: &Vec<Nuc>,
        variant_location: usize,
        transition_matrix: &TransitionMatrix,
        rng: &mut NeatRng
    ) -> Variant;

}

#[derive(Debug)]
struct SnpModel {
    // These are the 16 possible patterns for trinucleotides, each representing 4 trinucleotides for
    // a total of 64 trinucleotide combinations. The usize is the weight of that particular group
    // and the transition matrix is the chance of mutating the middle base from A, C, T, or G to a
    // different base (4x4 matrix with 0s on the diagonal).
    a_a: (usize, TransitionMatrix),
    a_c: (usize, TransitionMatrix),
    a_g: (usize, TransitionMatrix),
    a_t: (usize, TransitionMatrix),
    c_a: (usize, TransitionMatrix),
    c_c: (usize, TransitionMatrix),
    c_g: (usize, TransitionMatrix),
    c_t: (usize, TransitionMatrix),
    g_a: (usize, TransitionMatrix),
    g_c: (usize, TransitionMatrix),
    g_g: (usize, TransitionMatrix),
    g_t: (usize, TransitionMatrix),
    t_a: (usize, TransitionMatrix),
    t_c: (usize, TransitionMatrix),
    t_g: (usize, TransitionMatrix),
    t_t: (usize, TransitionMatrix),
}

impl SnpModel {
    pub fn new() -> Self {
        // Creating the default trinuc bias model for snps. In this model, all trinucleotides
        // mutate with equal probability and middle base mutates with the same probability no matter
        // the context (the default transition matrix).
        SnpModel {
            a_a: (1, TransitionMatrix::new()),
            a_c: (1, TransitionMatrix::new()),
            a_g: (1, TransitionMatrix::new()),
            a_t: (1, TransitionMatrix::new()),
            c_a: (1, TransitionMatrix::new()),
            c_c: (1, TransitionMatrix::new()),
            c_g: (1, TransitionMatrix::new()),
            c_t: (1, TransitionMatrix::new()),
            g_a: (1, TransitionMatrix::new()),
            g_c: (1, TransitionMatrix::new()),
            g_g: (1, TransitionMatrix::new()),
            g_t: (1, TransitionMatrix::new()),
            t_a: (1, TransitionMatrix::new()),
            t_c: (1, TransitionMatrix::new()),
            t_g: (1, TransitionMatrix::new()),
            t_t: (1, TransitionMatrix::new()),
        }
    }
    fn generate_bias(&self, input_sequence: &Vec<Nuc>) -> Vec<usize> {
        todo!()
        // We need some way to use this model to bias positions of SNPs, but it's not clear yet how.
    }
}

impl Mutate for SnpModel {
    fn generate_variant(
        &self,
        chromosome: &str,
        input_sequence: &Vec<Nuc>,
        variant_location: usize,
        transition_matrix: &TransitionMatrix,
        rng: &mut NeatRng
    ) -> Variant {
        todo!()
        // The sequence for this model must contain the nucelotide before and after for context.

    }

}

#[derive(Debug)]
struct IndelModel {
    // Based what was in the original NEAT
    insertion_probability: f64,
    ins_lengths: Vec<usize>,
    ins_weights: Vec<usize>,
    del_lengths: Vec<usize>,
    del_weights: Vec<usize>,
}

impl Mutate for IndelModel {
    fn generate_variant(
        &self,
        chromosome: &str,
        input_sequence: &Vec<Nuc>,
        variant_location: usize,
        transition_matrix: &TransitionMatrix,
        rng: &mut NeatRng
    ) -> Variant {
        todo!()
    }
}

impl IndelModel {
    pub fn new() -> Self {
        // Default Indel model from the original NEAT
        let insertion_probability = 0.6;
        let ins_lengths: Vec<usize> = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let ins_weights: Vec<usize> = vec![10, 10, 20, 5, 5, 5, 5, 5, 5, 5];
        let del_lengths: Vec<usize> = vec![1, 2, 3, 4, 5];
        let del_weights: Vec<usize> = vec![3, 2, 2, 2, 1];

        IndelModel {
            insertion_probability,
            ins_lengths,
            ins_weights,
            del_lengths,
            del_weights,
        }
    }
}