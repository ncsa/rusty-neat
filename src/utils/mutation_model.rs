use std::collections::HashMap;
use log::debug;
use rand::distributions::{Distribution, WeightedIndex};
use utils::neat_rng::NeatRng;
use utils::nucleotides::Nuc;
use utils::variants::VariantType;

pub struct MutationModel {
    // This is the model for mutations, the same construct used by the python version, basically.
    //
    // mutation rate is the average rate of mutation for the dataset. It is used to calculate how
    // many variants to add to the dataset.
    mutation_rate: f64,
    // Homozygous frequency is the fraction of mutations that are homozygous, meaning the alternate
    // was inherited from both parents (in the case of humans). The definition of "homozygous"
    // is ambiguous in polyploid organisms. We'll take it to mean "on all ploids"
    homozygous_frequency: f64,
    // If a variant occurs, this is the probability that it will be a SNP, the most common type
    // of variant. There are only 2 types at the moment, but this will get expanded out in time.
    variant_probs: HashMap<VariantType, usize>,
    // These hold all the statistical model data we need to apply the mutations with this model
    statistical_models: StatisticalModels,
}
impl MutationModel {
    pub fn new() -> Self {
        // Creating the default model based on the default for the original NEAT.
        let mutation_rate = 0.001;
        let homozygous_frequency = 0.01;
        // Originally this was expressed as indel_fraction 0.05. To make it easier to sample, we
        // will use weights. This would have been 0.95 and 0.05, so we divided both by 0.05 to get
        // 19 and 1 as our weights.
        let variant_probs = HashMap::from([
            (VariantType::SNP, 19),
            (VariantType::Indel, 1),
        ]);

        let statistical_models = StatisticalModels::new();

        MutationModel {
            mutation_rate,
            homozygous_frequency,
            variant_probs,
            statistical_models,
        }
    }
    // We may need several ways to create the mutation model. Here is one for testing.
    pub fn from_transition_matrix(transition_matrix: TransitionMatrix) -> Self {
        let mut model = MutationModel::new();
        model.statistical_models.transition_matrix = transition_matrix;
        model
    }

    pub fn generate_mutation(
        &self, variant_type: VariantType, input_sequence: &Vec<Nuc>, rng: &mut NeatRng
    ) -> Vec<Nuc> {
        self.statistical_models.get_variant(
            variant_type: VariantType, input_sequence: &Vec<Nuc>, rng: &mut NeatRng
        )
    }
}

// Statistical models are based on the models used in the original NEAT. NEAT made no attempt to
// distinguish between insertions or deletions, nor of different types of either, in terms of
// statistics. NOTE: any new statistical models would need to be created in the structs below and
// added to this list in order to be fully implemented.
pub struct StatisticalModels {
    // The 4 x 4 matrix that shows the probability of one nucleotide transitioning to another.
    transition_matrix: TransitionMatrix,
    // the model governing the single nucleotide polymorphisms for this run.
    snp_model: SnpModel,
    // The model governing indels for this run.
    indel_model: IndelModel,
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
        &self, variant_type: VariantType, input_sequence: &Vec<Nuc>, rng: &mut NeatRng
    ) -> Vec<Nuc> {
        match variant_type {
            VariantType::SNP => {
                self.snp_model.generate_variant(
                    input_sequence, &self.transition_matrix, rng
                )
            },
            VariantType::Indel => {
                self.indel_model.generate_variant(
                    input_sequence, &self.transition_matrix, rng
                )
            },
        }
    }
}

// The following section are the models for each type of variant. In order to create the variant,
// we need to model its statistical property. NEAT included two types of variants: SNPs and Indels.
// Note that indels are actually two types, insertions and deletions, but they are usually classed
// together in the literature and are usually short. Insertions are often slips that cause sections
// to be duplicated, but NEAT made no attempt to distinguish between the types of insertions or
// deletions, since they are treated similarly by variant calling software.
pub trait Mutate {
    fn generate_variant(
        &self,
        input_sequence: &Vec<Nuc>,
        transition_matrix: &TransitionMatrix,
        rng: &mut NeatRng
    ) -> Vec<Nuc>;

}

#[derive(Debug)]
pub struct TransitionMatrix {
    // Nucleotide transition matrix. Rows represent the base we are mutating and the weights are
    // in the standard nucleotide order (in the same a, c, g, t order)
    //
    // The model is a 4x4 matrix with zeros along the diagonal because, e.g., A can't "mutate" to A.
    // The model is usually symmetric, but technically, the probability for A -> G could be
    // different from the probability for G -> A, but in practice, this seems to not be the case.
    a_weights: Vec<usize>,
    c_weights: Vec<usize>,
    g_weights: Vec<usize>,
    t_weights: Vec<usize>,
}

impl TransitionMatrix {
    pub fn new() -> Self {
        // Default transition matrix for mutations from the original NEAT 2.0
        Self {
            a_weights: vec![0, 15, 70, 15],
            c_weights: vec![15, 0, 15, 70],
            g_weights: vec![70, 15, 0, 15],
            t_weights: vec![15, 70, 15, 0],
        }
    }
    #[allow(dead_code)]
    // todo, once we have numbers we can implement this.
    pub fn from(weights: Vec<Vec<usize>>) -> Self {
        // Supply a vector of 4 vectors that define the mutation chance
        // from the given base to the other 4 bases.

        // First some safety checks. This should be a 4x4 matrix defining mutation from
        // ACGT (top -> down) to ACGT (left -> right)
        if weights.len() != 4 {
            panic!("Weights supplied to TransitionMatrix is wrong size");
        }
        for weight_vec in &weights {
            if weight_vec.len() != 4 {
                panic!("Weights supplied to TransitionMatrix is wrong size");
            }
        }
        Self {
            a_weights: weights[0].clone(),
            c_weights: weights[1].clone(),
            g_weights: weights[2].clone(),
            t_weights: weights[3].clone(),
        }
    }
}

impl Mutate for TransitionMatrix {
    fn generate_variant(
        &self, input_sequence: &Vec<Nuc>, transition_matrix: TransitionMatrix, rng: &mut NeatRng
    ) -> Vec<Nuc> {
        // This is a basic mutation function for starting us off
        // Pick the weights list for the base that was input
        // We will use this simple model for sequence errors ultimately.
        debug!("Generating basic SNP variant");
        if input_sequence.len() > 1 {
            panic!("Basic variants can only be one base long.")
        }
        let base = input_sequence[0];
        let weights: &Vec<usize> = match base {
            Nuc::A => &transition_matrix.a_weights,
            Nuc::C => &transition_matrix.c_weights,
            Nuc::G => &transition_matrix.g_weights,
            Nuc::T => &transition_matrix.t_weights,
            // return the N value for N with no further computation.
            Nuc::N => { return vec![Nuc::N]; },
        };
        // Now we create a distribution from the weights and sample our choices.
        let dist = WeightedIndex::new(weights).unwrap();
        match dist.sample(rng) {
            0 => vec![Nuc::A],
            1 => vec![Nuc::C],
            2 => vec![Nuc::G],
            3 => vec![Nuc::T],
            _ => vec![Nuc::N],
        }
    }
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
        input_sequence: &Vec<Nuc>,
        transition_matrix: &TransitionMatrix,
        rng: &mut NeatRng
    ) -> Vec<Nuc> {
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
        input_sequence: &Vec<Nuc>,
        transition_matrix: &TransitionMatrix,
        rng: &mut NeatRng
    ) -> Vec<Nuc> {
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

#[cfg(test)]
mod tests {
    use rand_core::SeedableRng;
    use utils::mutation_model::{MutationModel, TransitionMatrix};
    use utils::neat_rng::NeatRng;
    use utils::nucleotides::Nuc::*;
    #[test]
    fn test_transition_matrix_build() {
        let a_weights = vec![0, 20, 1, 20];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];

        let model = TransitionMatrix {
            a_weights: a_weights.clone(),
            c_weights: c_weights.clone(),
            g_weights: g_weights.clone(),
            t_weights: t_weights.clone(),
        };

        assert_eq!(model.a_weights, a_weights);
    }

    #[test]
    fn test_transition_matrix_from_weights() {
        let a_weights = vec![0, 20, 1, 20];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];
        let mut rng = NeatRng::seed_from_u64(0);
        let matrix = TransitionMatrix::from(
            vec![a_weights, c_weights, g_weights, t_weights]
        );
        let test_model = MutationModel::from_transition_matrix(matrix);
        // It actually mutates the base
        assert_ne!(test_model.mutate(A, &mut rng), A);
        assert_ne!(test_model.mutate(C, &mut rng), C);
        assert_ne!(test_model.mutate(G, &mut rng), G);
        assert_ne!(test_model.mutate(T, &mut rng), T);
        // It gives back N when you give it N
        assert_eq!(test_model.mutate(N, &mut rng), N);
    }
    #[test]
    #[should_panic]
    fn test_transition_matrix_too_many_vecs() {
        let a_weights = vec![0, 20, 1, 20];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];
        let u_weights = vec![20, 1, 20, 0];
        TransitionMatrix::from(vec![a_weights, c_weights, g_weights, t_weights, u_weights]);
    }

    #[test]
    #[should_panic]
    fn test_transition_matrix_too_many_bases() {
        let a_weights = vec![0, 20, 1, 20, 1];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];
        TransitionMatrix::from(vec![a_weights, c_weights, g_weights, t_weights]);
    }
}