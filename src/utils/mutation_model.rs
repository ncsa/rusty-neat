use std::collections::HashMap;
use rand::distributions::{Distribution, WeightedIndex};
use utils::neat_rng::NeatRng;

#[derive(Debug)]
pub enum Variant {
    Indel,
    SNP,
}
#[derive(Debug)]
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
    // The probability given that a variant has occurred at a location that it is one of the types
    // of possible variants occurring.
    variant_probs: HashMap<Variant, f64>,
    // The 4 x 4 matrix that shows the probability of one nucleotide transitioning to another.
    transition_matrix: TransitionMatrix,
    // the model governing the single nucleotide polymorphisms for this run.
    snp_model: SnpModel,
    // The model governing indels for this run.
    indel_model: IndelModel,
}
impl MutationModel {
    pub fn mutate(&self, base: u8, mut rng: &mut NeatRng) -> u8 {
        // This is a basic mutation function for starting us off
        //
        // the canonical choices for DNA, as defined in the Nucleotides module
        let choices: [u8; 4] = [0, 1, 2, 3];
        // Pick the weights list for the base that was input
        let weights: &Vec<usize> = match base {
            0 => &self.transition_matrix.a_weights,
            1 => &self.transition_matrix.c_weights,
            2 => &self.transition_matrix.g_weights,
            3 => &self.transition_matrix.t_weights,
            // anything else we skip the hard part and return the N value of 4
            _ => { return 4; },
        };
        // Now we create a distribution from the weights and sample our choices.
        let dist = WeightedIndex::new(weights).unwrap();
        choices[dist.sample(&mut rng)]
    }
}
#[derive(Debug)]
struct TransitionMatrix {
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
#[derive(Debug)]
struct SnpModel {
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
        // mutate with equal probability and mutate with the same probability (the default
        // tranisition matrix).
        let mut trinuc_transition_matrices = Vec::new();
        let mut trinuc_mutation_weights = Vec::new();
        for _ in 0..16 {
            trinuc_transition_matrices.push(TransitionMatrix::new());
            trinuc_mutation_weights.push(1);
        }

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
    use utils::mutation_model::TransitionMatrix;
    use utils::neat_rng::NeatRng;
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

        let str = format!("{:?}", model);
        let str_repr = String::from("NucModel { a: [0, 20, 1, 20], c: [20, 0, 1, 1], g: [1, 1, 0, 20], t: [20, 1, 20, 0] }");
        assert_eq!(str, str_repr);
        assert_eq!(model.a_weights, a_weights);
    }

    #[test]
    fn test_transition_matrix_from_weights() {
        let a_weights = vec![0, 20, 1, 20];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];
        let mut rng = NeatRng::seed_from_u64(0);
        let test_model = TransitionMatrix::from(
            vec![a_weights, c_weights, g_weights, t_weights]
        );
        // It actually mutates the base
        assert_ne!(test_model.choose_new_nuc(0, &mut rng), 0);
        assert_ne!(test_model.choose_new_nuc(1, &mut rng), 1);
        assert_ne!(test_model.choose_new_nuc(2, &mut rng), 2);
        assert_ne!(test_model.choose_new_nuc(3, &mut rng), 3);
        // It gives back N when you give it N
        assert_eq!(test_model.choose_new_nuc(4, &mut rng), 4);
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