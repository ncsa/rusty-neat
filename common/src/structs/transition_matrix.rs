use simple_rng::DiscreteDistribution;

#[derive(Debug, Clone)]
pub struct TransitionMatrix {
    // Nucleotide transition matrix. Rows represent the base we are mutating and the weights are
    // in the standard nucleotide order (in the same a, c, g, t order). This structure is
    // fundamental to the others and to the mutation model in general.
    //
    // The model is a 4x4 matrix with zeros along the diagonal because, e.g., A can't "mutate" to A.
    // The model is usually symmetric, but technically, the probability for A -> G could be
    // different from the probability for G -> A, but in practice, this seems to not be the case.
    pub(crate) a_dist: DiscreteDistribution,
    pub(crate) c_dist: DiscreteDistribution,
    pub(crate) g_dist: DiscreteDistribution,
    pub(crate) t_dist: DiscreteDistribution,
}

impl TransitionMatrix where {
    pub fn default() -> Self {
        // Default transition matrix for mutations from the original NEAT 2.0
        Self {
            a_dist: DiscreteDistribution::new(&vec![0.0, 15.0, 70.0, 15.0]),
            c_dist: DiscreteDistribution::new(&vec![15.0, 0.0, 15.0, 70.0]),
            g_dist: DiscreteDistribution::new(&vec![70.0, 15.0, 0.0, 15.0]),
            t_dist: DiscreteDistribution::new(&vec![15.0, 70.0, 15.0, 0.0]),
        }
    }

    pub fn from(
        a_weights: Vec<f64>,
        c_weights: Vec<f64>,
        g_weights: Vec<f64>,
        t_weights: Vec<f64>
    ) -> Self {
        let weights_test = vec![
            a_weights.clone(), c_weights.clone(), g_weights.clone(), t_weights.clone()
        ];
        for vector in weights_test.iter() {
            if vector.len() != 4 {
                panic!("Weights must be of length 4")
            }
        }
        Self {
            a_dist: DiscreteDistribution::new(&a_weights),
            c_dist: DiscreteDistribution::new(&c_weights),
            g_dist: DiscreteDistribution::new(&g_weights),
            t_dist: DiscreteDistribution::new(&t_weights),
        }
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_transition_matrix_build() {
        // todo fix test
        let a_weights = vec![0.0, 20.0, 1.0, 20.0];
        let c_weights = vec![20.0, 0.0, 1.0, 1.0];
        let g_weights = vec![1.0, 1.0, 0.0, 20.0];
        let t_weights = vec![20.0, 1.0, 20.0, 0.0];

        let model = TransitionMatrix::from(
            a_weights,
            c_weights,
            g_weights,
            t_weights,
        );

        println!("{:?}", model);
        // assert_eq!(model.a_dist, a_weights);
    }

    #[test]
    #[should_panic]
    fn test_transition_matrix_too_many_bases() {
        let a_weights = vec![0.0, 20.0, 1.0, 20.0, 1.0];
        let c_weights = vec![20.0, 0.0, 1.0, 1.0];
        let g_weights = vec![1.0, 1.0, 0.0, 20.0];
        let t_weights = vec![20.0, 1.0, 20.0, 0.0];
        TransitionMatrix::from(a_weights, c_weights, g_weights, t_weights);
    }
}
