#[derive(Debug, Clone)]
pub struct TransitionMatrix {
    // Nucleotide transition matrix. Rows represent the base we are mutating and the weights are
    // in the standard nucleotide order (in the same a, c, g, t order). This structure is
    // fundamental to the others and to the mutation model in general.
    //
    // The model is a 4x4 matrix with zeros along the diagonal because, e.g., A can't "mutate" to A.
    // The model is usually symmetric, but technically, the probability for A -> G could be
    // different from the probability for G -> A, but in practice, this seems to not be the case.
    pub(crate) a_weights: Vec<u32>,
    pub(crate) c_weights: Vec<u32>,
    pub(crate) g_weights: Vec<u32>,
    pub(crate) t_weights: Vec<u32>,
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

    pub fn from(weights: Vec<Vec<u32>>) -> Self {
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

#[cfg(test)]
mod tests {
    use super::*;

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
