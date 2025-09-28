//! The transition matrix extends the notion of a distribution. Specifically, it denotes a transition
//! probability from one base to another. The trinucleotide SNP model requires 4 layers of these transition
//! matrices. Indexing is now implemented for the Tranisition matrix, which should make the code more straightforward.
use crate::structs::{
    distributions::{DiscreteDistribution, DistributionErrors}, 
    nucleotides::{Nucleotide, ALLOWED_NUCS}
};
use thiserror::Error;
use std::ops::Index;
use std::fmt::Debug;
use serde::{Serialize, Deserialize};

#[derive(Debug, Error)]
pub enum TransitionMatrixError {
    #[error("Transition matrix reported a distribution error: {0}")]
    DistributionError(DistributionErrors),
    #[error("Input weights and lengths were of unequal length")]
    UnequalWeightsError,
}

impl From<DistributionErrors> for TransitionMatrixError {
    fn from(error: DistributionErrors) -> Self {
        TransitionMatrixError::DistributionError(error)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransitionMatrix {
    // Nucleotide transition matrix. Rows represent the base we are mutating and the weights are
    // in the standard nucleotide order (in the same a, c, g, t order). This structure is
    // fundamental to the others and to the mutation model in general.
    //
    // This defines a transition from one Nucleotide to another.
    a: DiscreteDistribution<Nucleotide>,
    c: DiscreteDistribution<Nucleotide>,
    g: DiscreteDistribution<Nucleotide>,
    t: DiscreteDistribution<Nucleotide>,
}

impl Index<usize> for TransitionMatrix {
    type Output = DiscreteDistribution<Nucleotide>;
    fn index(&self, i: usize) -> &DiscreteDistribution<Nucleotide> {
        match i {
            0 => &self.a,
            1 => &self.c,
            2 => &self.g,
            3 => &self.t,
            _ => panic!("index out of range: {} with length 4", i)
        }
    }
}

impl Index<&Nucleotide> for TransitionMatrix {
    type Output = DiscreteDistribution<Nucleotide>;
    fn index(&self, n: &Nucleotide) -> &DiscreteDistribution<Nucleotide> {
        match n {
            Nucleotide::A => &self.a,
            Nucleotide::C => &self.c,
            Nucleotide::G => &self.g,
            Nucleotide::T => &self.t,
            _ => panic!("Nucleotide not found: {:?}", n)
        }
    }
}

impl TransitionMatrix where {
    pub fn default() -> Result<Self, TransitionMatrixError> {
        // Default transition matrix for mutations from the original NEAT 2.0
        Ok(Self {
            a: DiscreteDistribution::new(
                &vec![0.0, 0.16952785544157142, 0.6878218371885525,  0.1426503073698761],
                &Vec::from(ALLOWED_NUCS),
            )?,
            c: DiscreteDistribution::new(
                &vec![0.1615595177675533, 0.0, 0.1664600510853558, 0.6719804311470908],
                &Vec::from(ALLOWED_NUCS),
            )?,
            g: DiscreteDistribution::new(
                &vec![0.6704692217289652, 0.16706325641014994, 0.0, 0.1624675218608849],
                &Vec::from(ALLOWED_NUCS),
            )?,
            t: DiscreteDistribution::new(
                &vec![0.14281397280584493, 0.6885459433549382, 0.16864008383921686, 0.0],
                &Vec::from(ALLOWED_NUCS),
            )?,
        })
    }

    pub fn from(
        a_weights: [f64; 4],
        c_weights: [f64; 4],
        g_weights: [f64; 4],
        t_weights: [f64; 4],
    ) -> Result<Self, TransitionMatrixError> {
        let weights_test = vec![
            a_weights.clone(), c_weights.clone(), g_weights.clone(), t_weights.clone()
        ];
        for vector in weights_test.iter() {
            if vector.len() != 4 {
                return Err(TransitionMatrixError::UnequalWeightsError)
            }
        }
        Ok(Self {
            a: DiscreteDistribution::new(
                &a_weights.to_vec(),
                &Vec::from(ALLOWED_NUCS),
            )?,
            c: DiscreteDistribution::new(
                &c_weights.to_vec(),
                &Vec::from(ALLOWED_NUCS),
            )?,
            g: DiscreteDistribution::new(
                &g_weights.to_vec(),
                &Vec::from(ALLOWED_NUCS),
            )?,
            t: DiscreteDistribution::new(
                &t_weights.to_vec(),
                &Vec::from(ALLOWED_NUCS),
            )?,
        })
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_transition_matrix_build() {
        // todo fix test
        let a_weights = [0.0, 20.0, 1.0, 20.0];
        let c_weights = [20.0, 0.0, 1.0, 1.0];
        let g_weights = [1.0, 1.0, 0.0, 20.0];
        let t_weights = [20.0, 1.0, 20.0, 0.0];

        let model = TransitionMatrix::from(
            a_weights,
            c_weights,
            g_weights,
            t_weights,
        );

        println!("{:?}", model);
        // assert_eq!(model.a_dist, a_weights);
    }
}
