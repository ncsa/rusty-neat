//! The transition matrix extends the notion of a distribution. Specifically, it denotes a transition
//! probability from one base to another. The trinucleotide SNP model requires 4 layers of these transition
//! matrices. Indexing is now implemented for the Tranisition matrix, which should make the code more straightforward.
use crate::structs::{
    distributions::{DiscreteDistribution, DistributionErrors},
    nucleotides::{ALLOWED_NUCS, Nucleotide},
};
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::ops::Index;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum TransitionMatrixError {
    #[error("Transition matrix reported a distribution error: {0}")]
    DistributionError(DistributionErrors),
    #[error("Input weights and lengths were of unequal length")]
    UnequalWeightsError,
    #[error("I/O error reading transition matrix file: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Transition matrix file has {0} data rows (expected 4)")]
    InvalidRowCount(usize),
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
    pub a: DiscreteDistribution<Nucleotide>,
    pub c: DiscreteDistribution<Nucleotide>,
    pub g: DiscreteDistribution<Nucleotide>,
    pub t: DiscreteDistribution<Nucleotide>,
}

impl Index<usize> for TransitionMatrix {
    type Output = DiscreteDistribution<Nucleotide>;
    fn index(&self, i: usize) -> &DiscreteDistribution<Nucleotide> {
        match i {
            0 => &self.a,
            1 => &self.c,
            2 => &self.g,
            3 => &self.t,
            _ => panic!("index out of range: {} with length 4", i),
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
            _ => panic!("Nucleotide not found: {:?}", n),
        }
    }
}

impl TransitionMatrix {
    pub fn default() -> Result<Self, TransitionMatrixError> {
        // Default transition matrix for mutations from the original NEAT 2.0
        Ok(Self {
            a: DiscreteDistribution::new(
                &vec![
                    0.0,
                    0.16952785544157142,
                    0.6878218371885525,
                    0.1426503073698761,
                ],
                &Vec::from(ALLOWED_NUCS),
            )?,
            c: DiscreteDistribution::new(
                &vec![
                    0.1615595177675533,
                    0.0,
                    0.1664600510853558,
                    0.6719804311470908,
                ],
                &Vec::from(ALLOWED_NUCS),
            )?,
            g: DiscreteDistribution::new(
                &vec![
                    0.6704692217289652,
                    0.16706325641014994,
                    0.0,
                    0.1624675218608849,
                ],
                &Vec::from(ALLOWED_NUCS),
            )?,
            t: DiscreteDistribution::new(
                &vec![
                    0.14281397280584493,
                    0.6885459433549382,
                    0.16864008383921686,
                    0.0,
                ],
                &Vec::from(ALLOWED_NUCS),
            )?,
        })
    }

    /// Load a 4×4 SNP transition matrix from a whitespace-delimited TSV file.
    ///
    /// Rows and columns correspond to A/C/G/T (from-base and to-base).
    /// An optional header row is skipped if its first token is non-numeric.
    /// Diagonal values are zeroed so self-transitions are impossible.
    pub fn from_tsv(path: &PathBuf) -> Result<Self, TransitionMatrixError> {
        let content = std::fs::read_to_string(path)?;
        let mut rows: Vec<[f64; 4]> = Vec::new();

        for line in content.lines() {
            let tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.is_empty() {
                continue;
            }
            if rows.is_empty() && tokens[0].parse::<f64>().is_err() {
                continue; // skip header
            }
            let vals: Vec<f64> = tokens.iter().filter_map(|s| s.parse().ok()).collect();
            if vals.len() >= 4 {
                rows.push([vals[0], vals[1], vals[2], vals[3]]);
            }
        }

        if rows.len() < 4 {
            return Err(TransitionMatrixError::InvalidRowCount(rows.len()));
        }

        for (i, row) in rows.iter_mut().enumerate() {
            row[i] = 0.0; // zero diagonal so self-transitions are impossible
        }

        Self::from(rows[0], rows[1], rows[2], rows[3])
    }

    pub fn from(
        a_weights: [f64; 4],
        c_weights: [f64; 4],
        g_weights: [f64; 4],
        t_weights: [f64; 4],
    ) -> Result<Self, TransitionMatrixError> {
        let weights_test = vec![
            a_weights.clone(),
            c_weights.clone(),
            g_weights.clone(),
            t_weights.clone(),
        ];
        for vector in weights_test.iter() {
            if vector.len() != 4 {
                return Err(TransitionMatrixError::UnequalWeightsError);
            }
        }
        Ok(Self {
            a: DiscreteDistribution::new(&a_weights.to_vec(), &Vec::from(ALLOWED_NUCS))?,
            c: DiscreteDistribution::new(&c_weights.to_vec(), &Vec::from(ALLOWED_NUCS))?,
            g: DiscreteDistribution::new(&g_weights.to_vec(), &Vec::from(ALLOWED_NUCS))?,
            t: DiscreteDistribution::new(&t_weights.to_vec(), &Vec::from(ALLOWED_NUCS))?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::NeatRng;

    #[test]
    fn test_transition_matrix_build() {
        let a_weights = [0.0, 20.0, 1.0, 20.0];
        let c_weights = [20.0, 0.0, 1.0, 1.0];
        let g_weights = [1.0, 1.0, 0.0, 20.0];
        let t_weights = [20.0, 1.0, 20.0, 0.0];
        let model = TransitionMatrix::from(a_weights, c_weights, g_weights, t_weights).unwrap();
        // Index by usize and Nucleotide must reference the same distribution
        assert_eq!(
            model[0].values().unwrap(),
            model[&Nucleotide::A].values().unwrap()
        );
        assert_eq!(
            model[1].values().unwrap(),
            model[&Nucleotide::C].values().unwrap()
        );
        assert_eq!(
            model[2].values().unwrap(),
            model[&Nucleotide::G].values().unwrap()
        );
        assert_eq!(
            model[3].values().unwrap(),
            model[&Nucleotide::T].values().unwrap()
        );
        // Each row's values should be the four ACGT nucleotides
        assert_eq!(
            model[&Nucleotide::A].values().unwrap(),
            Vec::from(ALLOWED_NUCS)
        );
    }

    #[test]
    fn test_from_tsv_valid() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("matrix.tsv");
        std::fs::write(
            &tsv,
            "A\tC\tG\tT\n\
             0.0\t0.5\t0.3\t0.2\n\
             0.5\t0.0\t0.3\t0.2\n\
             0.4\t0.3\t0.0\t0.3\n\
             0.3\t0.3\t0.4\t0.0\n",
        )
        .unwrap();
        let tm = TransitionMatrix::from_tsv(&tsv).unwrap();
        // Check each row can sample one of the four ACGT nucleotides
        let mut rng = NeatRng::new_from_seed(&vec!["t".to_string()]).unwrap();
        for nuc in ALLOWED_NUCS {
            let sample = tm[&nuc].sample(rng.random().unwrap()).unwrap();
            assert!(ALLOWED_NUCS.contains(&sample));
            // Self-transition must be impossible: diagonal was zeroed
            assert_ne!(sample, nuc);
        }
    }

    #[test]
    fn test_from_tsv_too_few_rows_errors() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("short.tsv");
        std::fs::write(&tsv, "0.0\t0.5\t0.3\t0.2\n0.5\t0.0\t0.3\t0.2\n").unwrap();
        assert!(matches!(
            TransitionMatrix::from_tsv(&tsv),
            Err(TransitionMatrixError::InvalidRowCount(2))
        ));
    }

    #[test]
    fn test_transition_matrix_default() {
        let model = TransitionMatrix::default().unwrap();
        // Spot-check: sampling from A row should return one of [A, C, G, T]
        let mut rng = NeatRng::new_from_seed(&vec!["seed".to_string()]).unwrap();
        let sample = model[&Nucleotide::A].sample(rng.random().unwrap()).unwrap();
        assert!(ALLOWED_NUCS.contains(&sample));
    }
}
