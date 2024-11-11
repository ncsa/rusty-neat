use log::debug;

use structs::sequencing_errors::SequencingError;
use structs::sequencing_errors::{IndelErr, SnpErr};
use structs::transition_matrix::TransitionMatrix;
use structs::sequencing_errors::SequencingError::{IndelError, SNPError};
use simple_rng::{Rng, DiscreteDistribution};

#[derive(Clone)]
pub struct SequencingErrorModel {
    // Neat only dealt with 2 types of sequencing errors: snps and small indels.
    // We will retain that idea and assume it is accurate.
    error_rate: f64,
    lengths: [i8; 4],
    length_weights: [u32; 4],
    insertion_probability: f64,
    insertion_bias: [u32; 4],
    transition_matrix: TransitionMatrix,
}

impl SequencingErrorModel {
    pub fn new() -> Self {
        // This is the default sequencing error model employed by NEAT2, scaled to be weights
        // instead of percentages (each item multiplied by a factor of 10/(min(row)) and rounded)
        // Original:
        // [[0.0, 0.4918, 0.3377, 0.1705],
        //   [0.5238, 0.0, 0.2661, 0.2101],
        //   [0.3754, 0.2355, 0.0, 0.389],
        //   [0.2505, 0.2552, 0.4942, 0.0]]
        let default_transition_matrix = TransitionMatrix::from(vec![
            vec![0, 29, 20, 10],
            vec![25, 0, 13, 10],
            vec![15, 10, 0, 17],
            vec![10, 10, 20, 0],
        ]);
        let default_error_rate = 0.01;
        let default_lengths = [-2, -1, 1, 2];
        let default_length_weights = [1, 999, 999, 1];
        let default_insertion_probability = 0.4;
        // default is no bias
        let default_insertion_bias = [1, 1, 1, 1];

        SequencingErrorModel {
            error_rate: default_error_rate,
            lengths: default_lengths,
            length_weights: default_length_weights,
            insertion_probability: default_insertion_probability,
            insertion_bias: default_insertion_bias,
            transition_matrix: default_transition_matrix,
        }
    }
    pub fn generate_snp_error(&self, base: u8, rng: &mut Rng) -> SequencingError {
        // This is a basic mutation function for starting us off
        // Pick the weights list for the base that was input
        // We will use this simple model for sequence errors ultimately.
        debug!("Generating basic SNP variant");
        let weights: &Vec<u32> = match base {
            0 => &self.transition_matrix.a_weights,
            1 => &self.transition_matrix.c_weights,
            2 => &self.transition_matrix.g_weights,
            3 => &self.transition_matrix.t_weights,
            _ => panic!("Ve should not be trying to mutated unknown bases!")
        };
        // Now we create a distribution from the weights and sample our choices.
        let dist = DiscreteDistribution::new(weights);
        match dist.sample(rng) {
            0 => SNPError(SnpErr::new(0)),
            1 => SNPError(SnpErr::new(1)),
            2 => SNPError(SnpErr::new(2)),
            3 => SNPError(SnpErr::new(3)),
            // technically, the following part isn't reachable because of how we have constructed
            // things. Our weights vector has a max length of 4, but Rust doesn't know that.
            _ => panic!("Invalid index selected!"),
        }
    }

    pub fn generate_indel_error(&self, rng: &mut Rng) -> SequencingError {
        let dist = DiscreteDistribution::new(&Vec::from(self.length_weights));
        let length = self.lengths[dist.sample(rng)];
        if length > 0 {
            // insertion
            let mut sequence = Vec::new();
            let dist = DiscreteDistribution::new(&Vec::from(self.insertion_bias));
            for _ in 0..length {
                sequence.push(dist.sample(rng) as u8)
            }
            // Insertion of sequence
            IndelError(IndelErr::new(length, Some(sequence)))
        } else {
            // Deletion of length bases
            IndelError(IndelErr::new(length, None))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_sequencing_error_model() {
        println!("TODO");
        assert_eq!(1, 1)
    }
}