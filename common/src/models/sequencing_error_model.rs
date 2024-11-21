use log::debug;

use crate::structs::sequencing_errors::SequencingError;
use crate::structs::sequencing_errors::{IndelErr, SnpErr};
use crate::structs::transition_matrix::TransitionMatrix;
use crate::structs::sequencing_errors::SequencingError::{IndelError, SNPError};
use simple_rng::{NeatRng, DiscreteDistribution};

#[derive(Clone)]
pub struct SequencingErrorModel {
    // Neat only dealt with 2 types of sequencing errors: snps and small indels.
    // We will retain that idea and assume it is accurate.
    error_rate: f64,
    lengths: [i8; 4],
    length_distr: DiscreteDistribution,
    insertion_probability: f64,
    insertion_bias: DiscreteDistribution,
    transition_distros: TransitionMatrix,
}

impl SequencingErrorModel {
    pub fn new() -> Self {
        // This is the default sequencing error model employed by NEAT2
        let default_transition_distros = TransitionMatrix::from(
            vec![0.0, 0.4918, 0.3377, 0.1705],
            vec![0.5238, 0.0, 0.2661, 0.2101],
            vec![0.3754, 0.2355, 0.0, 0.389],
            vec![0.2505, 0.2552, 0.4942, 0.0],
        );
        let default_error_rate = 0.01;
        let default_lengths = [-2, -1, 1, 2];
        let default_length_distr = DiscreteDistribution::new(&vec![0.001, 0.999, 0.999, 0.001]);
        let default_insertion_probability = 0.4;
        // default is no bias
        let default_insertion_bias = DiscreteDistribution::new(&vec![1.0, 1.0, 1.0, 1.0]);

        SequencingErrorModel {
            error_rate: default_error_rate,
            lengths: default_lengths,
            length_distr: default_length_distr,
            insertion_probability: default_insertion_probability,
            insertion_bias: default_insertion_bias,
            transition_distros: default_transition_distros,
        }
    }
    pub fn generate_snp_error(&self, base: u8, rng: &mut NeatRng) -> SequencingError {
        // This is a basic mutation function for starting us off
        // Pick the weights list for the base that was input
        // We will use this simple model for sequence errors ultimately.
        debug!("Generating basic SNP variant");
        let weights: &DiscreteDistribution = match base {
            0 => &self.transition_distros.a_dist,
            1 => &self.transition_distros.c_dist,
            2 => &self.transition_distros.g_dist,
            3 => &self.transition_distros.t_dist,
            _ => panic!("Trying to mutate an unknown bases!")
        };
        // Now we create a distribution from the weights and sample our choices.
        match weights.sample(rng) {
            0 => SNPError(SnpErr::new(0)),
            1 => SNPError(SnpErr::new(1)),
            2 => SNPError(SnpErr::new(2)),
            3 => SNPError(SnpErr::new(3)),
            // technically, the following part isn't reachable because of how we have constructed
            // things. Our weights vector has a max length of 4, but Rust doesn't know that.
            _ => panic!("Invalid index selected!"),
        }
    }

    pub fn generate_indel_error(&self, rng: &mut NeatRng) -> SequencingError {
        let length = self.lengths[self.length_distr.sample(rng)];
        if length > 0 {
            // insertion
            let mut sequence = Vec::new();
            for _ in 0..length {
                sequence.push(self.insertion_bias.sample(rng) as u8)
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
        // todo needs tests
        println!("TODO");
        assert_eq!(1, 1)
    }
}