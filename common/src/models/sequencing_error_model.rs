use log::debug;
use rand::distributions::WeightedIndex;
use rand::prelude::Distribution;

use structs::nucleotides::Nuc;
use structs::sequencing_errors::SequencingError;
use structs::sequencing_errors::{IndelErr, SnpErr};
use structs::transition_matrix::TransitionMatrix;
use rand_chacha::ChaCha20Rng;
use structs::sequencing_errors::SequencingError::{IndelError, SNPError};

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
    pub fn generate_snp_error(&self, base: Nuc, rng: &mut ChaCha20Rng) -> SequencingError {
        // This is a basic mutation function for starting us off
        // Pick the weights list for the base that was input
        // We will use this simple model for sequence errors ultimately.
        debug!("Generating basic SNP variant");
        let weights: &Vec<usize> = match base {
            Nuc::A => &self.transition_matrix.a_weights,
            Nuc::C => &self.transition_matrix.c_weights,
            Nuc::G => &self.transition_matrix.g_weights,
            Nuc::T => &self.transition_matrix.t_weights,
            // return the N value for N with no further computation.
            Nuc::N => {
                return SNPError(SnpErr::new(Nuc::N));
            }
        };
        // Now we create a distribution from the weights and sample our choices.
        let dist = WeightedIndex::new(weights).unwrap();
        match dist.sample(rng) {
            0 => SNPError(SnpErr::new(Nuc::A)),
            1 => SNPError(SnpErr::new(Nuc::C)),
            2 => SNPError(SnpErr::new(Nuc::G)),
            3 => SNPError(SnpErr::new(Nuc::T)),
            // technically, the following part isn't reachable because of how we have constructed
            // things. Our weights vector has a max length of 4, but Rust doesn't know that.
            _ => SNPError(SnpErr::new(Nuc::N)),
        }
    }

    pub fn generate_indel_error(&self, rng: &mut ChaCha20Rng) -> SequencingError {
        let dist = WeightedIndex::new(&self.length_weights).unwrap();
        let length = self.lengths[dist.sample(rng)];
        return if length > 0 {
            // insertion
            let mut sequence = Vec::new();
            let dist = WeightedIndex::new(&self.insertion_bias).unwrap();
            for _ in 0..length {
                sequence.push(match dist.sample(rng) {
                    0 => Nuc::A,
                    1 => Nuc::C,
                    2 => Nuc::G,
                    3 => Nuc::T,
                    _ => Nuc::N,
                })
            }
            // Insertion of sequence
            IndelError(IndelErr::new(length, Some(sequence)))
        } else {
            // Deletion of length bases
            IndelError(IndelErr::new(length, None))
        };
    }
}
