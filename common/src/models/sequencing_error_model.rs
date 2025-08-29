use log::debug;
use thiserror::Error;
use crate::structs::transition_matrix::{TransitionMatrix, TransitionMatrixError};
use simple_rng::{DiscreteDistribution, NeatRng, NeatRngError};

#[derive(Error, Debug)]
pub enum SeqModelError{
    #[error("Error creating sequencing error model")]
    ModelCreationError,
    #[error("Error creating transition matrix: {0}")]
    TransMatrixError(TransitionMatrixError),
    #[error("Error sampling rng distribution: {0}")]
    DistributionError(NeatRngError),
    #[error("No RNG supplied for this model.")]
    MissingRngError,
}

impl From<NeatRngError> for SeqModelError {
    fn from(error: NeatRngError) -> Self {
        SeqModelError::DistributionError(error)
    }
}

impl From<TransitionMatrixError> for SeqModelError {
    fn from(error: TransitionMatrixError) -> Self {
        SeqModelError::TransMatrixError(error)
    }
}

pub enum SequencingErrorType {
    SnpError(u8),
    InsertionError(Vec<u8>),
    DeletionError(i8),
}

#[derive(Clone)]
pub struct SequencingErrorModel {
    // Neat only dealt with 2 types of sequencing errors: snps and small indels.
    // We will retain that idea and assume it is accurate.
    error_rate: f64,
    lengths: [i8; 4],
    length_distr: DiscreteDistribution,
    indel_probability: f64,
    insertion_bias: DiscreteDistribution,
    transition_distros: TransitionMatrix,
    rng: Option<NeatRng>,
}

impl SequencingErrorModel {
    pub fn default(rng: NeatRng) -> Result<Self, SeqModelError> {
        // This is the default sequencing error model employed by NEAT2
        let default_transition_distros = TransitionMatrix::from(
            vec![0.0, 0.4918, 0.3377, 0.1705],
            vec![0.5238, 0.0, 0.2661, 0.2101],
            vec![0.3754, 0.2355, 0.0, 0.389],
            vec![0.2505, 0.2552, 0.4942, 0.0],
        )?;
        let default_error_rate = 0.006638164688495656;
        let default_lengths = [-2, -1, 1, 2];
        let default_length_distr = DiscreteDistribution::new(&vec![0.001, 0.999, 0.999, 0.001])?;
        let default_indel_probability = 0.4;
        // default is no bias
        let default_insertion_bias = DiscreteDistribution::new(&vec![1.0, 1.0, 1.0, 1.0])?;

        Ok(SequencingErrorModel {
            error_rate: default_error_rate,
            lengths: default_lengths,
            length_distr: default_length_distr,
            indel_probability: default_indel_probability,
            insertion_bias: default_insertion_bias,
            transition_distros: default_transition_distros,
            rng: Some(rng),
        })
    }

    pub fn is_error(&self, score: usize) -> Result<bool, SeqModelError> {
        // Takes a quality score, converts it to a probability of error, and returns bool
        // true = this base is an error
        // false = this base is the same as reference
        let score = score as f64;
        match self.rng {
            Some(mut r) => Ok((&mut r).gen_bool(10.0_f64.powf(-score/10.0))?),
            None => return Err(SeqModelError::MissingRngError),
        }
    }

    pub fn generate_snp_error(&self, reference: u8) -> Result<u8, SeqModelError> {
        // This is a basic mutation function for starting us off
        // Pick the weights list for the base that was input
        // We will use this simple model for sequence errors ultimately.
        debug!("Generating basic SNP variant");
        let weights: &DiscreteDistribution = match reference {
            0 => &self.transition_distros.a_dist,
            1 => &self.transition_distros.c_dist,
            2 => &self.transition_distros.g_dist,
            3 => &self.transition_distros.t_dist,
            _ => panic!("Trying to mutate an unknown bases!")
        };
        // Now we create a distribution from the weights and sample our choices.
        // We have constructed things such that this will return a valid u8. But 
        // to be extra safe, we could mod by 4 and then convert
        match self.rng {
            Some(mut r) => Ok(weights.sample(&mut r)? as u8),
            None => return Err(SeqModelError::MissingRngError),
        }
    }

    pub fn generate_indel_error(&self) -> Result<SequencingErrorType, SeqModelError> {
        // Returns either an insertion (option 1) or a deletion (option 2) depending on a random selection from a list of potential
        // error lengths (-2..2). This makes an insertion of up to 2 bases as likely as a random deletion of up to 2 bases.
        match self.rng {
            Some(mut r) => {
                let index = self.length_distr.sample(&mut r)?;
                let length = self.lengths[index];
                match length {
                    1.. => {
                        // insertion
                        let mut sequence = Vec::new();
                        for _ in 0..length {
                            // We could mod this value by 4 to ensure it is a valid base. Or create a data structure.
                            sequence.push(self.insertion_bias.sample(&mut r)? as u8)
                        }
                        // Insertion of sequence
                        Ok(SequencingErrorType::InsertionError(sequence))
                    }

                    _ => {
                        // Deletion of length bases
                        // We'll return this as an error, and let the read generator interpret
                        Ok(SequencingErrorType::DeletionError(length))
                    }
                }
            },
            None => return Err(SeqModelError::MissingRngError),
        }
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_sequencing_error_model() {
        // todo needs tests
        println!("TODO");
        assert_eq!(1, 1)
    }
}