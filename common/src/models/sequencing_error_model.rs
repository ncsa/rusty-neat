use log::debug;
use std::io;
use serde::{Deserialize, Serialize};
use thiserror::Error;
use crate::{
    models::lib::{model_reader, model_writer}, 
    structs::{distributions::{DiscreteDistribution, DistributionErrors}, 
    nucleotides::{Nucleotide, ALLOWED_USIZE}, 
    transition_matrix::{TransitionMatrix, TransitionMatrixError}}
};
use simple_rng::{NeatRng, NeatRngError};

#[derive(Error, Debug)]
pub enum SeqModelError{
    #[error("Error creating sequencing error model")]
    ModelCreationError,
    #[error("Error creating transition matrix: {0}")]
    TransMatrixError(TransitionMatrixError),
    #[error("Error sampling distribution: {0}")]
    DistributionError(DistributionErrors),
    #[error("Error with rng: {0}")]
    RngError(NeatRngError),
    #[error("No RNG supplied for this model.")]
    MissingRngError,
    #[error("Sequencing Error model return an IO error: {0}")]
    IoError(io::Error),
}

impl From<io::Error> for SeqModelError {
    fn from(error: io::Error) -> Self {
        SeqModelError::IoError(error)
    }
}

impl From<DistributionErrors> for SeqModelError {
    fn from(error: DistributionErrors) -> Self {
        SeqModelError::DistributionError(error)
    }
}

impl From<NeatRngError> for SeqModelError {
    fn from(error: NeatRngError) -> Self {
        SeqModelError::RngError(error)
    }
}

impl From<TransitionMatrixError> for SeqModelError {
    fn from(error: TransitionMatrixError) -> Self {
        SeqModelError::TransMatrixError(error)
    }
}

pub enum SequencingErrorType {
    SnpError(Nucleotide),
    InsertionError(Vec<Nucleotide>),
    DeletionError(usize),
}

#[derive(Clone, Serialize, Deserialize)]
pub struct SequencingErrorModel {
    // Neat only dealt with 2 types of sequencing errors: snps and small indels.
    // We will retain that idea and assume it is accurate.
    error_rate: f64,
    del_length_distribution: DiscreteDistribution,
    ins_length_distribution: DiscreteDistribution,
    indel_probability: f64,
    insertion_bias: DiscreteDistribution,
    transition_distros: TransitionMatrix,
}

impl SequencingErrorModel {
    pub fn default() -> Result<Self, SeqModelError> {
        // This is the default sequencing error model employed by NEAT2
        // Note that this was originally in a file, and we could have done it the way we did
        // the other defaults, but it was so small, I just included it in full here.
        let default_transition_distros = TransitionMatrix::from(
            vec![0.0, 0.4918, 0.3377, 0.1705],
            vec![0.5238, 0.0, 0.2661, 0.2101],
            vec![0.3754, 0.2355, 0.0, 0.389],
            vec![0.2505, 0.2552, 0.4942, 0.0],
        )?;
        let default_error_rate = 0.006638164688495656;
        let default_lengths = vec![1, 2];
        let default_ins_distr = DiscreteDistribution::new(
            &vec![0.999, 0.001],
            &default_lengths,
        )?;
        let default_del_distr = default_ins_distr.clone();
        let default_indel_probability = 0.4;
        // default is no bias
        let default_insertion_bias = DiscreteDistribution::new(
            &vec![1.0, 1.0, 1.0, 1.0],
            &Vec::from(ALLOWED_USIZE.clone()),
        )?;

        Ok(SequencingErrorModel {
            error_rate: default_error_rate,
            del_length_distribution: default_del_distr,
            ins_length_distribution: default_ins_distr,
            indel_probability: default_indel_probability,
            insertion_bias: default_insertion_bias,
            transition_distros: default_transition_distros,
        })
    }

    pub fn from_file(filename: &str) -> Result<Self, SeqModelError> {
        let data: SequencingErrorModel = model_reader(filename).unwrap();
        Ok(data)
    }

    pub fn write_model(&self, filename: &str) -> Result<(), SeqModelError> {
        model_writer(self, filename)?;
        Ok(())
    }

    pub fn generate_sequencing_error(&self, reference: Nucleotide, rng: &mut NeatRng) -> Result<SequencingErrorType, SeqModelError> {
        // This method picks an error type and determines any additional data needed
        // for the current error, based on the statistical model
        if rng.random()? < self.indel_probability {
            // Indel error
            Ok(self.generate_indel_error(rng)?)
        } else {
            // SNP error
            Ok(SequencingErrorType::SnpError(self.generate_snp_error(reference, rng.random()?)?))
        }
    }

    pub fn convert_score(&self, score: usize) -> Result<f64, SeqModelError> {
        // Takes a quality score, converts it to a probability of error, and returns the result
        let score = score as f64;
        Ok(10.0_f64.powf(-score/10.0))
    }

    fn generate_snp_error(&self, reference: Nucleotide, rand: f64) -> Result<Nucleotide, SeqModelError> {
        // This is a basic mutation function for starting us off
        // Pick the weights list for the base that was input
        // We will use this simple model for sequence errors ultimately.
        debug!("Generating basic SNP variant");
        let distro = &self.transition_distros[&reference];
        // Now we create a distribution from the weights and sample our choices.
        // We have constructed things such that this will return a valid u8. But 
        // to be extra safe, we could mod by 4 and then convert
        Ok(Nucleotide::from(distro.sample(rand)?))
    }

    fn generate_indel_error(&self, rng: &mut NeatRng) -> Result<SequencingErrorType, SeqModelError> {
        // Returns either an insertion (option 1) or a deletion (option 2) depending on a random selection from a list of potential
        // error lengths (-2..2). This makes an insertion of up to 2 bases as likely as a random deletion of up to 2 bases.
        if rng.random()? < 0.5 {
            // We assume fifty-fifty chance of insertion v deletion
            // insertion
            let mut sequence = Vec::new();
            let length = self.ins_length_distribution.sample(rng.random()?)?;
            for _ in 0..length {
                // We could mod this value by 4 to ensure it is a valid base. Or create a data structure.
                sequence.push(Nucleotide::from(self.insertion_bias.sample(rng.random()?)?))
            }
            // Insertion of sequence
            Ok(SequencingErrorType::InsertionError(sequence))
        } else {
            // Deletion
            let length = self.del_length_distribution.sample(rng.random()?)?;
            Ok(SequencingErrorType::DeletionError(length))
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