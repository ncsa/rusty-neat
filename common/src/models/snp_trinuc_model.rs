//! The following section are the models for each type of variant. In order to create the variant,
//! we need to model its statistical property. NEAT included two types of variants: SNPs and Indels.
//! Note that indels are actually two types, insertions and deletions, but they are usually classed
//! together in the literature and are usually short. Insertions are often slips that cause sections
//! to be duplicated, but NEAT made no attempt to distinguish between the types of insertions or
//! deletions, since they are treated similarly by variant calling software.
use log::error;
use vectorize;
use serde;
use serde::{Deserialize, Serialize};
use thiserror::Error;
use std::collections::HashMap;
use std::slice::Iter;
use simple_rng::NeatRngError;
use crate::models::lib::{model_gzp_reader, model_writer};
use crate::structs::transition_matrix::{TransitionMatrix, TransitionMatrixError};
use crate::structs::distributions::{DiscreteDistribution, DistributionErrors};

#[derive(Error, Debug)]
pub enum SnpTrinucError {
    #[error("SNP variant model reported an RNG error: {0}")]
    RngError(NeatRngError),
    #[error("SNP variant model reported an error from Distributions: {0}")]
    DistributionError(DistributionErrors),
    #[error("SNP variant model reported an error with transition matrices: {0}")]
    TransMatrixError(TransitionMatrixError),
    #[error("SNP variant model reported an error generating a SNP.")]
    GenerateSnpError,
}

impl From<DistributionErrors> for SnpTrinucError {
    fn from(error: DistributionErrors) -> Self {
        SnpTrinucError::DistributionError(error)
    }
}

impl From<NeatRngError> for SnpTrinucError {
    fn from(error: NeatRngError) -> Self {
        SnpTrinucError::RngError(error)
    }
}

impl From<TransitionMatrixError> for SnpTrinucError {
    fn from(error: TransitionMatrixError) -> Self {
        SnpTrinucError::TransMatrixError(error)
    }
}

#[derive(Debug, Copy, Clone, Eq, Hash, PartialEq, Serialize, Deserialize)]
enum SnpTrinucFrame {
    // These SNP frames represent the 16 different trinculeotide contexts, with N representing
    // any of the 4 other bases. So each SnpFrame represents 4 trinucleotide combinations, for a
    // total of 64 possible trinucleotide combinations.
    ANA, ANC, ANG, ANT,
    CNA, CNC, CNG, CNT,
    GNA, GNC, GNG, GNT,
    TNA, TNC, TNG, TNT,
}

impl SnpTrinucFrame {
    // Set the order to iterate through these
    pub fn iterator() -> Iter<'static, SnpTrinucFrame> {
        static SNP_FRAMES: [SnpTrinucFrame; 16] = [
            SnpTrinucFrame::ANA, 
            SnpTrinucFrame::ANC, 
            SnpTrinucFrame::ANG, 
            SnpTrinucFrame::ANT, 
            SnpTrinucFrame::CNA, 
            SnpTrinucFrame::CNC, 
            SnpTrinucFrame::CNG, 
            SnpTrinucFrame::CNT, 
            SnpTrinucFrame::GNA, 
            SnpTrinucFrame::GNC,
            SnpTrinucFrame::GNG, 
            SnpTrinucFrame::GNT, 
            SnpTrinucFrame::TNA, 
            SnpTrinucFrame::TNC, 
            SnpTrinucFrame::TNG, 
            SnpTrinucFrame::TNT
        ];
        SNP_FRAMES.iter()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnpTrinucModel {
    // Relative weights given to each SNP frame. Ultimately this will be imputed from data.
    snp_distr: DiscreteDistribution,
    // The transition matrix is the chance of mutating the middle base from A, C, T, or G to a
    // different base (4x4 matrix with 0s on the diagonal).
    #[serde(with = "vectorize")]
    trinuc_distros: HashMap<SnpTrinucFrame, TransitionMatrix>,
}

impl SnpTrinucModel {
    pub fn default() -> Result<Self, SnpTrinucError> {
        // Creating the default trinuc bias model for snps. In this model, all trinucleotides
        // mutate with equal probability and middle base mutates with the same probability no matter
        // the context (the default transition matrix).

        // One weight for each snp frame.
        let mut snp_weights: Vec<f64> = Vec::with_capacity(16);
        let mut trinuc_distros: HashMap<SnpTrinucFrame, TransitionMatrix> = HashMap::new();
        for frame in SnpTrinucFrame::iterator() {
            snp_weights.push(1.0);
            trinuc_distros.insert(*frame, TransitionMatrix::default()?);
        }
        let snp_distr = DiscreteDistribution::new_index_only(&snp_weights)?;
        Ok(SnpTrinucModel {
            snp_distr,
            trinuc_distros,
        })
    }

    #[allow(unused)]
    /// need to think about the gc bias model of Neat2
    fn generate_gc_bias(&self, _input_sequence: &Vec<u8>) -> Vec<usize> {
        todo!()
        // We need some way to use this model to GC bias positions of SNPs, 
        // but it's not clear yet how.
    }

    #[allow(unused)]
    /// we will write utilities to use this, eventually
    fn write_out_quality_model(&self, filename: &str) -> std::io::Result<()> {
        // Uses the serde_json crate to write out the json form of the model. This will help us
        // create base datasets from old neat data, and give us a way to write out models that are
        // generated from user data.
        Ok(model_writer(self, filename).unwrap())
    }

    #[allow(unused)]
    /// we will write utilities to use this, eventually
    fn read_in_quality_model(&self, filename: &str) -> Result<Self, SnpTrinucError> {
        // Uses the serde_json crate to read a quality model from file
        let data: SnpTrinucModel = model_gzp_reader(filename).unwrap();
        Ok(data)
    }

    pub fn generate_snp(
        &self,
        rand: f64,
        trinuc_reference: &[u8; 3],
    ) -> Result<u8, SnpTrinucError> {
        // We shouldn't have N's here. Basically, this matches the correct trinuc from the enum,
        // then uses that as the index for the trinuc matrix of interest.
        // eliminating the need to copy the rng or a pointer to it, if possible
        let matrix = self.trinuc_distros.get(&match trinuc_reference[0] {
            0 => {
                match trinuc_reference[2] {
                    0 => SnpTrinucFrame::ANA,
                    1 => SnpTrinucFrame::ANC,
                    2 => SnpTrinucFrame::ANG,
                    3 => SnpTrinucFrame::ANT,
                    _ => { 
                        error!("Trying to use trinucleotide bias on an unknown base (N).");
                        return Err(SnpTrinucError::GenerateSnpError)
                    }
                }
            },
            1 => {
                match trinuc_reference[2] {
                    0 => SnpTrinucFrame::CNA,
                    1 => SnpTrinucFrame::CNC,
                    2 => SnpTrinucFrame::CNG,
                    3 => SnpTrinucFrame::CNT,
                    _ => { 
                        error!("Trying to use trinucleotide bias on an unknown base (N).");
                        return Err(SnpTrinucError::GenerateSnpError)
                    }
                }
            },
            2 => {
                match trinuc_reference[2] {
                    0 => SnpTrinucFrame::GNA,
                    1 => SnpTrinucFrame::GNC,
                    2 => SnpTrinucFrame::GNG,
                    3 => SnpTrinucFrame::GNT,
                    _ => { 
                        error!("Trying to use trinucleotide bias on an unknown base (N).");
                        return Err(SnpTrinucError::GenerateSnpError)
                    }
                }
            },
            3 => {
                match trinuc_reference[2] {
                    0 => SnpTrinucFrame::TNA,
                    1 => SnpTrinucFrame::TNC,
                    2 => SnpTrinucFrame::TNG,
                    3 => SnpTrinucFrame::TNT,
                    _ => { 
                        error!("Trying to use trinucleotide bias on an unknown base (N).");
                        return Err(SnpTrinucError::GenerateSnpError)
                    }
                }
            },
            _ => { 
                error!("Trying to use trinucleotide bias on an unknown base (N).");
                return Err(SnpTrinucError::GenerateSnpError)
            },
        }).unwrap();

        match trinuc_reference[1] {
            0 => Ok(matrix.a_dist.sample_index(rand)? as u8),
            1 => Ok(matrix.c_dist.sample_index(rand)? as u8),
            2 => Ok(matrix.g_dist.sample_index(rand)? as u8),
            3 => Ok(matrix.t_dist.sample_index(rand)? as u8),
            _ => { 
                error!("Trying to use trinucleotide bias on an unknown base (N).");
                return Err(SnpTrinucError::GenerateSnpError)
            }
        }
    }
}


#[cfg(test)]
mod tests {

    #[test]
    fn test_snp_model() {
        println!("TODO");
        assert_eq!(1, 1)
    }
}