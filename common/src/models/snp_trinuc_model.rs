//! The following section are the models for each type of variant. In order to create the variant,
//! we need to model its statistical property. NEAT included two types of variants: SNPs and Indels.
//! Note that indels are actually two types, insertions and deletions, but they are usually classed
//! together in the literature and are usually short. Insertions are often slips that cause sections
//! to be duplicated, but NEAT made no attempt to distinguish between the types of insertions or
//! deletions, since they are treated similarly by variant calling software.
use log::error;
use vectorize;
use serde;
use std::hash::Hash;
use std::io;
use std::ops::Index;
use serde::{Deserialize, Serialize};
use thiserror::Error;
use std::collections::HashMap;
use simple_rng::NeatRngError;
use lazy_static::lazy_static;
use crate::models::lib::{model_gzp_reader, model_writer};
use crate::structs::transition_matrix::{TransitionMatrix, TransitionMatrixError};
use crate::structs::distributions::{DiscreteDistribution, DistributionErrors};
use crate::structs::nucleotides::Nucleotide;

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
    #[error("SNP Trinuc model return an IO error: {0}")]
    IoError(io::Error),
}

impl From<io::Error> for SnpTrinucError {
    fn from(error: io::Error) -> Self {
        SnpTrinucError::IoError(error)
    }
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

/// These are all the trinucleotide combinations.
const SNP_TRINUCS: [(Nucleotide, Nucleotide, Nucleotide); 64] = [
    (Nucleotide::A, Nucleotide::A, Nucleotide::A), // AAA
    (Nucleotide::A, Nucleotide::A, Nucleotide::C), // AAC
    (Nucleotide::A, Nucleotide::A, Nucleotide::G), // AAG
    (Nucleotide::A, Nucleotide::A, Nucleotide::T), // AAT
    (Nucleotide::A, Nucleotide::C, Nucleotide::A), // ACA
    (Nucleotide::A, Nucleotide::C, Nucleotide::C), // ACC
    (Nucleotide::A, Nucleotide::C, Nucleotide::G), // ACG
    (Nucleotide::A, Nucleotide::C, Nucleotide::T), // ACT
    (Nucleotide::A, Nucleotide::G, Nucleotide::A), // AGA
    (Nucleotide::A, Nucleotide::G, Nucleotide::C), // AGC
    (Nucleotide::A, Nucleotide::G, Nucleotide::G), // AGG
    (Nucleotide::A, Nucleotide::G, Nucleotide::T), // AGT
    (Nucleotide::A, Nucleotide::T, Nucleotide::A), // ATA
    (Nucleotide::A, Nucleotide::T, Nucleotide::C), // ATC
    (Nucleotide::A, Nucleotide::T, Nucleotide::G), // ATG
    (Nucleotide::A, Nucleotide::T, Nucleotide::T), // ATT
    (Nucleotide::C, Nucleotide::A, Nucleotide::A), // CAA
    (Nucleotide::C, Nucleotide::A, Nucleotide::C), // CAC
    (Nucleotide::C, Nucleotide::A, Nucleotide::G), // CAG
    (Nucleotide::C, Nucleotide::A, Nucleotide::T), // CAT
    (Nucleotide::C, Nucleotide::C, Nucleotide::A), // CCA
    (Nucleotide::C, Nucleotide::C, Nucleotide::C), // CCC
    (Nucleotide::C, Nucleotide::C, Nucleotide::G), // CCG
    (Nucleotide::C, Nucleotide::C, Nucleotide::T), // CCT
    (Nucleotide::C, Nucleotide::G, Nucleotide::A), // CGA
    (Nucleotide::C, Nucleotide::G, Nucleotide::C), // CGC
    (Nucleotide::C, Nucleotide::G, Nucleotide::G), // CGG
    (Nucleotide::C, Nucleotide::G, Nucleotide::T), // CGT
    (Nucleotide::C, Nucleotide::T, Nucleotide::A), // CTA
    (Nucleotide::C, Nucleotide::T, Nucleotide::C), // CTC
    (Nucleotide::C, Nucleotide::T, Nucleotide::G), // CTG
    (Nucleotide::C, Nucleotide::T, Nucleotide::T), // CTT
    (Nucleotide::G, Nucleotide::A, Nucleotide::A), // GAA
    (Nucleotide::G, Nucleotide::A, Nucleotide::C), // GAC
    (Nucleotide::G, Nucleotide::A, Nucleotide::G), // GAG
    (Nucleotide::G, Nucleotide::A, Nucleotide::T), // GAT
    (Nucleotide::G, Nucleotide::C, Nucleotide::A), // GCA
    (Nucleotide::G, Nucleotide::C, Nucleotide::C), // GCC
    (Nucleotide::G, Nucleotide::C, Nucleotide::G), // GCG
    (Nucleotide::G, Nucleotide::C, Nucleotide::T), // GCT
    (Nucleotide::G, Nucleotide::G, Nucleotide::A), // GGA
    (Nucleotide::G, Nucleotide::G, Nucleotide::C), // GGC
    (Nucleotide::G, Nucleotide::G, Nucleotide::G), // GGG
    (Nucleotide::G, Nucleotide::G, Nucleotide::T), // GGT
    (Nucleotide::G, Nucleotide::T, Nucleotide::A), // GTA
    (Nucleotide::G, Nucleotide::T, Nucleotide::C), // GTC
    (Nucleotide::G, Nucleotide::T, Nucleotide::G), // GTG
    (Nucleotide::G, Nucleotide::T, Nucleotide::T), // GTT
    (Nucleotide::T, Nucleotide::A, Nucleotide::A), // TAA
    (Nucleotide::T, Nucleotide::A, Nucleotide::C), // TAC
    (Nucleotide::T, Nucleotide::A, Nucleotide::G), // TAG
    (Nucleotide::T, Nucleotide::A, Nucleotide::T), // TAT
    (Nucleotide::T, Nucleotide::C, Nucleotide::A), // TCA
    (Nucleotide::T, Nucleotide::C, Nucleotide::C), // TCC
    (Nucleotide::T, Nucleotide::C, Nucleotide::G), // TCG
    (Nucleotide::T, Nucleotide::C, Nucleotide::T), // TCT
    (Nucleotide::T, Nucleotide::G, Nucleotide::A), // TGA
    (Nucleotide::T, Nucleotide::G, Nucleotide::C), // TGC
    (Nucleotide::T, Nucleotide::G, Nucleotide::G), // TGG
    (Nucleotide::T, Nucleotide::G, Nucleotide::T), // TGT
    (Nucleotide::T, Nucleotide::T, Nucleotide::A), // TTA
    (Nucleotide::T, Nucleotide::T, Nucleotide::C), // TTC
    (Nucleotide::T, Nucleotide::T, Nucleotide::G), // TTG
    (Nucleotide::T, Nucleotide::T, Nucleotide::T), // TTT
];

lazy_static! {
    static ref ALL_FRAMES: Vec<(Nucleotide, Nucleotide, Nucleotide)> = {
        Vec::from(SNP_TRINUCS)
    };

    static ref ALIAS_MAP: HashMap<TrinucFrame, TrinucFrame> = {
        // This builds an alias map. For each frame in all frames, 
        // We assign a value of the context featureing N in the middle position.
        //      AAA => ANA ([0, 0, 0]: [0, 5, 0])
        //      ATA => ANA ([0, 3, 0]: [0, 5, 0])
        //      GTC => GNC ([2, 3, 1]: [2, 5, 1])
        // and so on. Building this to save time
        let all_frames = ALL_FRAMES.clone();
        let mut alias_map = HashMap::new();
        for frame in all_frames {
            alias_map.insert(TrinucFrame::from(frame), TrinucFrame::from((frame.0, Nucleotide::N, frame.2)));
        }
        alias_map
    };

    static ref ALL_CONTEXTS: Vec<TrinucFrame> = {
        let alias_map = ALIAS_MAP.clone();
        let frames: Vec<TrinucFrame> = alias_map.keys().map(|s| s.convert()).collect();
        frames
    };

    static ref CONTEXT_FRAME_MAP: HashMap<TrinucFrame, Vec<TrinucFrame>> = {
        // This builds the contexts for the trinucleotides. The 5 in the middle
        // represents any one of 4 trinculeotides (unknown). 
        //    ANA => [AAA, ACA, AGA, ATA]
        //    ANT => [AAT, ACT, AGT, ATT]
        // and so on;
        let mut frame_map: HashMap<TrinucFrame, Vec<TrinucFrame>> = HashMap::new();
        for i in 0..4 {
            for k in 0..4 {
                let mut frame_list = Vec::new();
                for frame in ALL_FRAMES.clone().into_iter() {
                    let frame_0: usize = frame.0.into();
                    let frame_2: usize = frame.2.into();
                    if frame_0 == i && frame_2 == k {
                        frame_list.push(TrinucFrame::from(frame.clone()))
                    }
                }
                frame_map.insert(
                    TrinucFrame::from((Nucleotide::from(i), Nucleotide::N, Nucleotide::from(k))), 
                    frame_list
                );
            }
        }
        frame_map
    };

}

// We need all these derivations to write these out to file
#[derive(Debug, Eq, PartialEq, Hash, Clone, Copy, Serialize, Deserialize)]
enum TrinucFrame {
    Frame(Nucleotide, Nucleotide, Nucleotide)
}

impl From<(Nucleotide, Nucleotide, Nucleotide)> for TrinucFrame {
    fn from(
        (n1, n2, n3): (Nucleotide, Nucleotide, Nucleotide)
    ) -> Self {
        Self::Frame(n1, n2, n3)
    }
}

impl TrinucFrame {
    pub fn convert(self) -> TrinucFrame {
        TrinucFrame::from(self)
    }
}

impl From<&[Nucleotide; 3]> for TrinucFrame {
    fn from(trinuc: &[Nucleotide; 3]) -> Self {
        Self::Frame(trinuc[0], trinuc[1], trinuc[2])
    }
}

impl Index<usize> for TrinucFrame {
    type Output = Nucleotide;
    fn index(&self, i: usize) -> &Nucleotide {
        match i {
            0 => &self[0],
            1 => &self[1],
            2 => &self[2],
            _ => panic!("Index out of range: {}", i)
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnpTrinucModel {
    // Relative weights given to each SNP frame. Ultimately this will be imputed from data.
    snp_distr: DiscreteDistribution,
    // The transition matrix is the chance of mutating the middle base from A, C, T, or G to a
    // different base (4x4 matrix with 0s on the diagonal).
    #[serde(with = "vectorize")]
    trinuc_distros: HashMap<TrinucFrame, TransitionMatrix>,
}

impl SnpTrinucModel {
    pub fn default_minimal() -> Result<Self, SnpTrinucError> {
        // Creating the default trinuc bias model for snps. In this model, all trinucleotides
        // mutate with equal probability and middle base mutates with the same probability no matter
        // the context (the default transition matrix).

        // One weight for each snp frame, making them all equally
        let mut snp_weights: Vec<f64> = Vec::with_capacity(16);
        let mut trinuc_distros: HashMap<TrinucFrame, TransitionMatrix> = HashMap::new();
        let all_contexts = ALL_CONTEXTS.clone();
        let all_frames = ALL_FRAMES.clone();
        for frame in &all_contexts {
            trinuc_distros.insert(frame.clone() ,TransitionMatrix::default()?);
            snp_weights.push(1.0);
        }
        let snp_distr = DiscreteDistribution::new(
            &snp_weights,
            // Instead of enumerating all those frames, we will just use the index here
            &(0..all_frames.len()).collect(),
        )?;
        Ok(SnpTrinucModel {
            snp_distr,
            trinuc_distros,
        })
    }

    pub fn default() -> Result<Self, SnpTrinucError> {
        // This default will read in data from the original NEAT trinuc model to create a more realistic trinuc
        // bias.
        todo!()
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
        trinuc_reference: &[Nucleotide; 3],
    ) -> Result<Nucleotide, SnpTrinucError> {
        // We shouldn't have N's here. Basically, this matches the correct trinuc from the enum,
        // then uses that as the index for the trinuc matrix of interest.
        // eliminating the need to copy the rng or a pointer to it, if possible
        let trinuc = TrinucFrame::from(trinuc_reference);
        let alias_map = &ALIAS_MAP;
        let alias_trinuc = alias_map[&trinuc].clone();
        let matrix = &self.trinuc_distros[&alias_trinuc];

        let nuc: Nucleotide = matrix[&trinuc[1]].sample(rand)?.into();
        Ok(nuc)
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