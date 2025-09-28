//! The following section are the models for each type of variant. In order to create the variant,
//! we need to model its statistical property. NEAT included two types of variants: SNPs and Indels.
//! Note that indels are actually two types, insertions and deletions, but they are usually classed
//! together in the literature and are usually short. Insertions are often slips that cause sections
//! to be duplicated, but NEAT made no attempt to distinguish between the types of insertions or
//! deletions, since they are treated similarly by variant calling software.
use flate2::read::GzDecoder;
use itertools::Itertools;
use log::{error, debug};
use vectorize;
use serde;
use std::hash::Hash;
use std::{
    fmt::{self, Debug}, 
    io, 
    path::PathBuf
};
use std::ops::Index;
use serde::{Deserialize, Serialize};
use thiserror::Error;
use std::collections::HashMap;
use simple_rng::NeatRngError;
use lazy_static::lazy_static;
use crate::models::lib::{model_reader, model_writer};
use crate::structs::transition_matrix::{TransitionMatrix, TransitionMatrixError};
use crate::structs::distributions::{DiscreteDistribution, DistributionErrors};
use crate::structs::nucleotides::Nucleotide::{A, C, G, T, N};
use crate::structs::nucleotides::Nucleotide;


#[derive(Error, Debug)]
pub enum SnpTrinucError {
    #[error("SNP variant model reported an RNG error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("SNP variant model reported an error from Distributions: {0}")]
    DistributionError(#[from] DistributionErrors),
    #[error("SNP variant model reported an error with transition matrices: {0}")]
    TransMatrixError(#[from] TransitionMatrixError),
    #[error("SNP variant model reported an error generating a SNP.")]
    GenerateSnpError,
    #[error("SNP Trinuc model return an IO error: {0}")]
    IoError(#[from] io::Error),
    #[error("Error Indexing the trinucleotide: {0}")]
    IndexError(usize),
    #[error("Serde error building default model: {0}")]
    SerdeError(#[from] serde_json::Error),
}

/// These are all the trinucleotide combinations.
const SNP_TRINUCS: [(Nucleotide, Nucleotide, Nucleotide); 64] = [
    (A, A, A), // AAA
    (A, A, C), // AAC
    (A, A, G), // AAG
    (A, A, T), // AAT
    (A, C, A), // ACA
    (A, C, C), // ACC
    (A, C, G), // ACG
    (A, C, T), // ACT
    (A, G, A), // AGA
    (A, G, C), // AGC
    (A, G, G), // AGG
    (A, G, T), // AGT
    (A, T, A), // ATA
    (A, T, C), // ATC
    (A, T, G), // ATG
    (A, T, T), // ATT
    (C, A, A), // CAA
    (C, A, C), // CAC
    (C, A, G), // CAG
    (C, A, T), // CAT
    (C, C, A), // CCA
    (C, C, C), // CCC
    (C, C, G), // CCG
    (C, C, T), // CCT
    (C, G, A), // CGA
    (C, G, C), // CGC
    (C, G, G), // CGG
    (C, G, T), // CGT
    (C, T, A), // CTA
    (C, T, C), // CTC
    (C, T, G), // CTG
    (C, T, T), // CTT
    (G, A, A), // GAA
    (G, A, C), // GAC
    (G, A, G), // GAG
    (G, A, T), // GAT
    (G, C, A), // GCA
    (G, C, C), // GCC
    (G, C, G), // GCG
    (G, C, T), // GCT
    (G, G, A), // GGA
    (G, G, C), // GGC
    (G, G, G), // GGG
    (G, G, T), // GGT
    (G, T, A), // GTA
    (G, T, C), // GTC
    (G, T, G), // GTG
    (G, T, T), // GTT
    (T, A, A), // TAA
    (T, A, C), // TAC
    (T, A, G), // TAG
    (T, A, T), // TAT
    (T, C, A), // TCA
    (T, C, C), // TCC
    (T, C, G), // TCG
    (T, C, T), // TCT
    (T, G, A), // TGA
    (T, G, C), // TGC
    (T, G, G), // TGG
    (T, G, T), // TGT
    (T, T, A), // TTA
    (T, T, C), // TTC
    (T, T, G), // TTG
    (T, T, T), // TTT
];

lazy_static! {
    pub static ref ALL_FRAMES: Vec<TrinucFrame> = {
        let mut all_frames = Vec::new();
        for trinuc in SNP_TRINUCS {
            all_frames.push(TrinucFrame::from(trinuc))
        }
        all_frames
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
            alias_map.insert(
                TrinucFrame::from(frame), 
                TrinucFrame::from((frame[0], N, frame[2]))
            );
        }
        alias_map
    };

    static ref ALL_CONTEXTS: Vec<TrinucFrame> = {
        let alias_map = ALIAS_MAP.clone();
        let frames: Vec<TrinucFrame> = alias_map
            .values()
            .map(|s| s.convert())
            .collect();
        let frames = frames.iter().unique().cloned().collect::<Vec<_>>();
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
                    let frame_0: usize = frame[0].into();
                    let frame_2: usize = frame[2].into();
                    if frame_0 == i && frame_2 == k {
                        frame_list.push(TrinucFrame::from(frame.clone()))
                    }
                }
                frame_map.insert(
                    TrinucFrame::from((i.into(), N, k.into())), 
                    frame_list
                );
            }
        }
        frame_map
    };

}

// We need all these derivations to write these out to file
#[derive(Clone, Debug, Eq, PartialEq, Hash, Copy, Serialize, Deserialize)]
pub enum TrinucFrame {
    Frame(Nucleotide, Nucleotide, Nucleotide),
}

impl fmt::Display for TrinucFrame {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut display: String = String::new();
        for i in 0..3 {
            display.push(self[i].into())
        }
        write!(f, "{}", display)
    }
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

impl From<[&Nucleotide; 3]> for TrinucFrame {
    fn from(trinuc: [&Nucleotide; 3]) -> Self {
        Self::Frame(*trinuc[0], *trinuc[1], *trinuc[2])
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
        let TrinucFrame::Frame(b0, b1, b2) = self;
        match i {
            0 => b0,
            1 => b1,
            2 => b2,
            _ => panic!("Index out of range: {}", i)
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnpTrinucModel {
    // Relative weights given to each SNP frame. Ultimately this will be imputed from data.
    snp_trinuc_distro: DiscreteDistribution<TrinucFrame>,
    // The transition matrix is the chance of mutating the middle base from A, C, T, or G to a
    // different base (4x4 matrix with 0s on the diagonal).
    // Each Trinuc index is a context index, with a middle "N".
    // Each row is a from base, and each corresponding column is the to base, with the intersecting
    //    f64 being the probability of transitioning from...to.
    #[serde(with = "vectorize")]
    trinuc_context_distros: HashMap<TrinucFrame, TransitionMatrix>,
}

static DATA_FILE: &'static [u8] = include_bytes!("model_data/default_trinuc_model.json.gz");

impl SnpTrinucModel {
    pub fn from_raw_data(
        snp_trinuc_weights_raw: HashMap<TrinucFrame, f64>,
        trinuc_transition_frequency: HashMap<(TrinucFrame, TrinucFrame), f64>,
    ) -> Result<Self, SnpTrinucError> {
        // Inputs:
        // snp_trinuc_weights: This is the probability that the trinucleotide combination mutates
        //    at all. A 0.0 value here indicates the trinucleotide was not observed in the data.
        // trinuc_transition_frequency: This is the probability that Trinucleotide combination A
        //    transitions into trinucleotide combination B. There may be 0s here where we did not
        //    observe one transitioning to the other.
        // To account for potential zeroes in the data, we will add in a 0.001 in place of each
        //    zero. This should get ensmallated when normalization happens, making them even less
        //    likely. We can adjust this epsilon in response to user input.
        let epsilon = 0.001;
        let all_contexts = ALL_CONTEXTS.clone();
        let context_frame_map = CONTEXT_FRAME_MAP.clone();
        let mut trinuc_context_distros: HashMap<TrinucFrame, TransitionMatrix> = HashMap::new();
        for context in &all_contexts {
            let mut temp_trans_matrix: [[f64; 4] ;4] = [
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
            ];
            let frames = &context_frame_map[context];
            for from_frame in frames {
                for nuc2 in 0..4 {
                    let to_frame = TrinucFrame::from(
                        &[context[0], Nucleotide::from(nuc2), context[2]]
                    );
                    temp_trans_matrix[from_frame[1] as usize][nuc2] = trinuc_transition_frequency[
                        &(*from_frame, to_frame)
                    ];
                }
            }
            let trans_matrix = TransitionMatrix::from(
                temp_trans_matrix[Nucleotide::A as usize],
                temp_trans_matrix[Nucleotide::C as usize],
                temp_trans_matrix[Nucleotide::G as usize],
                temp_trans_matrix[Nucleotide::T as usize],
            )?;

            trinuc_context_distros.insert(context.clone(), trans_matrix);
        }
        let (snp_trinuc_weights, snp_trinuc_values) = {
            let mut temp_weights: Vec<f64> = Vec::new();
            let all_frames = ALL_FRAMES.clone();
            for frame in &all_frames {
                if snp_trinuc_weights_raw.contains_key(&frame) {
                    temp_weights.push(snp_trinuc_weights_raw[&frame] + epsilon);
                } else {
                    temp_weights.push(epsilon);
                }
            }
            (temp_weights, all_frames)
        };
        let snp_trinuc_distro = DiscreteDistribution::new(
            &snp_trinuc_weights,
            &snp_trinuc_values,
        )?;
        Ok(Self {
            snp_trinuc_distro,
            trinuc_context_distros,
        })
    }

    pub fn default_minimal() -> Result<Self, SnpTrinucError> {
        // Creating the default trinuc bias model for snps. In this model, all trinucleotides
        // mutate with equal probability and middle base mutates with the same probability no matter
        // the context (the default transition matrix).

        // One weight for each snp frame, making them all equally
        let mut snp_weights: Vec<f64> = Vec::with_capacity(16);
        let mut trinuc_distros: HashMap<TrinucFrame, TransitionMatrix> = HashMap::new();
        let all_contexts = ALL_CONTEXTS.clone();
        debug!("All contexts: {:?}", &all_contexts);
        let all_frames = ALL_FRAMES.clone();
        for frame in &all_contexts {
            trinuc_distros.insert(frame.clone() ,TransitionMatrix::default()?);
            snp_weights.push(1.0);
        }
        let snp_distr = DiscreteDistribution::new(
            &snp_weights,
            // Instead of enumerating all those frames, we will just use the index here
            &all_frames,
        )?;
        Ok(SnpTrinucModel {
            snp_trinuc_distro: snp_distr,
            trinuc_context_distros: trinuc_distros,
        })
    }

    pub fn default() -> Result<Self, SnpTrinucError> {
        // This default will read in data from the original NEAT trinuc model to create a more realistic trinuc
        // bias.
        let reader = GzDecoder::new(DATA_FILE);
        let data: SnpTrinucModel = serde_json::from_reader(reader)
            .map_err(SnpTrinucError::SerdeError)?;
        Ok(data)
    }

    #[allow(unused)]
    /// we will write utilities to use this, eventually
    fn from_file(&self, filename: &PathBuf) -> Result<Self, SnpTrinucError> {
        // Uses the serde_json crate to read a quality model from file
        let data: SnpTrinucModel = model_reader(filename).unwrap();
        Ok(data)
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
    fn write_out(&self, filename: &PathBuf) -> std::io::Result<()> {
        // Uses the serde_json crate to write out the json form of the model. This will help us
        // create base datasets from old neat data, and give us a way to write out models that are
        // generated from user data.
        Ok(model_writer(self, filename).unwrap())
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
        let alias_map = &ALIAS_MAP.clone();
        let alias_trinuc = alias_map[&trinuc].clone();
        let matrix = self.trinuc_context_distros[&alias_trinuc].clone();

        let nuc: Nucleotide = matrix[&trinuc[1]].sample(rand)?.into();
        Ok(nuc)
    }
}


#[cfg(test)]
mod tests {
    use std::fs;
    use super::*;

    #[test]
    fn test_model_write_read() {
        let output_file: PathBuf = PathBuf::from("test.json.gz");
        let model = SnpTrinucModel::default().unwrap();
        let frame = TrinucFrame::from((A, N, G));
        assert!(model.trinuc_context_distros.contains_key(&frame));
        let result = model.write_out(&output_file);
        assert_eq!(result.unwrap(), ());
        fs::remove_file(output_file).unwrap();
    }
}