use std::collections::HashMap;
use std::slice::Iter;
use rand::distributions::{Distribution, WeightedIndex};
use structs::nucleotides::Nuc;
use structs::transition_matrix::TransitionMatrix;
use rand_chacha::ChaCha20Rng;
use self::SnpFrame::*;
// The following section are the models for each type of variant. In order to create the variant,
// we need to model its statistical property. NEAT included two types of variants: SNPs and Indels.
// Note that indels are actually two types, insertions and deletions, but they are usually classed
// together in the literature and are usually short. Insertions are often slips that cause sections
// to be duplicated, but NEAT made no attempt to distinguish between the types of insertions or
// deletions, since they are treated similarly by variant calling software.
#[derive(Debug, Copy, Clone, Eq, Hash, PartialEq)]
enum SnpFrame {
    // These SNP frames represent the 16 different trinculeotide contexts, with N representing
    // any of the 4 other bases. So each SnpFrame represents 4 trinucleotide combinations, for a
    // total of 64 possible trinucleotide combinations.
    ANA,
    ANC,
    ANG,
    ANT,
    CNA,
    CNC,
    CNG,
    CNT,
    GNA,
    GNC,
    GNG,
    GNT,
    TNA,
    TNC,
    TNG,
    TNT,
}

impl SnpFrame {
    pub fn iterator() -> Iter<'static, SnpFrame> {
        static SNP_FRAMES: [SnpFrame; 16] = [
            ANA, ANC, ANG, ANT, CNA, CNC, CNG, CNT, GNA, GNC, GNG, GNT, TNA, TNC, TNG, TNT
        ];
        SNP_FRAMES.iter()
    }
}

#[derive(Debug, Clone)]
pub struct SnpModel {
    // Relative weights given to each SNP frame.
    snp_weights: HashMap<SnpFrame, usize>,
    // The transition matrix is the chance of mutating the middle base from A, C, T, or G to a
    // different base (4x4 matrix with 0s on the diagonal).
    trinuc_matrix: HashMap<SnpFrame, TransitionMatrix>,
}

impl SnpModel {
    pub fn new() -> Self {
        // Creating the default trinuc bias model for snps. In this model, all trinucleotides
        // mutate with equal probability and middle base mutates with the same probability no matter
        // the context (the default transition matrix).
        let mut snp_weights: HashMap<SnpFrame, usize> = HashMap::new();
        let mut trinuc_matrix: HashMap<SnpFrame, TransitionMatrix> = HashMap::new();
        for frame in SnpFrame::iterator() {
            snp_weights.insert(*frame, 1);
            trinuc_matrix.insert(*frame, TransitionMatrix::new());
        }
        SnpModel {
            snp_weights,
            trinuc_matrix,
        }
    }
    fn generate_bias(&self, input_sequence: &Vec<Nuc>) -> Vec<usize> {
        todo!()
        // We need some way to use this model to bias positions of SNPs, but it's not clear yet how.
    }

    pub fn generate_snp(
        &self,
        trinuc_reference: &[Nuc],
        rng: &mut ChaCha20Rng
    ) -> Nuc {
        // We shouldn't have N's here. Basically, this matches the correct trinuc from the enum,
        // then uses that as the index for the trinuc matrix of interest.
        let matrix = self.trinuc_matrix.get(&match trinuc_reference[0] {
            Nuc::A => {
                match trinuc_reference[2] {
                    Nuc::A => ANA,
                    Nuc::C => ANC,
                    Nuc::G => ANG,
                    Nuc::T => ANT,
                    _ => { panic!("Trying to use trinucleotide bias on an unknown base (N).") }
                }
            },
            Nuc::C => {
                match trinuc_reference[2] {
                    Nuc::A => CNA,
                    Nuc::C => CNC,
                    Nuc::G => CNG,
                    Nuc::T => CNT,
                    _ => { panic!("Trying to use trinucleotide bias on an unknown base (N).") }
                }
            },
            Nuc::G => {
                match trinuc_reference[2] {
                    Nuc::A => GNA,
                    Nuc::C => GNC,
                    Nuc::G => GNG,
                    Nuc::T => GNT,
                    _ => { panic!("Trying to use trinucleotide bias on an unknown base (N).") }
                }
            },
            Nuc::T => {
                match trinuc_reference[2] {
                    Nuc::A => TNA,
                    Nuc::C => TNC,
                    Nuc::G => TNG,
                    Nuc::T => TNT,
                    _ => { panic!("Trying to use trinucleotide bias on an unknown base (N).") }
                }
            },
            _ => { panic!("Trying to use trinucleotide bias on an unknown base (N).") },
        }).unwrap();

        let weights = match trinuc_reference[1] {
            Nuc::A => matrix.a_weights.clone(),
            Nuc::C => matrix.c_weights.clone(),
            Nuc::G => matrix.g_weights.clone(),
            Nuc::T => matrix.t_weights.clone(),
            _ => { panic!("Trying to use trinucleotide bias on an unknown base (N).") },
        };

        let dist = WeightedIndex::new(weights).unwrap();
        match dist.sample(rng) {
            0 => Nuc::A,
            1 => Nuc::C,
            2 => Nuc::G,
            3 => Nuc::T,
            _ => { panic!("Invalid nucleotide reference.") },
        }
    }
}