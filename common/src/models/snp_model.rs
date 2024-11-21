use std::collections::HashMap;
use std::slice::Iter;
use simple_rng::{DiscreteDistribution, NeatRng};
use crate::structs::transition_matrix::TransitionMatrix;
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
    ANA, ANC, ANG, ANT,
    CNA, CNC, CNG, CNT,
    GNA, GNC, GNG, GNT,
    TNA, TNC, TNG, TNT,
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
    // Relative weights given to each SNP frame. Ultimately this will be imputed from data.
    snp_distr: DiscreteDistribution,
    // The transition matrix is the chance of mutating the middle base from A, C, T, or G to a
    // different base (4x4 matrix with 0s on the diagonal).
    trinuc_distros: HashMap<SnpFrame, TransitionMatrix>,
}

impl SnpModel {
    pub fn new() -> Self {
        // Creating the default trinuc bias model for snps. In this model, all trinucleotides
        // mutate with equal probability and middle base mutates with the same probability no matter
        // the context (the default transition matrix).

        // One weight for each snp frame.
        let mut snp_weights: Vec<f64> = Vec::with_capacity(16);
        let mut trinuc_distros: HashMap<SnpFrame, TransitionMatrix> = HashMap::new();
        for frame in SnpFrame::iterator() {
            snp_weights.push(1.0);
            trinuc_distros.insert(*frame, TransitionMatrix::default());
        }
        let snp_distr = DiscreteDistribution::new(&snp_weights);
        SnpModel {
            snp_distr,
            trinuc_distros,
        }
    }
    fn generate_bias(&self, _input_sequence: &Vec<u8>) -> Vec<usize> {
        todo!()
        // We need some way to use this model to bias positions of SNPs, but it's not clear yet how.
    }

    pub fn generate_snp(
        &self,
        trinuc_reference: &[u8; 3],
        rng: &mut NeatRng
    ) -> u8 {
        // We shouldn't have N's here. Basically, this matches the correct trinuc from the enum,
        // then uses that as the index for the trinuc matrix of interest.
        let matrix = self.trinuc_distros.get(&match trinuc_reference[0] {
            0 => {
                match trinuc_reference[2] {
                    0 => ANA,
                    1 => ANC,
                    2 => ANG,
                    3 => ANT,
                    _ => { panic!("Trying to use trinucleotide bias on an unknown base (N).") }
                }
            },
            1 => {
                match trinuc_reference[2] {
                    0 => CNA,
                    1 => CNC,
                    2 => CNG,
                    3 => CNT,
                    _ => { panic!("Trying to use trinucleotide bias on an unknown base (N).") }
                }
            },
            2 => {
                match trinuc_reference[2] {
                    0 => GNA,
                    1 => GNC,
                    2 => GNG,
                    3 => GNT,
                    _ => { panic!("Trying to use trinucleotide bias on an unknown base (N).") }
                }
            },
            3 => {
                match trinuc_reference[2] {
                    0 => TNA,
                    1 => TNC,
                    2 => TNG,
                    3 => TNT,
                    _ => { panic!("Trying to use trinucleotide bias on an unknown base (N).") }
                }
            },
            _ => { panic!("Trying to use trinucleotide bias on an unknown base (N).") },
        }).unwrap();

        match trinuc_reference[1] {
            0 => matrix.a_dist.sample(rng) as u8,
            1 => matrix.c_dist.sample(rng) as u8,
            2 => matrix.g_dist.sample(rng) as u8,
            3 => matrix.t_dist.sample(rng) as u8,
            _ => { panic!("Trying to use trinucleotide bias on an unknown base (N).") },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_snp_model() {
        println!("TODO");
        assert_eq!(1, 1)
    }
}