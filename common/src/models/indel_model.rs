use rand::distributions::{WeightedIndex, Distribution};
use rand::Rng;
// The following section are the models for each type of variant. In order to create the variant,
// we need to model its statistical property. NEAT included two types of variants: SNPs and Indels.
// Note that indels are actually two types, insertions and deletions, but they are usually classed
// together in the literature and are usually short. Insertions are often slips that cause sections
// to be duplicated, but NEAT made no attempt to distinguish between the types of insertions or
// deletions, since they are treated similarly by variant calling software.
use structs::nucleotides::Nuc;
use rand_chacha::ChaCha20Rng;

#[derive(Debug, Clone)]
pub struct IndelModel {
    // Based what was in the original NEAT
    insertion_probability: f64,
    ins_lengths: Vec<usize>,
    ins_weights: Vec<usize>,
    del_lengths: Vec<usize>,
    del_weights: Vec<usize>,
}

impl IndelModel {
    pub fn new() -> Self {
        // Default Indel model from the original NEAT
        let insertion_probability = 0.6;
        let ins_lengths: Vec<usize> = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let ins_weights: Vec<usize> = vec![10, 10, 20, 5, 5, 5, 5, 5, 5, 5];
        let del_lengths: Vec<usize> = vec![1, 2, 3, 4, 5];
        let del_weights: Vec<usize> = vec![3, 2, 2, 2, 1];

        IndelModel {
            insertion_probability,
            ins_lengths,
            ins_weights,
            del_lengths,
            del_weights,
        }
    }

    pub fn generate_new_indel_length(&self, mut rng: &mut ChaCha20Rng) -> i64 {
        let rand_num = rng.gen::<f64>();
        println!("number selected was {}", &rand_num);
        let is_insertion = { rand_num < self.insertion_probability };
        if is_insertion {
            let mut dist = WeightedIndex::new(&self.ins_weights).unwrap();
            self.ins_lengths[dist.sample(&mut rng)] as i64
        } else {
            let mut dist = WeightedIndex::new(&self.del_weights).unwrap();
            -(self.del_lengths[dist.sample(&mut rng)] as i64)
        }
    }

    pub fn generate_random_insertion(&self, length: usize, rng: &mut ChaCha20Rng) -> Vec<Nuc> {
        // We could refine this with a nucleotide bias matrix. Maybe it would make a difference,
        // but probably not, since the presence of the insertion is more important than it's content,
        // for this use. If there was a call for it, maybe.
        let mut insertion_vec = Vec::new();
        for _ in 0..length {
            insertion_vec.push(
                // gen_range is equally weighted, which works fine for now.
                match rng.gen_range(0..=3) {
                0 => Nuc::A,
                1 => Nuc::C,
                2 => Nuc::G,
                3 => Nuc::T,
                // Since our range is 0..=3, this is unreachable, but Rust needs it.
                _ => panic!("gen_range generated an invalid number."),
                }
            )
        }
        insertion_vec
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use structs::nucleotides::Nuc::*;
    use structs::variants::{VariantType, Variant};
    use structs::transition_matrix::TransitionMatrix;
    use rand_core::SeedableRng;

    #[test]
    fn test_generate_variant() {
        let test_indel_model = IndelModel::new();
        let input_sequence = vec![A, C, G, G, A];
        let transition_matrix = TransitionMatrix::new();
        let mut rng = ChaCha20Rng::seed_from_u64(0);
        let test_variant: Variant = test_indel_model.generate_variant(
            "chr1",
            &input_sequence,
            0,
            &transition_matrix,
            &mut rng,
        );
        assert_eq!(test_variant.variant_type, VariantType::Indel);
        assert_eq!(test_variant.reference, input_sequence);
    }
}
