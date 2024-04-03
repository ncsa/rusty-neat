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
use rand_core::RngCore;

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
        // Default Indel model from the original NEAT, scaled by the lowest value and rounded.
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
        let is_insertion = { rand_num < self.insertion_probability };
        if is_insertion {
            let dist = WeightedIndex::new(&self.ins_weights).unwrap();
            self.ins_lengths[dist.sample(&mut rng)] as i64
        } else {
            let dist = WeightedIndex::new(&self.del_weights).unwrap();
            -(self.del_lengths[dist.sample(&mut rng)] as i64)
        }
    }
}

pub fn generate_random_insertion(length: usize, rng: &mut ChaCha20Rng) -> Vec<Nuc> {
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
                // Since our range is `0..=3`, this is unreachable, but Rust needs it.
                _ => panic!("gen_range generated an invalid number."),
            }
        )
    }
    insertion_vec
}

#[cfg(test)]
mod tests {
    use super::*;
    use structs::nucleotides::Nuc::*;
    use rand_core::SeedableRng;

    #[test]
    fn test_indel_model() {
        let indel_model = IndelModel {
            insertion_probability: 0.75,
            ins_lengths: vec![1, 2],
            ins_weights: vec![100, 1],
            del_lengths: vec![3, 4],
            del_weights: vec![1, 1000],
        };
        let mut rng = ChaCha20Rng::seed_from_u64(0);
        let length = indel_model.generate_new_indel_length(&mut rng);
        if length > 0 {
            assert!((length == 1) || (length == 2));
        } else {
            assert!((length == -3) || (length == -4));
        }
    }

    #[test]
    fn test_generate_insertion () {
        let mut rng = ChaCha20Rng::seed_from_u64(0);

        let insertion = generate_random_insertion(10, &mut rng);
        assert_eq!(insertion.len(), 10);
        assert!(!insertion.contains(&N));

        let insertion2 = generate_random_insertion(12, &mut rng);
        assert_eq!(insertion2.len(), 12);
        assert!(!insertion2.contains(&N));

        println!("{:?} and {:?}", insertion, insertion2);
    }
}
