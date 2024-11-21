use std::collections::HashMap;
use common::structs::variants::{Variant, VariantType};
use simple_rng::{NeatRng, DiscreteDistribution};

pub fn apply_variants(
    raw_sequence: &[u8],
    relevant_variants: HashMap<usize, &Variant>,
    start: usize,
    mut rng: NeatRng,
) -> Vec<u8> {
    // since this is a fragment, our 0 index is relative to our start position. So an index of 3
    // is at position 3 + start relative to the entire reference.
    // todo modify so that only half the reads of a heterzygote get the mutation.
    let mut positions_to_mutate = Vec::new();
    for location in relevant_variants.keys() {
        positions_to_mutate.push((*location).clone() - start)
    }
    let mut mutated_sequence = Vec::new();
    for mut position in 0..raw_sequence.len() {
        if positions_to_mutate.contains(&position) {
            let variant = relevant_variants.get(&(position + start)).unwrap();
            let mut probability = 0.5;
            if variant.is_homozygous() {
                probability = 1.0;
            }
            // if it's heterozygous, we only want the variant half the time.
            if (&mut rng).gen_bool(probability) {
                match variant.variant_type {
                    VariantType::Insertion => {
                        for base in &variant.alternate {
                            mutated_sequence.push(base.clone());
                        }
                    },
                    VariantType::Deletion => {
                        // we keep the initial base
                        mutated_sequence.push(raw_sequence[position].clone());
                        // since we always add one below for position, we only need to skip
                        // length - 1 bases to capture the full deletion
                        position += variant.reference.len() - 1;
                    },
                    VariantType::SNP => {
                        // We store the snp as a vector of length 1, so the 0 index is always correct.
                        mutated_sequence.push(variant.alternate[0].clone())
                    },
                }
            } else {
                // half the reads will have the ref still, since it's heterozygous.
                mutated_sequence.push(raw_sequence[position].clone())
            }
        } else {
            mutated_sequence.push(raw_sequence[position].clone())
        }
        position += 1
    }
    mutated_sequence
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_apply_variants() {
        let raw_sequence = [0, 1, 2, 3, 3, 0, 3, 2];
        let relevant_variants = HashMap::from([
            (1, &Variant::new(
                VariantType::SNP,
                1,
                &vec![3],
                &vec![1],
                &vec![1, 0],
            )),
            (3, &Variant::new(
                VariantType::Deletion,
                3,
                &vec![4, 4],
                &vec![4],
                &vec![0,1],
            ))
        ]);
    }
}