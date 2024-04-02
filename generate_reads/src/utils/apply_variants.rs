use std::collections::HashMap;
use rand::Rng;
use rand_chacha::ChaChaRng;
use common::structs::nucleotides::Nuc;
use common::structs::variants::{Variant, VariantType};

pub fn apply_variants(
    raw_sequence: &[Nuc],
    relevant_variants: HashMap<&usize, &Variant>,
    start: usize,
    mut rng: ChaChaRng,
) -> Vec<Nuc> {
    // since this is a fragment, our 0 index is relative to our start position. So an index of 3
    // is at position 3 + start relative to the entire reference.
    let mut positions_to_mutate = Vec::new();
    for location in relevant_variants.keys().into_iter() {
        positions_to_mutate.push(location.clone() - start)
    }
    let mut mutated_sequence = Vec::new();
    for mut position in 0..raw_sequence.len() {
        if positions_to_mutate.contains(&position) {
            let variant = relevant_variants.get(&(position + start)).unwrap();
            let probability = 0.5;
            if variant.is_homozygous {
                let probability = 1.0;
            }
            // if it's heterozygous, we only want the variant half the time.
            if (&mut rng).gen_bool(probability) {
                match variant.variant_type {
                    VariantType::Indel => {
                        if variant.is_insertion() {
                            // insertion
                            for base in &variant.alternate {
                                mutated_sequence.push(base.clone());
                            }
                        } else {
                            // deletion
                            // we keep the initial base
                            mutated_sequence.push(raw_sequence[position].clone());
                            // since we always add one below for position, we only need to skip
                            // length - 1 bases to capture the full deletion
                            position += variant.reference.len() - 1;
                        }
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