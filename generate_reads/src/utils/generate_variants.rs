use std::collections::HashMap;
use itertools::Itertools;
use rand::distributions::{WeightedIndex, Distribution};
use rand_chacha::ChaCha20Rng;
use common::models::mutation_model::MutationModel;
use common::structs::variants::Variant;
use common::structs::nucleotides::Nuc;
use common::structs::fasta_map::SequenceBlock;
use common::structs::fasta_map::BlockType::*;

pub fn generate_variants(
    contig_map: (usize, Vec<SequenceBlock>),
    reference_sequence: &Vec<Nuc>,
    mutation_model: &mut MutationModel,
    ploidy: usize,
    mut rng: ChaCha20Rng,
) -> (usize, Vec<SequenceBlock>) {
    // todo work in the reference map, fewer calls to the ref sequence itself. Maybe we don't even
    //   need it for this.
    // Variants for this contig, organized by location. Each location may have only one variant.
    // And we may need to implement a further check that there aren't other types of overlaps.
    // This will simply overwrite any existing value, so we just capture the last variant at each
    // location.
    let reference_length = contig_map.0;
    // Determine a base number of mutations
    let mut number_of_mutations = (*reference_length as f64 * mutation_model.mutation_rate)
        .round() as usize;
    // Because of how N's are represented, we will not add variants to any non-n section
    let mut blocks: Vec<SequenceBlock> = contig_map.1;

    let weights: Vec<usize> = get_weights(&blocks);
    let positions: Vec<usize> = (0..reference_length).collect();

    let dist = WeightedIndex::new(&weights).unwrap();
    while number_of_mutations > 0 {
        let location: usize = positions[dist.sample(&mut rng)];
        blocks = mutation_model.generate_mutation(
            reference_sequence,
            &blocks,
            location,
            ploidy,
            rng.clone(),
        );
        // For now, we'll count duplicates just to get through it.
        // basically this means we successfully added a mutation.
        number_of_mutations -= 1;
    }
    contig_map.clone()
}

fn get_weights(fasta_map: &Vec<SequenceBlock>) -> Vec<usize> {
    // todo incorporate GC-bias and trinuc bias.
    let mut weights = Vec::new();
    for block in fasta_map {
        match block.block_type {
            // The only place we want to insert variants.
            // Todo this is where we'd want to insert some calculated weights.
            Standard => {
                for _ in block.start..block.end {
                    weights.push(1);
                }
            }
            // no variants in N blocks
            NBlock => {
                for _ in block.start..block.end {
                    weights.push(0)
                }
            },
            // We won't insert overlapping variants.
            Variant => {
                for _ in block.start..block.end {
                    weights.push(0)
                }
            },
        }
    }
    debug_assert_eq!(fasta_map[fasta_map.len()-1].end, weights.len());
    weights
}