use std::collections::HashMap;
use itertools::Itertools;
use rand::distributions::{WeightedIndex, Distribution};
use rand_chacha::ChaCha20Rng;
use common::models::mutation_model::MutationModel;
use common::structs::variants::Variant;
use common::structs::nucleotides::Nuc;

pub fn generate_variants(
    reference_sequence: &Vec<Nuc>,
    mutation_model: &mut MutationModel,
    ploidy: usize,
    mut rng: ChaCha20Rng,
) -> HashMap<usize, Variant> {
    // Variants for this contig, organized by location. Each location may have only one variant.
    // And we may need to implement a further check that there aren't other types of overlaps.
    // This will simply overwrite any existing value, so we just capture the last variant at each
    // location.

    // todo We need to incorporate this idea

    // let non_n_positions: Vec<usize> = mutated_record
    //     .iter()
    //     .enumerate()
    //     .filter(|&(_, y)| *y != 4) // Filter out the N's
    //     .map(|(x, _)| x)
    //     .collect();
    // We don't want variants in our N zones, I think

    let reference_length = reference_sequence.len();
    let mut contig_variants: HashMap<usize, Variant> = HashMap::new();
    let mut number_of_mutations = (reference_length as f64 * mutation_model.mutation_rate)
        .round() as usize;
    let locations: Vec<usize> = (0..reference_length).collect();
    // Todo: Only look in non-n sections of the reference
    let weights: Vec<usize> = vec![1; reference_length];
    // Todo: Add weights due to gc-bias and trinucleotide bias
    let dist = WeightedIndex::new(&weights).unwrap();
    while number_of_mutations > 0 {
        let location = locations[dist.sample(&mut rng)];
        let mutation: Variant = mutation_model.generate_mutation(
            reference_sequence,
            location,
            ploidy,
            rng.clone(),
        );
        if contig_variants.keys().contains(&location) {
            // We'll just throw out the previous one and use this one.
            *contig_variants.get_mut(&location).unwrap() = mutation;
        } else {
            contig_variants.insert(location, mutation);
        }
        // For now, we'll count duplicates just to get through it.
        // basically this means we successfully added a mutation.
        number_of_mutations -= 1;
    }
    contig_variants
}