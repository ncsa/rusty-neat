use std::collections::HashMap;
use itertools::Itertools;
use rand::distributions::{WeightedIndex, Distribution};
use common::models::mutation_model::MutationModel;
use common::structs::variants::Variant;
use common::structs::nucleotides::Nuc;

pub fn generate_variants(
    reference_sequence: &Vec<Nuc>,
    mutation_model: &MutationModel,
    ploidy: usize,
) -> HashMap<usize, Variant> {
    // Variants for this contig, organized by location. Each location may have only one variant.
    // And we may need to implement a further check that there aren't other types of overlaps.
    // This will simply overwrite any existing value, so we just capture the last variant at each
    // location.
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
        let location = locations[dist.sample(mutation_model.get_mut_rng())];
        let mutation: Variant = mutation_model.generate_mutation(
            reference_sequence,
            location,
            ploidy,
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