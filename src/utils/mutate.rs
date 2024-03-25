// This is a basic mutation with SNPs using a basic mutation model.
// mutate_fasta takes a fasta Hashmap and returns a mutated version and the locations of the
// mutations introduced
//
// mutate_sequence adds actual mutations to the fasta sequence

use std::collections::HashMap;
use rand::prelude::*;
use rand::Rng;
use log::{debug, error};
use rand::distributions::WeightedIndex;
use utils::mutation_model::MutationModel;
use utils::neat_rng::NeatRng;
use utils::nucleotides::Nuc;

pub fn mutate_fasta(
    file_struct: &HashMap<String, Vec<Nuc>>,
    minimum_mutations: Option<usize>,
    mut rng: &mut NeatRng
) -> (Box<HashMap<String, Vec<Nuc>>>, Box<HashMap<String, Vec<(usize, Nuc, Nuc)>>>) {
    // Takes:
    // file_struct: a hashmap of contig names (keys) and a vector
    // representing the reference sequence.
    // minimum_mutations is a usize or None that indicates if there is a requested minimum.
    //      The default is for rusty-neat to allow 0 mutations.
    // ploidy: The number of copies of the genome within an organism's cells
    // rng: random number generator for the run
    //
    // Returns:
    // A tuple with pointers to:
    // A hashmap with keys that are contig names and a vector with the mutated sequence
    // A vector of tuples containing the location and alt of each variant
    //
    // This function performs a basic calculation (length x mutation rate +/- a random amount)
    // and chooses that many positions along the sequence to mutate. It then builds a return
    // string that represents the altered sequence and stores all the variants.
    const MUT_RATE: f64 = 0.01; // will update this with something more elaborate later.
    let mut return_struct: HashMap<String, Vec<Nuc>> = HashMap::new(); // the mutated sequences
    // hashmap with keys of the contig names with a list of positions and alts under the contig.
    let mut all_variants: HashMap<String, Vec<(usize, Nuc, Nuc)>> = HashMap::new();
    // For each sequence, figure out how many variants it should get and add them
    for (name, sequence) in file_struct {
        // The length of this sequence
        let sequence_length = sequence.len();
        debug!("Sequence {} is {} bp long", name, sequence_length);
        // Clone the reference to create mutations
        // Calculate how many mutations to add
        let mut rough_num_positions: f64 = sequence_length as f64 * MUT_RATE;
        // Add or subtract a few extra positions.
        rough_num_positions += {
            // A random amount up to 10% of the reads
            let factor: f64 = rng.gen::<f64>() * 0.10;
            // 25% of the time subtract, otherwise we'll add.
            let sign: f64 = if rng.gen_bool(0.25) { -1.0 } else { 1.0 };
            // add or subtract up to 10% of the reads.
            rough_num_positions * (sign * factor)
        };
        let rounded_num_positions = rough_num_positions.round() as usize;
        // Round the number of positions to the nearest usize.
        // If mininum_mutations have been entered, we'll use that, else we'll set that to 0.
        let mut num_positions = 0;
        if !minimum_mutations.is_none() {
            // if a minimum mutations value was entered, then that is the minimum per contig.
            if Some(rounded_num_positions) < minimum_mutations {
                num_positions = minimum_mutations.unwrap();
            } else {
                num_positions = rounded_num_positions;
            }
        } else {
            // Else 0 is our minimum
            if rough_num_positions.round() as usize > 0 {
                num_positions = rough_num_positions.round() as usize;
            }
        }
        // Mutates the sequence, using the original
        let (mutated_record, contig_mutations) = mutate_sequence(
            &sequence, num_positions, &mut rng
        );
        // Add to the return struct and variants map.
        return_struct.entry(name.clone()).or_insert(mutated_record.clone());
        all_variants.entry(name.clone()).or_insert(contig_mutations);
    }

    (Box::new(return_struct), Box::new(all_variants))
}

fn mutate_sequence(
    sequence: &Vec<Nuc>,
    num_positions: usize,
    mut rng: &mut NeatRng
) -> (Vec<Nuc>, Vec<(usize, Nuc, Nuc)>) {
    // Takes:
    // sequence: A u8 vector representing a sequence of DNA
    // num_positions: The number of mutations to add to this sequence
    // rng: random number generator for the run
    //
    // returns a tuple with:
    // Vec<Nuc> is the sequence itself
    // Vec(usize, u8, u8) is the position of the snp and the alt and ref alleles for that snp.
    //
    // Takes a vector of u8's and mutate a few positions at random. Returns the mutated sequence and
    // a list of tuples with the position and the alts of the SNPs.
    debug!("Adding {} mutations", num_positions);
    let mut mutated_record = sequence.clone();
    // find all non n positions.
    let non_n_positions: Vec<usize> = mutated_record
        .iter()
        .enumerate()
        .filter(|&(_, y)| *y != Nuc::N) // Filter out the N's
        .map(|(x, _)| x)
        .collect();
    // Randomly select num_positions from the non-N positions, weighted by gc and trinuc bias
    // (eventually). For now all he weights are just equal.
    let weights = vec![1; non_n_positions.len()];
    // create the distribution
    let dist = WeightedIndex::new(weights).unwrap();
    // now choose a random selection of num_positions without replacement
    let mut indexes_to_mutate: Vec<usize> = Vec::new();
    for _ in 0..num_positions {
        let pos = non_n_positions[dist.sample(&mut rng)];
        indexes_to_mutate.push(pos);
    }
    // Build the default mutation model
    // todo incorporate custom models
    let mutation_model = MutationModel::new();
    // Will hold the variants added to this sequence
    let mut sequence_variants: Vec<(usize, Nuc, Nuc)> = Vec::new();
    // for each index, picks a new base
    for index in indexes_to_mutate {
        // remember the reference for later.
        let reference_base = sequence[index];
        // pick a new base and assign the position to it.
        mutated_record[index] = mutation_model.mutate(reference_base, &mut rng);
        // This check simply ensures that our model actually mutated the base.
        debug_assert_ne!(mutated_record[index], reference_base);
        // add the location and alt base for the variant
        sequence_variants.push((index, mutated_record[index], reference_base))
    }
    (mutated_record, sequence_variants)
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use utils::neat_rng::NeatRng;
    use super::*;
    use utils::nucleotides::Nuc::*;

    #[test]
    fn test_mutate_sequence() {
        let seq1: Vec<Nuc> = vec![N, N, A, A, A, C, C, G, A, T, C, C, C];
        let num_positions = 2;
        let mut rng = NeatRng::seed_from_u64(0);
        let mutant = mutate_sequence(&seq1, num_positions, &mut rng);
        assert_eq!(mutant.0.len(), seq1.len());
        assert!(!mutant.1.is_empty());
        assert_eq!(mutant.0[0], N);
        assert_eq!(mutant.0[1], N);
    }

    #[test]
    fn test_mutate_fasta() {
        let seq: Vec<Nuc> = vec![N, N, A, A, A, C, C, G, A, T, C, C, C];
        let file_struct: HashMap<String, Vec<Nuc>> = HashMap::from([
                ("chr1".to_string(), seq.clone())
            ]);
        let mut rng = NeatRng::seed_from_u64(0);
        let mutations = mutate_fasta(
            &file_struct,
            Some(1),
            &mut rng,
        );
        assert!(mutations.0.contains_key("chr1"));
        assert!(mutations.1.contains_key("chr1"));
        let mutation_location = mutations.1["chr1"][0].0;
        let mutation_alt = mutations.1["chr1"][0].1;
        let mutation_ref = mutations.1["chr1"][0].2;
        assert_eq!(mutation_ref, seq[mutation_location]);
        assert_ne!(mutation_alt, mutation_ref)
    }

    #[test]
    fn test_mutate_fasta_no_mutations() {
        let seq: Vec<Nuc> = vec![N, N, A, A, A, C, C, G, A, T, C, C, C];
        let file_struct: HashMap<String, Vec<Nuc>> = HashMap::from([
            ("chr1".to_string(), seq.clone())
        ]);
        // if a random mutation suddenly pops up in a build, it's probably the seed for this.
        let mut rng = NeatRng::seed_from_u64(0);
        let mutations = mutate_fasta(
            &file_struct,
            None,
            &mut rng,
        );
        assert!(mutations.0.contains_key("chr1"));
        assert!(mutations.1.contains_key("chr1"));
        assert!(mutations.1["chr1"].is_empty());
    }
}