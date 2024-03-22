// This is a basic mutation with SNPs using a basic mutation model.
// mutate_fasta takes a fasta Hashmap and returns a mutated version and the locations of the
// mutations introduced
//
// mutate_sequence adds actual mutations to the fasta sequence

use std::cmp::max;
use std::collections::HashMap;
use rand::prelude::IndexedRandom;
use rand::Rng;
use log::debug;
use itertools::izip;
use utils::nucleotides::NucModel;
use utils::neat_rng::NeatRng;

pub fn mutate_fasta(
    file_struct: &HashMap<String, Vec<u8>>,
    mut rng: &mut NeatRng
) -> (Box<HashMap<String, Vec<u8>>>, Box<HashMap<String, Vec<(usize, u8, u8)>>>) {
    /*
    Takes:
        file_struct: a hashmap of contig names (keys) and a vector
            representing the reference sequence.
        ploidy: The number of copies of the genome within an organism's cells
        rng: random number generator for the run

    Returns:
        A tuple with pointers to:
          A hashmap with keys that are contig names and a vector with the mutated sequence
          A vector of tuples containing the location and alt of each variant

    This function performs a basic calculation (length x mutation rate +/- a random amount)
    and chooses that many positions along the sequence to mutate. It then builds a return
    string that represents the altered sequence and stores all the variants.
     */
    const MUT_RATE: f64 = 0.01; // will update this with something more elaborate later.
    let mut return_struct: HashMap<String, Vec<u8>> = HashMap::new(); // the mutated sequences

    // hashmap with keys of the contig names with a list of positions and alts under the contig.
    let mut all_variants: HashMap<String, Vec<(usize, u8, u8)>> = HashMap::new();

    for (name, sequence) in file_struct {
        // Mutations for this contig
        let contig_mutations: Vec<(usize, u8, u8)>;
        // The length of this sequence
        let sequence_length = sequence.len();
        debug!("Sequence {} is {} bp long", name, sequence_length);
        // Clone the reference to create mutations
        let mut mutated_record: Vec<u8> = sequence.clone();
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
        // Round the number of positions to the nearest usize.
        // If negative or no reads, we still want at least 1 mutation per contig.
        let num_positions = max(1, rough_num_positions.round() as usize);

        (mutated_record, contig_mutations) = mutate_sequence(
            &mutated_record, num_positions, &mut rng
        );
        // Add to the return struct and variants map.
        return_struct.entry(name.clone()).or_insert(mutated_record.clone());
        all_variants.entry(name.clone()).or_insert(contig_mutations);
    }

    (Box::new(return_struct), Box::new(all_variants))
}

fn mutate_sequence(
    sequence: &Vec<u8>,
    num_positions: usize,
    mut rng: &mut NeatRng
) -> (Vec<u8>, Vec<(usize, u8, u8)>) {
    /*
    Takes:
        sequence: A u8 vector representing a sequence of DNA
        num_positions: The number of mutations to add to this sequence
        rng: random number generator for the run

    returns a tuple with:
        Vec<u8> = the sequence itself
        Vec(usize, u8) = the position of the snp and the alt allele for that snp.

    Takes a vector of u8's and mutate a few positions at random. Returns the mutated sequence and
    a list of tuples with the position and the alts of the SNPs.
     */
    debug!("Adding {} mutations", num_positions);
    let mut mutated_record = sequence.clone();
    // Randomly select num_positions from positions, weighted by gc bias and whatever. For now
    // all he weights are just equal.
    let weights = vec![1; mutated_record.len()];
    // zip together weights and positions
    let mut weighted_positions: Vec<(usize, i32)> = Vec::new();
    // find all non n positions.
    let non_n_positions: Vec<usize> = mutated_record
        .iter()
        .enumerate()
        .filter(|&(_, y)| *y != 4u8)
        .map(|(x, _)| x)
        .collect();
    // izip!() accepts iterators and/or values with IntoIterator.
    for (x, y) in izip!(&non_n_positions, &weights) {
        weighted_positions.push((*x, *y))
    }
    // now choose a random selection of num_positions without replacement
    let mutated_elements: Vec<&(usize, i32)> = weighted_positions
        .choose_multiple_weighted(&mut rng, num_positions, |x| x.1)
        .unwrap()
        .collect::<Vec<_>>();
    // Build the default mutation model
    // todo incorporate custom models
    let nucleotide_mutation_model = NucModel::new();

    // Will hold the variants added to this sequence
    let mut sequence_variants: Vec<(usize, u8, u8)> = Vec::new();
    // for each index, picks a new base
    for (index, _weight) in mutated_elements {
        // Skip the N's
        if sequence[*index] == 4 {
            continue
        }
        // remeber the reference for later.
        let reference_base = sequence[*index];
        // pick a new base and assign the position to it.
        mutated_record[*index] = nucleotide_mutation_model.choose_new_nuc(reference_base, &mut rng);
        // This check simply ensures that our model actually mutated the base.
        if mutated_record[*index] == reference_base {
            panic!("BUG: Mutation model failed to mutate the base. This should not happen.")
        }
        // add the location and alt base for the variant
        sequence_variants.push((*index, mutated_record[*index], reference_base))
    }
    (mutated_record, sequence_variants)
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use utils::neat_rng::NeatRng;
    use super::*;

    #[test]
    fn test_mutate_sequence() {
        let seq1: Vec<u8> = vec![4, 4, 0, 0, 0, 1, 1, 2, 0, 3, 1, 1, 1];
        let num_positions = 2;
        let mut rng = NeatRng::seed_from_u64(0);
        let mutant = mutate_sequence(&seq1, num_positions, &mut rng);
        assert_eq!(mutant.0.len(), seq1.len());
        assert!(!mutant.1.is_empty());
        // N's stay N's
        assert_eq!(mutant.0[0], 4);
        assert_eq!(mutant.0[1], 4);
    }

    #[test]
    fn test_mutate_fasta() {
        let seq = vec![4, 4, 0, 0, 0, 1, 1, 2, 0, 3, 1, 1, 1];
        let file_struct: HashMap<String, Vec<u8>> = HashMap::from([
                ("chr1".to_string(), seq.clone())
            ]);
        let mut rng = NeatRng::seed_from_u64(0);
        let mutations = mutate_fasta(
            &file_struct,
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
}