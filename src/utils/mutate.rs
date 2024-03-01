extern crate log;

use std::collections::HashMap;
use rand::distributions::{Distribution, Uniform};
use rand::prelude::ThreadRng;
use rand::{Rng, thread_rng};
use self::log::debug;

fn mutate_nucleotide(nucleotide: u8, rng: &mut ThreadRng) -> u8 {
    let mut new_nuc = nucleotide.clone();
    let mut nuc_choices = Uniform::new_inclusive(0, 3).unwrap();
    while new_nuc == nucleotide {
        new_nuc = nuc_choices.sample(rng);
    }
    new_nuc
}

pub fn mutate_fasta(file_struct: &HashMap<String, Vec<u8>>) -> Box<HashMap<String, Vec<u8>>> {
    // This takes a hashmap of contig names (keys) and a vector representing the reference sequence.
    // It performs a basic calculation (length x mutation rate) and chooses that many positions
    // along the sequence to mutate. It then builds a return string that represents the altered
    // sequence.
    const MUT_RATE: f64 = 0.01; // will update this with something more elaborate later.
    let mut rng = thread_rng(); // basic random number generator
    let mut return_struct: HashMap<String, Vec<u8>> = HashMap::new(); // the mutated sequences
    for (name, sequence) in file_struct {
        let sequence_length = sequence.len();
        debug!("Sequence {} is {} bp long", name, sequence_length);

        // This builds a uniform distribution across the range, curious how we can adapt this idea
        let position_range = Uniform::new(0, sequence_length).unwrap();

        // Calculate how many mutations to add
        let num_positions = (sequence_length as f64 * MUT_RATE).round() as usize;
        debug!("Adding {} mutations", num_positions);

        // Start by cloning the reference
        let mut mutated_record = sequence.clone();
        // Grab a set of random positions
        // todo: check that this is without replacement
        let mutated_elements: Vec<usize> = position_range.sample_iter(&mut rng).take(num_positions).collect();
        // Debug check
        for index in mutated_elements {
            mutated_record[index] = mutate_nucleotide(sequence[index], &mut rng);
        }
        // Add to the new hashmap
        return_struct.entry(name.clone()).or_insert(mutated_record);
    }
    Box::new(return_struct)
}