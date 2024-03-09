extern crate log;

use std::collections::HashMap;
use rand::distributions::{Distribution, Uniform};
use rand::prelude::ThreadRng;
use rand::{Rng, thread_rng};
use self::log::{debug, error};

fn mutate_nucleotide(nucleotide: u8, rng: &mut ThreadRng) -> u8 {
    let mut new_nuc = nucleotide.clone();
    let mut nuc_choices = Uniform::new_inclusive(0, 3).unwrap();
    while new_nuc == nucleotide {
        new_nuc = nuc_choices.sample(rng);
    }
    new_nuc
}

pub fn mutate_fasta(file_struct: &HashMap<String, Vec<u8>>, ploidy: usize, mut rng: ThreadRng) -> Box<HashMap<String, Vec<Vec<u8>>>> {
    // This takes a hashmap of contig names (keys) and a vector representing the reference sequence.
    // It performs a basic calculation (length x mutation rate) and chooses that many positions
    // along the sequence to mutate. It then builds a return string that represents the altered
    // sequence.

    // Todo: potentially re-write this to simple record the positions. Maybe position and length.
    // and which ploid. Need to do this once for each ploid.
    const MUT_RATE: f64 = 0.01; // will update this with something more elaborate later.
    // Generate random proportions to give the ploids, to keep things interesting.
    let mut proportions: Vec<f64> = Vec::with_capacity(ploidy);
    (&mut rng).try_fill(&mut proportions[..]).unwrap();
    // Normalize my randos
    let total: f64 = proportions.iter().sum();
    proportions = proportions.iter().map(|x| x/total).collect(); // This leaves a list of numbers that add up to 1
    // since these are proportions of 1, will ensure our mutation rate of each ploid adds up to the total MUT_RATE
    let ploid_mut_rate: Vec<f64> = proportions.iter().map(|x| x*MUT_RATE).collect();

    let mut return_struct: HashMap<String, Vec<Vec<u8>>> = HashMap::new(); // the mutated sequences
    for (name, sequence) in file_struct {
        let mut strands: Vec<Vec<u8>> = Vec::new();
        for ploid in 0..ploidy {
            let sequence_length = sequence.len();
            debug!("Sequence {} is {} bp long", name, sequence_length);

            // Clone the reference
            let mut mutated_record: Vec<u8> = sequence.clone();
            let new = rng.gen::<f64>();
            debug!("Here's a random number, randomly: {}", new);

            // Calculate how many mutations to add
            let mut num_positions = (sequence_length as f64 * ploid_mut_rate[ploid]).round() as usize;
            if num_positions == 0 {
                if rng.gen_bool(1.0/(ploidy as f64)) == true {
                    mutated_record = mutate_sequence(mutated_record, 1, sequence_length, mut rng);
                }
            } else {
                mutated_record = mutate_sequence(mutated_record, num_positions, sequence_length, mut rng);
            }
            strands.push(mutated_record);
        // Add to the new hashmap
        return_struct.entry(name.clone()).or_insert(strands.clone());
        }
    }

    Box::new(return_struct)
}

fn mutate_sequence(sequence: Vec<u8>, num_positions: usize, length: usize, mut rng: ThreadRng) -> Vec<u8> {
    debug!("Adding {} mutations", num_positions);
    let mutated_record = sequence.clone();
    // This builds a uniform distribution across the range, curious how we can adapt this idea
    let position_range = Uniform::new(0, length).unwrap();
    // Grab a set of random positions
    // todo: check that this is without replacement
    let mutated_elements: Vec<usize> = position_range.sample_iter(&mut rng).take(num_positions).collect();
    let mut check = mutated_elements.clone();
    check.sort();
    check.dedup();
    if check.len() != mutated_elements.len() {
        error!("sampl_iter is pulling doubles on us!");
        panic!("Fix sample_iter")
    }
    // Debug check
    for index in mutated_elements {
        mutated_record[index] = mutate_nucleotide(sequence[index], &mut rng);
    }

    mutated_record
}