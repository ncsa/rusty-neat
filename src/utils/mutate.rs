extern crate log;
extern crate itertools;

use std::collections::HashMap;
use rand::distributions::Distribution;
use rand::prelude::{ThreadRng, IndexedRandom};
use rand::seq::SliceRandom;
use rand::Rng;
use self::log::debug;
use self::itertools::izip;

#[derive(Debug, Ord, PartialOrd, Eq, PartialEq)]
struct NucModel {
    // This simple nucleotide model simple tracks the base to mutate from
    // and the weights of each are turned into a vector. For example, if the weight count for "A"
    // is 8, then the vector will be [0, 0, 0, 0, 0, 0, 0, 0].
    // 0: A, 1: C, 2: G, 3: T
    base: u8,
    // vector of 0s, the length of which is the weight of the A in the vector.
    a: Vec<u8>,
    // vector of 1s, the length of which is the weight of the C in the vector.
    c: Vec<u8>,
    // vector of 2s, the length of which is the weight of the G in the vector.
    g: Vec<u8>,
    // vector of 3s, the length of which is the weight of the T in the vector.
    t: Vec<u8>,
}

impl NucModel {
    fn from(base: u8, weights: Vec<usize>) -> Self {
        Self {
            base,
            a: vec![0; weights[0]],
            c: vec![1; weights[1]],
            g: vec![2; weights[2]],
            t: vec![3; weights[3]],
        }
    }

    fn choose_new_nuc(&self, mut rng: &mut ThreadRng) -> u8 {
        // concatenate the vectors together
        let mut choices = [&self.a[..], &self.c[..], &self.g[..], &self.t[..]].concat();
        // shuffle the weighted list, based on whatever statistics we come up with.
        choices.shuffle(&mut rng);
        // start with the first choice
        let mut i = 0;
        // Look for the first nuc that isn't the old nuc.
        while choices[i] == self.base {
            i += 1;
        }
        choices[i]
    }
}

pub fn mutate_fasta(
    file_struct: &HashMap<String, Vec<u8>>,
    ploidy: usize,
    mut rng: &mut ThreadRng
) -> Box<HashMap<String, Vec<Vec<u8>>>> {
    // This takes a hashmap of contig names (keys) and a vector representing the reference sequence.
    // It performs a basic calculation (length x mutation rate) and chooses that many positions
    // along the sequence to mutate. It then builds a return string that represents the altered
    // sequence.

    // Todo: potentially re-write this to simple record the positions. Maybe position and length.
    // and which ploid. Need to do this once for each ploid.
    const MUT_RATE: f64 = 0.01; // will update this with something more elaborate later.
    // Generate random proportions to give the ploids, to keep things interesting.
    let mut proportions: Vec<f64> = Vec::with_capacity(ploidy);
    // For ploidy = N, we generate N random values, then divide them by the sum of the values, this
    // gives N random numbers whose sum = 1, to randomly assign proportion to the ploids. We may
    // find there is a better model later.
    for _ in 0..ploidy {
        proportions.push((&mut rng).gen::<f64>());
    }
    // Normalize my randos (divide by their total)
    let total: f64 = proportions.iter().sum();
    proportions = proportions.iter().map(|x| x/total).collect();
    // This leaves a list of numbers that add up to 1, since these are proportions of 1, will
    // ensure our mutation rate of each ploid adds up to the total MUT_RATE
    let ploid_mut_rate: Vec<f64> = proportions.iter().map(|x| x*MUT_RATE).collect();

    let mut return_struct: HashMap<String, Vec<Vec<u8>>> = HashMap::new(); // the mutated sequences
    for (name, sequence) in file_struct {
        let mut strands: Vec<Vec<u8>> = Vec::new();
        for ploid in 0..ploidy {
            let sequence_length = sequence.len();
            debug!("Sequence {} is {} bp long", name, sequence_length);

            // Clone the reference
            let mut mutated_record: Vec<u8> = sequence.clone();

            // Calculate how many mutations to add
            let num_positions = (
                sequence_length as f64 * ploid_mut_rate[ploid]
            ).round() as usize;
            if num_positions == 0 {
                if rng.gen_bool(1.0/(ploidy as f64)) == true {
                    mutated_record = mutate_sequence(
                        mutated_record, 1, sequence_length, &mut rng
                    );
                }
            } else {
                mutated_record = mutate_sequence(
                    mutated_record, num_positions, sequence_length, &mut rng
                );
            }
            strands.push(mutated_record.clone());
        }
        // Add to the new hashmap
        return_struct.entry(name.clone()).or_insert(strands.clone());
    }

    Box::new(return_struct)
}

fn mutate_sequence(
    sequence: Vec<u8>,
    num_positions: usize,
    length: usize,
    mut rng: &mut ThreadRng
) -> Vec<u8> {
    /*
    Takes a vector of u8's and mutate a few positions at random.
     */
    debug!("Adding {} mutations", num_positions);
    let mut mutated_record = sequence.clone();
    // Randomly select num_positions from positions, weighted by gc bias and whatever. For now
    // all he weights are just equal.
    let weights = vec![1; length];
    // zip together weights and positions
    let mut weighted_positions: Vec<(usize, i32)> = Vec::new();
    // izip!() accepts iterators and/or values with IntoIterator.
    for (x, y) in izip!(&mut (0..length), &weights) {
        weighted_positions.push((x, *y))
    }
    // now choose a random selection of num_positions without replacement
    let mutated_elements: Vec<&(usize, i32)> = weighted_positions
        .choose_multiple_weighted(&mut rng, num_positions, |x| x.1)
        .unwrap()
        .collect::<Vec<_>>();
    // Build mutation model
    let mut_a = NucModel::from(0, vec![0, 1, 1, 1]);
    let mut_c = NucModel::from(1, vec![1, 0, 1, 1]);
    let mut_g = NucModel::from(2, vec![1, 1, 0, 1]);
    let mut_t = NucModel::from(3, vec![1, 1, 1, 0]);

    for (index, _) in mutated_elements {
        mutated_record[*index] = match sequence[*index] {
            0 => mut_a.choose_new_nuc(&mut rng),
            1 => mut_c.choose_new_nuc(&mut rng),
            2 => mut_g.choose_new_nuc(&mut rng),
            3 => mut_t.choose_new_nuc(&mut rng),
            _ => 4
        }
    }
    mutated_record
}