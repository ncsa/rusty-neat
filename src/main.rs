extern crate log;
mod logger;

use rand::prelude::IteratorRandom;
use rand::thread_rng;

fn find_random_non_n(positions: &Vec<bool>) -> Option<usize> {
    /*
    This function is intended to find a non-n position from a vector
    An improvement is would be to use the gc-bias as a weighting vector, path them both in
    and then use weighted sample.
     */
    let mut rng = thread_rng();
    let indexes = 0..positions.len();
    let mut index = indexes.choose_stable(&mut rng);
    while positions[index.unwrap()] == false {
        index = find_random_non_n(positions)
    };
    index
}

fn map_non_n_regions(sequence: &String) -> Vec<bool> {
    let mut check_vec: Vec<bool> = vec![true; sequence.len()];
    for (index, char) in sequence.chars().enumerate() {
        if char == 'N' {
            check_vec[index] = false;
        }
    }

    check_vec
}

fn main() {
    logger::init();
    logger::info!("Begin processing");
    let sequence = String::from("AACATACATACACATGAGCACGCGACTACATGAGACGTAGCA");
    logger::info!("Input sequence: {}", sequence);
    let non_n_map = map_non_n_regions(&sequence);
    logger::debug!("Non-n-map = {:?}", non_n_map);
    let position = find_random_non_n(&non_n_map).unwrap();
    logger::info!("Found a position! It is {}",position);
    let code = &sequence.as_bytes()[position];
    let nucleotide: char = *code as char;
    logger::info!("Nucleotide to mutate: {:?}", nucleotide);
    println!("End");
}

