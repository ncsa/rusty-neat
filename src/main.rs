#![allow(dead_code)]
mod logger;
use rusty_neat::{Config};
use std::{env, process};
use rand::{thread_rng};
use rand::prelude::*;

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
        match char {
            'A' => (),
            'C' => (),
            'G' => (),
            'T' => (),
            _ => { check_vec[index] = false },
        }
    }
    check_vec
}

fn generate_new_nucleotide(nuc: &char) -> char {
    let mut rng = thread_rng();
    let nucleotides = ['A', 'C', 'G', 'T'];
    let mut new = *nuc;
    while new == *nuc {
        new = nucleotides[rng.gen_range(0..4)];
    }
    new
}

fn main() {
    logger::init();
    logger::info!("Begin processing");

    // Argument parser
    let args: Vec<String> = env::args().collect();
    let config = Config::build(&args).unwrap_or_else(|err| {
        println!("Problem parsing arguments: {err}");
        process::exit(1)
    });

    println!("In file {}", config.fasta_path);

    if let Err(e) = rusty_neat::run(config) {
        println!("Application error: {e}");
        process::exit(1);
    }

    let sequence = String::from("ACATCTACACGGAGCACGACGCACCNNNNNNccccNNNNCACTGCAGTGCATGACG");
    logger::info!("Input sequence: {}", sequence);
    let non_n_map = map_non_n_regions(&sequence);
    if !non_n_map.contains(&true) {
        logger::info!("Input sequence contained no non-n values. Quitting.");
        std::process::exit(0)
    }
    logger::debug!("Non-n-map = {:?}", non_n_map);
    let position = find_random_non_n(&non_n_map).unwrap();
    logger::info!("Found a position! It is {}",position);
    let code = &sequence.as_bytes()[position];
    let nucleotide: char = *code as char;
    logger::info!("Nucleotide to mutate: {:?}", nucleotide);
    let new_nucleotide = generate_new_nucleotide(&nucleotide);
    let result: String = sequence.chars().enumerate().map(|(i, c)| if i == position { new_nucleotide } else { c }).collect();
    logger::info!("New sequence = {}", result);
    logger::info!("End program");
}

