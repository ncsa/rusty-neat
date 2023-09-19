mod logger;

use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Lines};
use std::path::Path;
use rand::prelude::IteratorRandom;
use rand::{Rng, thread_rng};

pub struct Config {
    pub fasta_path: String,
}

impl Config {
    pub fn build(args: &[String]) -> Result<Config, &'static str> {
        if args.len() < 2 {
            return Err("not enough arguments");
        }
        let file_path = args[1].clone();

        Ok(Config { fasta_path: file_path })
    }
}

pub fn run(config: Config) -> Result<(), Box<dyn Error>>{
    if let Ok(lines) = read_lines(config.fasta_path) {
        // Consumes the iterator, returns an (Optional) String
        process_lines(&lines);
    }


    // let sequence = String::from("ACATCTACACGGAGCACGACGCACCAGACACACTGCAGTGCATGACG");
    // logger::info!("Input sequence: {}", sequence);
    // let non_n_map = map_non_n_regions(&sequence);
    // if !non_n_map.contains(&true) {
    //     logger::info!("Input sequence contained no non-n values. Quitting.");
    //     std::process::exit(0)
    // }
    // logger::debug!("Non-n-map = {:?}", non_n_map);
    // let position = find_random_non_n(&non_n_map).unwrap();
    // logger::info!("Found a position! It is {}",position);
    // let code = &sequence.as_bytes()[position];
    // let nucleotide: char = *code as char;
    // logger::info!("Nucleotide to mutate: {:?}", nucleotide);
    // let new_nucleotide = generate_new_nucleotide(&nucleotide);
    // let result: String = sequence.chars().enumerate().map(|(i, c)| if i == position { new_nucleotide } else { c }).collect();
    // logger::info!("New sequence = {}", result);
    // logger::info!("End program");
    Ok(())
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn process_lines(&lines: Lines<BufReader<File>>) -> HashMap<String, String> {
    let mut file_map: HashMap<String, String> = HashMap::new();
    let mut sequence: String = String::new();
    let mut prev_key = "";
    let mut current_key = "";

    for line in lines {
        if let Ok(ip) = line {
            if ip.starts_with(">") {
                // strip off the > and the \n
                ip.retain(|c| !r#"(),",;:'"#.contains(c))
            }
        }
    }

        if stripped_line.starts_with(">") {
            let mut line_chars = stripped_line.chars();
            line_chars.next();
            let new_key = line_chars.as_str();
            if prev_key.is_empty() {
                if current_key.is_empty() {
                    current_key = &new_key;
                } else {
                    file_map.insert(current_key.to_string(), sequence);
                    prev_key = &current_key;
                    current_key = &new_key;
                    sequence = String::new();
                }
            }
        } else {
            sequence += &stripped_line;
        }
    }
    file_map
}

#[derive(Debug, PartialEq, PartialOrd)]
pub enum Nucleotide {
    A,
    C,
    G,
    T
}

#[derive(Debug)]
pub struct Base;

pub trait Reference {
    fn get_reference(
        self: &Self,
        reference: Nucleotide,
    ) -> Result<usize, String>;
}

impl Reference for Base {
    fn get_reference(self: &Self, reference: Nucleotide) -> Result<usize, String> {
        Ok(0)
    }
}

pub fn find_random_non_n(positions: &Vec<bool>) -> Option<usize> {
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

pub fn map_non_n_regions(sequence: &String) -> Vec<bool> {
    let mut check_vec: Vec<bool> = vec![true; sequence.len()];
    for (index, char) in sequence.chars().enumerate() {
        if char == 'N' {
            check_vec[index] = false;
        }
    }

    check_vec
}

pub fn generate_new_nucleotide(nuc: &char) -> char {
    let mut rng = thread_rng();
    let nucleotides = ['A', 'C', 'G', 'T'];
    let mut new = *nuc;
    while new == *nuc {
        new = nucleotides[rng.gen_range(0..4)];
    }
    new
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_result() {
        let fasta_path = "C:\\Users\\joshf\\data\\H1N1.fa";
        let file_map = parse_file(fasta_path);
        let mut first_key = String::new();
        for key in file_map.keys() {
            first_key = String::from(key);
            break;
        }

        assert_eq!("H1N1_HA".to_string(), first_key.to_string());
    }
}