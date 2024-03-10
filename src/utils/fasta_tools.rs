extern crate rand;
extern crate log;

use rand::prelude::{IteratorRandom, ThreadRng};
use std::collections::HashMap;
use serde_yaml::Value::String;
use crate::utils::file_tools::read_lines;
use self::log::info;
use std::fs::File;

pub fn read_fasta(fasta_path: &str) -> Box<HashMap<String, Vec<u8>>> {
    info!("Reading fasta: {}", fasta_path);

    let mut fasta_map: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_key = String::new();

    if let Ok(lines) = read_lines(fasta_path) {
        let mut temp_seq: Vec<u8> = vec![];
        lines.for_each(|line| match line {
            Ok(l) => {
                if l.starts_with('>') {
                    if !current_key.is_empty() {
                        fasta_map.entry(current_key.clone()).or_insert(temp_seq.clone());
                    }
                    current_key = String::from(l.strip_prefix('>').unwrap());
                    temp_seq = vec![];
                } else {
                    for char in l.chars() {
                        temp_seq.push(char_to_num(char));
                    }
                }
            },
            Err(_) => (),
        });
        // Need to pick up the last one
        fasta_map.entry(current_key.clone()).or_insert(temp_seq.clone());
    }
    Box::new(fasta_map)
}

pub fn char_to_num(char_of_interest: char) -> u8 {
    return match char_of_interest {
        'A' | 'a' => 0,
        'C' | 'c' => 1,
        'G' | 'g' => 2,
        'T' | 't' => 3,
        _ => 4
    }
}

pub fn num_to_char(nuc_num: u8) -> &'static str {
    return match nuc_num {
        0 => "A",
        1 => "C",
        2 => "G",
        3 => "T",
        _ => "N",
    }
}

pub fn sequence_array_to_string(input_array: &Vec<u8>) -> String {
    let mut return_string = String::new();
    for num in input_array {
        return_string += num_to_char(*num);
    }
    return_string
}

pub fn map_non_n_regions(sequence: &Vec<u8>) -> Vec<bool> {
    let mut check_vec: Vec<bool> = vec![true; sequence.len()];
    for (index, nuc) in sequence.iter().enumerate() {
        match nuc {
            4 => { check_vec[index] = false; },
            _ => continue,
        }
    }

    check_vec
}

pub fn find_random_non_n(positions: &Vec<bool>, rng: &mut ThreadRng) -> Option<usize> {
    /*
    This function is intended to find a non-n position from a vector
    An improvement is would be to use the gc-bias as a weighting vector, path them both in
    and then use weighted sample.
     */
    let indexes = 0..positions.len();
    let mut index = indexes.choose_stable(rng);
    while positions[index.unwrap()] == false {
        index = find_random_non_n(positions, rng);
    };
    index
}

pub fn write_fasta(fasta_output: &Box<HashMap<String, Vec<Vec<u8>>>>, output_file: &str, ploidy: usize) {
    // setup files to write.
    let mut fasta_files: Vec<String> = Vec::new();
    for ploid in 0..ploidy {
        let this_fasta = format!("{}p{}.fasta", output_file, ploid+1);
        fasta_files.push(&this_fasta);
        let mut outfile = File::options().create(true).append(true).open(this_fasta);
        for (contig, sequences) in &fasta_output {
            // Write contig name
            writeln!(&mut outfile, ">{}", contig)?;
            // write sequences[ploid] to this_fasta
            let mut i = 0;
            let sequence_to_write = sequences[ploid];
            while i < sequence_to_write.len() {
                let line = String::new();
                let mut max: usize = 60;
                let this_length: usize = sequence_to_write[i..].len();
                if this_length < 60 {
                    max = this_length
                }
                for j in 0..max {
                    line += num_to_char(sequence_to_write[i + j]);
                }
                i += 60;
                writeln!(&mut outfile, "{}", line)
            }
        }
    }
}