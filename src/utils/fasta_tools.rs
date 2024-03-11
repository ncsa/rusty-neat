extern crate rand;
extern crate log;
extern crate assert_fs;

use std::collections::HashMap;
use crate::utils::file_tools::read_lines;
use self::log::info;
use std::fs::File;
use std::{fs, io};
use std::io::Write;
use std::*;
use std::os::unix::fs::MetadataExt;

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

#[repr(u8)]
#[allow(unused)]
pub enum NucRepr {
    // I just discovered this repr type thing and need to investigate.
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    N = 4,
}

pub fn sequence_array_to_string(input_array: &Vec<u8>) -> String {
    let mut return_string = String::new();
    for num in input_array {
        return_string += num_to_char(*num);
    }
    return_string
}

pub fn read_fasta(
    fasta_path: &str
) -> (Box<HashMap<String, Vec<u8>>>, Vec<String>) {
    // Reads a fasta and turns it into a HashMap and puts it in the heap
    info!("Reading fasta: {}", fasta_path);

    let mut fasta_map: HashMap<String, Vec<u8>> = HashMap::new();
    let mut fasta_order: Vec<String> = Vec::new();
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
                    fasta_order.push(current_key.clone());
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
    (Box::new(fasta_map), fasta_order)
}

pub fn write_fasta(
    fasta_output: &Box<HashMap<String, Vec<Vec<u8>>>>,
    fasta_order: &Vec<String>,
    output_file: &str,
    ploidy: usize
) -> io::Result<()> {
    // writing fasta output to files
    let mut fasta_files: Vec<String> = Vec::new();
    for ploid in 0..ploidy {
        let this_fasta = format!("{}_p{}.fasta", output_file, ploid+1);
        fasta_files.push(this_fasta.clone());
        info!("Writing file: {}", this_fasta);
        let mut outfile = File::options().create(true).append(true).open(this_fasta)?;
        for contig in fasta_order {
            let sequences = &fasta_output[contig];
            // Write contig name
            writeln!(&mut outfile, ">{}", contig)?;
            // write sequences[ploid] to this_fasta
            let mut i = 0;
            let sequence_to_write: &Vec<u8> = &sequences[ploid];
            while i < sequence_to_write.len() {
                let mut line = String::new();
                let mut max: usize = 70;
                let this_length: usize = sequence_to_write[i..].len();
                // If we don't have 70 characters, write out whatever is left.
                if this_length < 70 {
                    max = this_length
                }
                for j in 0..max {
                    line += num_to_char(sequence_to_write[i + j]);
                }
                writeln!(&mut outfile, "{}", line)?;
                i += 70;
            }
        }
    };
    Ok(())
}

#[test]
fn test_conversions() {
    let initial_sequence = "AAAANNNNGGGGCCCCTTTTAAAA";
    let test_map: Vec<u8> = vec![0,0,0,0,4,4,4,4,2,2,2,2,1,1,1,1,3,3,3,3,0,0,0,0];
    let remap: Vec<u8> = initial_sequence.chars().map(|x| char_to_num(x)).collect();
    assert_eq!(remap, test_map);
    assert_eq!(sequence_array_to_string(&test_map), initial_sequence);
}

#[test]
fn test_read_fasta() {
    let test_fasta = "data/H1N1.fa";
    let (_test_map, map_order) = read_fasta(test_fasta);
    assert_eq!(map_order[0], "H1N1_HA".to_string())
}

#[test]
fn test_write_fasta() -> Result<(), Box<dyn error::Error>> {

    let seq1: Vec<u8> = vec![0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1];
    let seq2: Vec<u8> = vec![0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,2];
    let sequences = vec![seq1, seq2];
    let fasta_output: HashMap<String, Vec<Vec<u8>>> = HashMap::from([
        (String::from("H1N1_HA"), sequences)
    ]);
    let fasta_pointer = Box::new(fasta_output);
    let fasta_order = vec![String::from("H1N1_HA")];
    let output_file = "test";
    let ploidy: usize = 2;
    let test_write = write_fasta(
        &fasta_pointer,
        &fasta_order,
        output_file,
        ploidy
    );
    let file_name_p1 = "test_p1.fasta";
    let file_name_p2 = "test_p2.fasta";
    assert_eq!(test_write.unwrap(), ());
    let attr = fs::metadata(file_name_p1).unwrap();
    assert!(attr.size() > 0);
    let attr = fs::metadata(file_name_p2).unwrap();
    assert!(attr.size() > 0);
    fs::remove_file(file_name_p1)?;
    fs::remove_file(file_name_p2)?;
    Ok(())
}