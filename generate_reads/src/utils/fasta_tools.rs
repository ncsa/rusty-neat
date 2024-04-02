// This library contains tools needed to process fasta files as input and output.
use common;

use log::info;
use std::{io, str, usize};
use std::io::Write;
use std::*;
use std::collections::{HashMap, VecDeque};
use common::file_tools::open_file;
use common::file_tools::read_lines;
use common::models::mutation_model::MutationModel;
use common::structs::nucleotides::Nuc;
use common::structs::nucleotides::base_to_nuc;
use common::structs::variants::{Variant, VariantType};

pub fn read_fasta(
    fasta_path: &str,
) -> Result<(Box<HashMap<String, Vec<Nuc>>>, VecDeque<String>), io::Error> {
    // Reads a fasta file and turns it into a HashMap and puts it in the heap
    info!("Reading fasta: {}", fasta_path);

    let mut fasta_map: HashMap<String, Vec<Nuc>> = HashMap::new();
    let mut fasta_order: VecDeque<String> = VecDeque::new();
    let mut current_key = String::new();

    let lines = read_lines(fasta_path).unwrap();
    let mut temp_seq: Vec<Nuc> = vec![];
    lines.for_each(|line| match line {
        Ok(l) => {
            if l.starts_with('>') {
                if !current_key.is_empty() {
                    fasta_map
                        .entry(current_key.clone())
                        .or_insert(temp_seq.clone());
                }
                current_key = String::from(l.strip_prefix('>').unwrap());
                fasta_order.push_back(current_key.clone());
                temp_seq = vec![];
            } else {
                for char in l.chars() {
                    temp_seq.push(base_to_nuc(char));
                }
            }
        }
        Err(error) => panic!("Problem reading fasta file: {}", error),
    });
    // Need to pick up the last one
    fasta_map
        .entry(current_key.clone())
        .or_insert(temp_seq.clone());
    Ok((Box::new(fasta_map), fasta_order))
}

pub fn write_fasta(
    reference_fasta: &Box<HashMap<String, Vec<Nuc>>>,
    variants: &HashMap<String, HashMap<usize, Variant>>,
    fasta_order: &VecDeque<String>,
    overwrite_output: bool,
    output_file: &str,
) -> io::Result<()> {
    // reference_fasta: the hashmap of mutated sequences with contig names
    // variants: the location, type, and size of the variants that have been inserted in the fasta.
    // fasta_lengths: A vector with the proper order for the fasta elements and lengths of the
    //      sequences.
    // overwrite_output: Boolean that controls whether an existing file is overwritten.
    // output_file: the prefix for the output file name
    const LINE_LENGTH: usize = 70; // fasta is usually writes with 70 characters per line
    let mut output_fasta = format!("{}.fasta", output_file);
    let mut outfile = open_file(&mut output_fasta, overwrite_output)
        .expect(&format!("Error opening {}", output_fasta));
    for contig in fasta_order {
        let reference_sequence = &reference_fasta[contig];
        let contig_length = reference_sequence.len();
        // Write contig name
        writeln!(&mut outfile, ">{}", contig)?;
        // write sequences[ploid] to this_fasta
        let mut i = 0;
        let contig_variants = variants.get(contig).unwrap();
        let variant_locations: Vec<&usize> = contig_variants.keys().collect();
        while i < contig_length {
            let mut mutated_line = String::new();
            for j in 0..LINE_LENGTH {
                // Check to make sure we are not out of bounds, which will happen if we
                // reach the end and there are not exactly 70 bases left (basically, every time)
                if (i + j) >= contig_length {
                    // ensure that the outer loop breaks after we perform the last writeln.
                    i += j;
                    break
                }
                let base_to_write: &Nuc = reference_sequence.get(i).unwrap();
                if *base_to_write == Nuc::N || variant_locations.contains(&&(i + j)) {
                    let variant_to_apply = contig_variants.get(&(i+j)).unwrap();
                    mutated_line += &match variant_to_apply.variant_type {
                        VariantType::Indel => {
                            if variant_to_apply.is_insertion() { // insertion
                                let mut insertion_string = String::new();
                                for nuc in &variant_to_apply.alternate {
                                    insertion_string += &nuc.to_str();
                                }
                                insertion_string
                            } else { // Deletion
                                // The plus one for the base will get added after this if check
                                i += variant_to_apply.reference.len();
                                base_to_write.to_str().to_string()
                            }
                        },
                        VariantType::SNP => {
                            variant_to_apply.alternate[0].to_str().to_string()
                        },
                    }
                } else {
                    mutated_line += &base_to_write.to_str();
                }
                i += 1
            }
            writeln!(&mut outfile, "{}", mutated_line)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use common::structs::nucleotides::Nuc::*;
    use common::structs::nucleotides::base_to_nuc;
    use common::structs::nucleotides::sequence_array_to_string;

    #[test]
    fn test_conversions() {
        let initial_sequence = "AAAANNNNGGGGCCCCTTTTAAAA";
        let test_map: Vec<Nuc> = vec![
            A, A, A, A, N, N, N, N, G, G, G, G, C, C, C, C, T, T, T, T, A, A, A, A,
        ];
        let remap: Vec<Nuc> = initial_sequence
            .chars()
            .map(|x| base_to_nuc(x))
            .collect();
        assert_eq!(remap, test_map);
        assert_eq!(sequence_array_to_string(&test_map), initial_sequence);
    }

    #[test]
    fn test_read_fasta() {
        let test_fasta = "test_data/H1N1.fa";
        let (_test_map, map_order) = read_fasta(test_fasta).unwrap();
        assert_eq!(map_order[0], "H1N1_HA".to_string())
    }

    #[test]
    #[should_panic]
    fn test_read_bad_fasta() {
        let test_fasta = "test_data/fake.fasta";
        read_fasta(test_fasta).unwrap();
    }

    #[test]
    fn test_write_fasta() -> Result<(), Box<dyn error::Error>> {
        let reference_seq = vec![A, G, T, A, C, T, C, A, G, T, G, T, T, C, C, T];
        let fasta_map = Box::new(HashMap::from([
            ("chr1".to_string(), reference_seq.clone())
        ]));
        let reads_dataset = Box::new(HashMap::from([
            ("chr1".to_string(), vec![(0, 4), (6, 10)])
        ]));
        let mutations = Box::new(HashMap::from([
            ("chr1".to_string(), HashMap::from([
                (0, Variant::new(VariantType::SNP, 0, &vec![A], &vec![T],
                                 vec![0,1], false)),
                (7, Variant::new(VariantType::Indel, 7, &vec![T],
                                 &vec![T, A, C], vec![1, 1], true)),
                (12, Variant::new(VariantType::SNP, 12, &vec![A], &vec![T],
                                  vec![0,1], false))
            ]))
        ]));

        let fasta_output: HashMap<String, Vec<Nuc>> =
            HashMap::from([(String::from("chr1"), reference_seq)]);
        let fasta_pointer = Box::new(fasta_output);
        let fasta_order = VecDeque::from([String::from("chr1")]);
        let output_file = "test";
        let test_write = write_fasta(
            &fasta_map,
            &mutations,
            &fasta_order,
            false, output_file
        ).unwrap();
        let file_name = "test.fasta";
        assert_eq!(test_write, ());
        let attr = fs::metadata(file_name).unwrap();
        assert!(attr.len() > 0);
        fs::remove_file(file_name)?;
        Ok(())
    }
}
