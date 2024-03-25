// This library contains tools needed to process fasta files as input and output.

use log::info;
use std::io;
use std::io::Write;
use std::*;
use HashMap;
use utils::file_tools::read_lines;
use utils::file_tools::open_file;
use utils::nucleotides::Nuc;

pub fn sequence_array_to_string(input_array: &Vec<Nuc>) -> String {
    // Converts a sequence vector into a string representing the DNA sequence
    let mut return_string = String::new();
    for nuc in input_array {
        return_string += &nuc.nuc_to_char().to_string();
    }
    return_string
}

pub fn read_fasta(
    fasta_path: &str
) -> Result<(Box<HashMap<String, Vec<Nuc>>>, Vec<String>), io::Error> {
    // Reads a fasta file and turns it into a HashMap and puts it in the heap
    info!("Reading fasta: {}", fasta_path);

    let mut fasta_map: HashMap<String, Vec<Nuc>> = HashMap::new();
    let mut fasta_order: Vec<String> = Vec::new();
    let mut current_key = String::new();

    let lines = read_lines(fasta_path).unwrap();
    let mut temp_seq: Vec<Nuc> = vec![];
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
                    temp_seq.push(Nuc::base_to_nuc(char));
                }
            }
        },
        Err(error) => panic!("Problem reading fasta file: {}", error)
    });
    // Need to pick up the last one
    fasta_map.entry(current_key.clone()).or_insert(temp_seq.clone());
    Ok((Box::new(fasta_map), fasta_order))
}

pub fn write_fasta(
    fasta_output: &Box<HashMap<String, Vec<Nuc>>>,
    fasta_order: &Vec<String>,
    overwrite_output: bool,
    output_file: &str,
) -> io::Result<()> {
    /*
    Takes:
        fasta_output: the hashmap of mutated sequences with contig names
        fasta_order: A vector with the proper order for the fasta elements
        output_file: the prefix for the output file name
        ploidy: the number of copies of each chromosome the simulation is using
    Returns:
        Errors if there is a problem writing the file, otherwise it returns nothing.
     */
    // writing fasta output to files
    let mut output_fasta = format!("{}.fasta", output_file);
    let mut outfile = open_file(&mut output_fasta, overwrite_output)
        .expect(&format!("Error opening {}", output_fasta));
    for contig in fasta_order {
        let sequence = &fasta_output[contig];
        // Write contig name
        writeln!(&mut outfile, ">{}", contig)?;
        // write sequences[ploid] to this_fasta
        let mut i = 0;
        let sequence_to_write: &Vec<Nuc> = &sequence;
        while i < sequence_to_write.len() {
            let mut line = String::new();
            let mut max: usize = 70;
            let this_length: usize = sequence_to_write[i..].len();
            // If we don't have 70 characters, write out whatever is left.
            if this_length < 70 {
                max = this_length
            }
            for j in 0..max {
                line += &(sequence_to_write[i + j].nuc_to_char().to_string());
            }
            writeln!(&mut outfile, "{}", line)?;
            i += 70;
        }
    };
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use utils::nucleotides::Nuc::*;

    #[test]
    fn test_conversions() {
        let initial_sequence = "AAAANNNNGGGGCCCCTTTTAAAA";
        let test_map: Vec<Nuc> = vec![
            A, A, A, A, N, N, N, N, G, G, G, G, C, C, C, C, T, T, T, T, A, A, A, A
        ];
        let remap: Vec<Nuc> = initial_sequence.chars().map(|x| Nuc::base_to_nuc(x)).collect();
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
        let seq1: Vec<Nuc> = vec![A, A, A, A, A, A, A, A, C, C, C, C, C, C, C, C];
        let fasta_output: HashMap<String, Vec<Nuc>> = HashMap::from([
            (String::from("H1N1_HA"), seq1)
        ]);
        let fasta_pointer = Box::new(fasta_output);
        let fasta_order = vec![String::from("H1N1_HA")];
        let output_file = "test";
        let test_write = write_fasta(
            &fasta_pointer,
            &fasta_order,
            true,
            output_file
        ).unwrap();
        let file_name = "test.fasta";
        assert_eq!(test_write, ());
        let attr = fs::metadata(file_name).unwrap();
        assert!(attr.len() > 0);
        fs::remove_file(file_name)?;
        Ok(())
    }
}