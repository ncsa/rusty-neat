use std::io::Write;
use std::fs::File;
use std::{fs, io};

use utils::fasta_tools::sequence_array_to_string;

fn complement(nucleotide: u8) -> u8 {
    /*
    0 = A, 1 = C, 2 = G, 3 = T,
    matches with the complement of each nucleotide.
     */
    return match nucleotide {
        0 => 3,
        1 => 2,
        2 => 1,
        3 => 0,
        _ => 4,
    }
}

fn reverse_complement(sequence: &Vec<u8>) -> Vec<u8> {
    /*
    Returns the reverse complement of a vector of u8's representing a DNA sequence.
     */
    let length = sequence.len();
    let mut rev_comp = Vec::new();
    for i in (0..length).rev() {
        rev_comp.push(complement(sequence[i]))
    }
    rev_comp
}

pub fn write_fastq(
    fastq_filename: &str,
    paired_ended: bool,
    dataset: Vec<&Vec<u8>>,
) -> io::Result<()> {
    /*
    Takes:
        fastq_filename: prefix for the output fastq files.
        paired_ended: boolean to set paired ended mode on or off.
        dataset: List of u8 vectors representing dna sequences.
    returns:
        Error if there is a problem or else nothing.

    Writes fastq files. At the moment, it only writes out single r1 file, but will eventually write
    out r1 and r2 files.
     */

    // name_prefix is for the prefix for the read names. Reads are numbered in output order
    // (Although this feature is currently untested and unknown).
    // (May need sorting.)
    let name_prefix = "neat_generated_".to_string();
    let filename1 = String::from(fastq_filename) + "_r1.fastq";
    // open the file and append lines
    let mut outfile1 = File::options().create_new(true).append(true).open(&filename1)?;
    // setting up pairend ended reads For single ended reads, this will go unused.
    let filename2 = String::from(fastq_filename) + "_r2.fastq";
    // open the second file and append lines
    let mut outfile2 = File::options().create_new(true).append(true).open(&filename2)?;
    // write out sequence by sequence using index to track the numbering
    let mut index = 1;
    for sequence in dataset {
        let sequence_len = sequence.len();
        // sequence name
        writeln!(&mut outfile1, "@{}/1", name_prefix.clone() + &index.to_string())?;
        // Array as a string
        writeln!(&mut outfile1, "{}", sequence_array_to_string(&sequence))?;
        // The stupid plus sign
        writeln!(&mut outfile1, "+")?;
        // Qual score of all F's for the whole thing.
        writeln!(&mut outfile1, "{}", std::iter::repeat("F").take(sequence_len).collect::<String>())?;
        index += 1;
        if paired_ended {
            // sequence name
            writeln!(&mut outfile2, "@{}/2", name_prefix.clone() + &index.to_string())?;
            // Array as a string
            writeln!(&mut outfile2, "{}", sequence_array_to_string(&reverse_complement(&sequence)))?;
            // The stupid plus sign
            writeln!(&mut outfile2, "+")?;
            // Qual score of all F's for the whole thing.
            writeln!(&mut outfile2, "{}", std::iter::repeat("F").take(sequence_len).collect::<String>())?;
        }
    };
    if !paired_ended {
        fs::remove_file(filename2)?;
    }
    Ok(())
}