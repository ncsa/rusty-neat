use std::io::Write;
use std::{fs, io};

use utils::fasta_tools::sequence_array_to_string;
use utils::file_tools::open_file;

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
    overwrite_output: bool,
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
    let mut filename1 = String::from(fastq_filename) + "_r1.fastq";
    // open the file and append lines
    let mut outfile1 = open_file(&mut filename1, overwrite_output)
        .expect(&format!("Error opening output {}", filename1));
    // setting up pairend ended reads For single ended reads, this will go unused.
    let mut filename2 = String::from(fastq_filename) + "_r2.fastq";
    // open the second file and append lines
    let mut outfile2 = open_file(&mut filename2, overwrite_output)
        .expect(&format!("Error opening output {}", filename2));
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_complement() {
        let nuc1 = 0;
        let nuc2 = 1;
        let nuc3 = 2;
        let nuc4 = 3;
        let nuc5 = 4;

        assert_eq!(complement(nuc1), 3);
        assert_eq!(complement(nuc2), 2);
        assert_eq!(complement(nuc3), 1);
        assert_eq!(complement(nuc4), 0);
        assert_eq!(complement(nuc5), 4);
    }

    #[test]
    fn test_reverse_complement() {
        let read: Vec<u8> = vec![0, 0, 0, 0, 1, 1, 1, 1];
        let revcomp: Vec<u8> = vec![2, 2, 2, 2, 3, 3, 3, 3];
        assert_eq!(reverse_complement(&read), revcomp);
    }

    #[test]
    fn test_write_fastq_single() {
        let fastq_filename = "test_single";
        let overwrite_output = true;
        let paired_ended = false;
        let seq1 = vec![0, 0, 0, 0, 1, 1, 1, 1];
        let seq2 = vec![2, 2, 2, 2, 3, 3, 3, 3,];
        let dataset = vec![&seq1, &seq2];
        write_fastq(
            fastq_filename,
            overwrite_output,
            paired_ended,
            dataset
        ).unwrap();
        let outfile1 = Path::new("test_single_r1.fastq");
        let outfile2 = Path::new("test_single_r2.fastq");
        assert!(outfile1.exists());
        assert!(!outfile2.exists());
        fs::remove_file(outfile1).unwrap();
    }

    #[test]
    fn test_write_fastq_paired() {
        let fastq_filename = "test_paired";
        // might as well test the o_o function as well
        let overwrite_output = false;
        let paired_ended = true;
        let seq1 = vec![0, 0, 0, 0, 1, 1, 1, 1];
        let seq2 = vec![2, 2, 2, 2, 3, 3, 3, 3,];
        let dataset = vec![&seq1, &seq2];
        write_fastq(
            fastq_filename,
            overwrite_output,
            paired_ended,
            dataset
        ).unwrap();
        let outfile1 = Path::new("test_paired_r1.fastq");
        let outfile2 = Path::new("test_paired_r2.fastq");
        assert!(outfile1.exists());
        assert!(outfile2.exists());
        fs::remove_file(outfile1).unwrap();
        fs::remove_file(outfile2).unwrap();
    }
}