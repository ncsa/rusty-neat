// This library writes either single ended or paired-ended fastq files.

use std::io::Write;
use std::{fs, io};

use utils::fasta_tools::sequence_array_to_string;
use utils::file_tools::open_file;
use utils::neat_rng::NeatRng;
use utils::quality_scores::QualityScoreModel;

fn complement(nucleotide: u8) -> u8 {
    // 0 = A, 1 = C, 2 = G, 3 = T,
    // matches with the complement of each nucleotide.
    return match nucleotide {
        0 => 3,
        1 => 2,
        2 => 1,
        3 => 0,
        _ => 4,
    }
}

fn reverse_complement(sequence: &Vec<u8>) -> Vec<u8> {
    // Returns the reverse complement of a vector of u8's representing a DNA sequence.
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
    mut rng: &mut NeatRng,
) -> io::Result<()> {
    // Takes:
    // fastq_filename: prefix for the output fastq files.
    // paired_ended: boolean to set paired ended mode on or off.
    // dataset: List of u8 vectors representing dna sequences.
    // returns:
    // Error if there is a problem or else nothing.
    //
    // Writes fastq files. At the moment, it only writes out single r1 file, but will eventually write
    // out r1 and r2 files.

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
    // set up the quality score creator
    let quality_score_model = QualityScoreModel::new();
    // write sequences
    for sequence in dataset {
        // This assumes that the sequence length is the correct length at this point.
        let read_length = sequence.len();
        // Need to convert the raw scores to a string
        let quality_scores = quality_score_model.generate_quality_scores(
            read_length, &mut rng
        );
        // sequence name
        writeln!(&mut outfile1, "@{}/1", name_prefix.clone() + &index.to_string())?;
        // Array as a string
        writeln!(&mut outfile1, "{}", sequence_array_to_string(&sequence))?;
        // The stupid plus sign
        writeln!(&mut outfile1, "+")?;
        // Qual score of all F's for the whole thing.
        writeln!(&mut outfile1, "{}", quality_scores_to_str(quality_scores))?;
        index += 1;
        if paired_ended {
            // Need a quality score for this read as well
            let quality_scores = quality_score_model.generate_quality_scores(
                read_length, &mut rng
            );
            // sequence name
            writeln!(&mut outfile2, "@{}/2", name_prefix.clone() + &index.to_string())?;
            // Array as a string
            writeln!(&mut outfile2, "{}", sequence_array_to_string(&reverse_complement(&sequence)))?;
            // The stupid plus sign
            writeln!(&mut outfile2, "+")?;
            // Qual score of all F's for the whole thing.
            writeln!(&mut outfile2, "{}", quality_scores_to_str(quality_scores))?;
        }
    };
    if !paired_ended {
        fs::remove_file(filename2)?;
    }
    Ok(())
}

fn quality_scores_to_str(array: Vec<usize>) -> String {
    let mut score_text = String::new();
    for score in array {
        score_text += &(((score + 33) as u8) as char).to_string();
    }
    score_text
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;
    use rand_core::SeedableRng;

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
        let mut rng = NeatRng::seed_from_u64(0);
        let dataset = vec![&seq1, &seq2];
        write_fastq(
            fastq_filename,
            overwrite_output,
            paired_ended,
            dataset,
            &mut rng,
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
        let mut rng = NeatRng::seed_from_u64(0);
        let dataset = vec![&seq1, &seq2];
        write_fastq(
            fastq_filename,
            overwrite_output,
            paired_ended,
            dataset,
            &mut rng,
        ).unwrap();
        let outfile1 = Path::new("test_paired_r1.fastq");
        let outfile2 = Path::new("test_paired_r2.fastq");
        assert!(outfile1.exists());
        assert!(outfile2.exists());
        fs::remove_file(outfile1).unwrap();
        fs::remove_file(outfile2).unwrap();
    }
}