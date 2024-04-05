// This library writes either single ended or paired-ended fastq files.
use common;
use crate::utils;

use std::io::Write;
use std::{fs, io};
use std::collections::HashMap;
use rand::seq::SliceRandom;
use rand_chacha::ChaCha20Rng;

use common::structs::nucleotides::sequence_array_to_string;
use common::file_tools::open_file;
use common::structs::nucleotides::Nuc;
use common::models::quality_scores::QualityScoreModel;
use common::structs::variants::Variant;
use utils::apply_variants::apply_variants;

fn reverse_complement(sequence: &[Nuc]) -> Vec<Nuc> {
    // Returns the reverse complement of a vector of u8's representing a DNA sequence.
    let length = sequence.len();
    let mut rev_comp = Vec::new();
    for i in (0..length).rev() {
        rev_comp.push(sequence[i].complement())
    }
    rev_comp
}

pub fn write_fastq(
    mutated_fasta_map: &Box<HashMap<String, Vec<Nuc>>>,
    reads_dataset: &mut Box<Vec<(String, usize, usize)>>,
    read_length: usize,
    overwrite_output: bool,
    fastq_filename: &str,
    paired_ended: bool,
    quality_score_model: QualityScoreModel,
    mut rng: ChaCha20Rng,
) -> io::Result<()> {
    // Takes:
    // mutated_fasta_map: The sequences with mutations applied.
    // reads_dataset: coordinates for the reads.
    // read_length: How long a read to grab
    // overwrite_output: fail if output exists (false) else overwrite (true)
    // fastq_filename: the prefix for the fastq file(s) to write
    // paired_ended: boolean to set paired ended mode on or off.
    // quality_score_model: The statistical model we will use to generate quality scores
    // todo: sequencing_error_model: The statistical model we'll need to generate read errors.
    // rng: the random number generator for the run
    //
    // returns:
    // Error if there is a problem or else it successfully wrote the file(s)
    let read_name_prefix = "neat_generated_".to_string();
    let mut filename1 = String::from(fastq_filename) + "_r1.fastq";
    // open the file and append lines
    let mut outfile1 = open_file(&mut filename1, overwrite_output)
        .expect(&format!("Error opening output {}", filename1));
    // Setting up paired ended reads. For single ended reads, this will go unused.
    let mut filename2 = String::from(fastq_filename) + "_r2.fastq";
    // open the second file and append lines
    let mut outfile2 = open_file(&mut filename2, overwrite_output)
        .expect(&format!("Error opening output {}", filename2));
    // write out sequence by sequence using index to track the numbering
    let mut index = 1;
    reads_dataset.shuffle(&mut rng);
    for (contig, position1, position2) in reads_dataset.iter() {
        let quality_scores = quality_score_model
            .generate_quality_scores(read_length, &mut rng);
        // Todo apply sequencing errors
        let fragment = mutated_fasta_map
            .get(contig)
            .expect(&format!("Invalid contig: {contig}"))
            .get(*position1..*position2)
            .unwrap(); // Trim to length of initial fragment

        // sequence name
        writeln!(
            &mut outfile1,
            "@{}/1",
            read_name_prefix.clone() + &index.to_string()
        )?;
        // Array as a string
        writeln!(
            &mut outfile1,
            // First read is the first read_length amount of the fragment
            "{}", sequence_array_to_string(&fragment.get(..read_length).unwrap())
        )?;
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
            // The second read is the reverse complement at the end of the fragment
            let second_read = &fragment.get((fragment.len() - read_length)..).unwrap();
            // sequence name
            writeln!(
                &mut outfile2,
                "@{}/2",
                read_name_prefix.clone() + &index.to_string()
            )?;
            // Array as a string
            writeln!(
                &mut outfile2,
                "{}",
                sequence_array_to_string(&reverse_complement(second_read))
            )?;
            // The stupid plus sign
            writeln!(&mut outfile2, "+")?;
            // Qual score of all F's for the whole thing.
            writeln!(&mut outfile2, "{}", quality_scores_to_str(quality_scores))?;
        }
    }
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
    use rand_core::SeedableRng;
    use std::path::Path;
    use common::structs::nucleotides::Nuc::*;
    use common::structs::variants::VariantType;

    #[test]
    fn test_reverse_complement() {
        let read: Vec<Nuc> = vec![A, A, A, A, C, C, C, C];
        let revcomp: Vec<Nuc> = vec![G, G, G, G, T, T, T, T];
        assert_eq!(reverse_complement(&read), revcomp);
    }

    #[test]
    fn test_write_fastq() {
        let mutated_reference_seq = vec![A, G, T, A, C, T, C, A, G, T, G, T, T, C, C, T];
        let mutated_fasta_map = Box::new(HashMap::from([
            ("chr1".to_string(), mutated_reference_seq)
        ]));
        let mut reads_dataset = Box::new(vec![
            ("chr1".to_string(), 0, 4),
            ("chr1".to_string(), 6, 10),
        ]);
        let read_length = 4;
        let overwrite_output = true;
        let fastq_filename = "test_single";
        let paired_ended = false;
        let rng = ChaCha20Rng::seed_from_u64(0);
        let quality_score_model = QualityScoreModel::new();
        write_fastq(
            &mutated_fasta_map,
            &mut reads_dataset,
            read_length,
            overwrite_output,
            fastq_filename,
            paired_ended,
            quality_score_model.clone(),
            rng.clone()
        )
        .unwrap();
        let outfile1 = Path::new("test_single_r1.fastq");
        let outfile2 = Path::new("test_single_r2.fastq");
        assert!(outfile1.exists());
        assert!(!outfile2.exists());
        fs::remove_file(outfile1).unwrap();

        let mut reads_dataset = Box::new(vec![
            ("chr1".to_string(), 0, 4),
            ("chr1".to_string(), 6, 14),
        ]);

        let read_length = 4;
        let overwrite_output = false;
        let fastq_filename = "test_paired";
        let paired_ended = true;

        write_fastq(
            &mutated_fasta_map,
            &mut reads_dataset,
            read_length,
            overwrite_output,
            fastq_filename,
            paired_ended,
            quality_score_model,
            rng.clone()
        )
            .unwrap();

        let outfile1 = Path::new("test_paired_r1.fastq");
        let outfile2 = Path::new("test_paired_r2.fastq");
        assert!(outfile1.exists());
        assert!(outfile2.exists());
        fs::remove_file(outfile1).unwrap();
        fs::remove_file(outfile2).unwrap();
    }

}
