// This library writes either single ended or paired-ended fastq files.
use common;
use crate::utils;

use std::io::Write;
use std::{fs, io};
use std::cmp::min;
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
    fasta_map: &Box<HashMap<String, Vec<Nuc>>>,
    reads_dataset: &Box<HashMap<String, Vec<(usize, usize)>>>,
    mutations: &Box<HashMap<String, HashMap<usize, Variant>>>,
    read_length: usize,
    overwrite_output: bool,
    fastq_filename: &str,
    paired_ended: bool,
    quality_score_model: QualityScoreModel,
    mut rng: ChaCha20Rng,
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
    //
    // name_prefix is for the prefix for the read names. Reads are numbered in output order
    // (Although this feature is currently untested and unknown).
    // (May need sorting.)
    // todo, fold in the variants to the output reads,
    let name_prefix = "neat_generated_".to_string();
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
    // I think we need to produce the mutated fasta at this point.
    for (contig, sequence) in fasta_map.iter() {
        let contig_variants = mutations.get(contig).unwrap();
        // This is just lists of tuples of numbers, so it should be small in memory.
        let mut contig_reads: Vec<(usize, usize)> = reads_dataset.get(contig).unwrap().clone();
        contig_reads.shuffle(&mut rng);

        for (position1, position2) in contig_reads {
            // Prepping the read
            // Grab any relevant mutations
            let relevant_mutations  =
                contig_variants
                    .iter()
                    .filter(
                        |(location, _)|
                        (position1 < **location) && (**location < position2)
                    )
                    .collect::<HashMap<_, _>>()
                    .clone();
            // grabbing the raw sequence plus a buffer
            let end_plus_buffer = min(position2 + 50, sequence.len());
            let raw_sequence = sequence.get(position1..end_plus_buffer).unwrap();
            let mutated_sequence = apply_variants(
                raw_sequence, relevant_mutations, position1, rng.clone()
            );
            let quality_scores = quality_score_model
                .generate_quality_scores(read_length, &mut rng);
            // Todo apply sequencing errors
            let mutated_sequence = mutated_sequence
                .get(..(position2 - position1))
                .unwrap(); // Trim to length of initial fragment

            // sequence name
            writeln!(
                &mut outfile1,
                "@{}/1",
                name_prefix.clone() + &index.to_string()
            )?;
            // Array as a string
            writeln!(
                &mut outfile1,
                "{}", sequence_array_to_string(&mutated_sequence[..read_length])
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
                // The second read is at the end of the fragment
                let second_read = &mutated_sequence[(mutated_sequence.len() - read_length)..];
                // sequence name
                writeln!(
                    &mut outfile2,
                    "@{}/2",
                    name_prefix.clone() + &index.to_string()
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
        let reference_seq = vec![A, G, T, A, C, T, C, A, G, T, G, T, T, C, C, T];
        let fasta_map = Box::new(HashMap::from([
            ("chr1".to_string(), reference_seq)
        ]));
        let reads_dataset = Box::new(HashMap::from([
            ("chr1".to_string(), vec![(0, 4), (6, 10)])
        ]));
        let mutations = Box::new(HashMap::from([
            ("chr1".to_string(), HashMap::from([
                (0, Variant::new(VariantType::SNP, &vec![A], &vec![T],
                                 vec![0,1])),
                (7, Variant::new(VariantType::Indel, &vec![T],
                                 &vec![T, A, C], vec![1, 1])),
                (12, Variant::new(VariantType::SNP, &vec![A], &vec![T],
                                  vec![0,1]))
            ]))
        ]));
        let read_length = 4;
        let overwrite_output = true;
        let fastq_filename = "test_single";
        let paired_ended = false;
        let rng = ChaCha20Rng::seed_from_u64(0);
        let quality_score_model = QualityScoreModel::new();
        write_fastq(
            &fasta_map,
            &reads_dataset,
            &mutations,
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

        let reads_dataset = Box::new(HashMap::from([
            ("chr1".to_string(), vec![(0, 8), (6, 14)])
        ]));

        let read_length = 4;
        let overwrite_output = false;
        let fastq_filename = "test_paired";
        let paired_ended = true;

        write_fastq(
            &fasta_map,
            &reads_dataset,
            &mutations,
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
