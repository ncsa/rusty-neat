//! This library writes either single ended or paired-ended fastq files.
//! Need to update this method. We want to use the data structures and we want to make sure
//! this function is generic enough to work with the fragmented method we are implementing
//! This one needs a major overhaul, it is autogenerating quality scores etc. 
//! Will wait to get other things set up first
use std::io::Write;
use bincode::de::read;
use log::{debug, error, warn, info};
use std::fs::File;
use std::path::PathBuf;
use std::collections::{HashMap, VecDeque};
use simple_rng::{NeatRng, NeatRngError};
use flate2::Compression;
use flate2::write::GzEncoder;
use thiserror::Error;
use indicatif::ProgressBar;

use crate::structs::mutated_map::{MutatedMap, MutatedMapError};
use crate::structs::fasta_map::{FastaMapError, SequenceBlock};
use crate::models::quality_scores::{self, QualityScoreModel};
use crate::structs::nucleotides::Nucleotide;
use crate::file_tools::file_io::{create_output_file, append_to_file};
use crate::models::sequencing_error_model::{SeqModelError, SequencingErrorModel, SequencingErrorType};
use crate::structs::nucleotides::Nucleotide::N;
use crate::structs::variants::{Genotypes, Variant};

#[derive(Error, Debug)]
pub enum FastqToolsError {
    #[error("Error writing bgzip fastq block {0}")]
    FastqWriteError(String),
    #[error("Quality scores missing for read 2")]
    MissingQScores,
    #[error("Error reading bgzip fastq block {0}")]
    FastqReadError(String),
    #[error("Mismatch between indexing and reads set for block {0}")]
    InvalidFastqBlock(String),
    #[error("Fastq Tools reported a FastaMap error: {0}")]
    FastaMapError(#[from] FastaMapError),
    #[error("Fastq tools reported a error model error: {0}")]
    ErrorModelError(#[from] SeqModelError),
    #[error("Fastq tools reported an IO error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Fastq tools reported an Rng Error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("Error from MutatedMap: {0}")]
    MutatedMapError(#[from] MutatedMapError),
    #[error("Error locating a matching mutated map")]
    FindMapError,
    #[error("Paired ended declared, but r2 buffer empty")]
    BufferInitError,
    #[error("Malformed read during pair-ended analysis: {0}")]
    MalformedReadError(String),
    #[error("Strand value must either be 1 or 2, received {0}")]
    StrandError(usize),
}

fn reverse_complement(sequence: Vec<Nucleotide>) -> Vec<Nucleotide> {
    // Returns the reverse complement of a vector of u8's representing a DNA sequence.
    let length = sequence.len();
    let mut rev_comp = Vec::new();
    for i in (0..length).rev() {
        rev_comp.push(sequence[i].complement())
    }
    rev_comp
}

fn find_map (interval: (usize, usize), mutated_maps: &Vec<MutatedMap>) -> Result<&MutatedMap, FastqToolsError> {
    for map in mutated_maps {
        if map.contains(interval)? {
            return Ok(map);
        }
    }
    Err(FastqToolsError::FindMapError)
}

pub fn write_block_fastq<T: Write, W: Write> (
    block_fragments: Vec<(usize, usize)>,
    block_map: &MutatedMap,
    paired_ended: bool,
    buffer1: &mut GzEncoder<T>,
    buffer2: &mut GzEncoder<W>,
    read_length: usize,
    read_name_prefix: &str,
    quality_score_model: &QualityScoreModel,
    sequencing_error_model: &SequencingErrorModel,
    rng: &mut NeatRng,
) -> Result<(), FastqToolsError> {
    // The idea here is to write a temporary fastq file that we can later recombine into
    // a single fastq file that will the output of NEAT. These should be the coordinates relative to the reference
    // block reads: We'll expect these to be sorted, because we think that will make this 
    //     part eaiser, but this could be an erroneous assumption.
    // block_maps: pairs the sequence block with the variants on that block
    // output_buffer: the buffer where we are writing data to, should be writing to a file
    // read_length: The length of the reads in the output data set
    // reverse_strand: Will determine whether these reads get a /1 or /2 at the end
    // quality_score_model: Generate quality scores for a read.
    // sequencing_error_model: Model to generate sequencing errors for a read.
    // rng: The NeatRng random number generator for the this run.
    let bar: ProgressBar = ProgressBar::new(block_fragments.len() as u64);
    let mut counter: usize = 1;
    let sequence_block: SequenceBlock = SequenceBlock::from(&block_map.sequence_block)?;
    // These reads should be in order of sequence block, but if we find one out of order, we can tuck it
    // in the back and maybe it matches later. If not it'll just get skipped.
    for (start, end) in block_fragments {
        bar.tick();
        let fragment = sequence_block.get_subseq(start, end)?;
        let mut read1_variants: HashMap<usize, &Variant> = HashMap::new();
        let mut reads1_flagged: Vec<usize> = Vec::new();
        let mut read2_variants: HashMap<usize, &Variant> = HashMap::new();
        let mut reads2_flagged: Vec<usize> = Vec::new();
        // So now we find the variants that are within the read window
        for (pos, variant) in &block_map.variant_map {
            if (start..(start+read_length)).contains(&pos) {
                // This gives is the index within the sequence vector
                let var_pos = pos - start;
                read1_variants.insert(var_pos, variant);
                reads1_flagged.push(var_pos);
            }
            // Not else if because we need to check both indpendently, they could overlap.
            if paired_ended == true && ((end-read_length)..end).contains(&pos) {
                // This gives the index within the reversed sequence vector
                // This calculates the distance from the last element. The math seemed to work.
                let var_pos = (end - 1) - pos;
                read2_variants.insert(var_pos, variant);
                reads2_flagged.push(var_pos);
            }
        }

        let quality_scores_1 = quality_score_model.generate_quality_scores(read_length, rng)?;
        apply_variants_and_write_sequence(
            &fragment,
            &reads1_flagged,
            &read1_variants,
            read_length,
            &read_name_prefix,
            counter,
            1,
            buffer1,
            quality_scores_1,
            sequencing_error_model,
            rng,
        )?;

        if paired_ended {
            // open second buffer
            let quality_scores_2 = quality_score_model.generate_quality_scores(read_length, rng)?;
            apply_variants_and_write_sequence(
                &reverse_complement(fragment),
                &reads2_flagged,
                &read2_variants,
                read_length,
                &read_name_prefix,
                counter,
                2,
                buffer2,
                quality_scores_2,
                sequencing_error_model,
                rng,
            )?;
        }
        counter += 1;
    }
    Ok(())
}

pub fn combine_temp_fastqs(files: Vec<PathBuf>, shuffle: bool) -> Result<(), FastqToolsError> {
    // This will take all the temporary fastqs we wrote out on a per-block level and 
    // combine them into one fastq. An improvement here might be to randomize the order:
    // We could do this with a hashmap from the original name to a number in a shuffled vector
    // of indexes. We'll make that optional.
    if shuffle {
        warn!("Fastq shuffle implementation comming")
    }
    for file in files {
        println!("{:?}", file)
    }
    Ok(())
}

fn apply_variants_and_write_sequence<T: Write> (
    sequence: &Vec<Nucleotide>,
    flagged_positions: &Vec<usize>,
    variant_map: &HashMap<usize, &Variant>,
    read_length: usize,
    read_name_prefix: &str,
    read_num: usize,
    read_strand: usize,
    buffer: &mut GzEncoder<T>,
    quality_scores: Vec<usize>,
    sequencing_error_model: &SequencingErrorModel,
    rng: &mut NeatRng,
) -> Result<(), FastqToolsError> {
    // This process will apply variants and sequencing errors to a fragment, then write the reads,
    // as applicable.
    // we start by writing the read name:
    // Write read name
    buffer.write(
        format!("@{}{}:{}\n", 
        &read_name_prefix, 
        &read_num, 
        read_strand
    ).as_bytes())?;
    // Use a while loop so we can skip for deletions
    // To figure out where we are on the contig, we use the reference_read_position. This helps us find relevant variants
    let mut reference_position = 0;
    let mut bases_written = 0;
    let mut seq = String::new();
    let mut quality_del_offset = 0;
    while (reference_position < sequence.len()) && (bases_written < read_length) {
        let reference_base = sequence[reference_position];
        let mut base_to_write: Vec<Nucleotide> = vec![reference_base];
        // Figure out what to actually write at this position
        if reference_base == N {
            // skip this block
        } else if flagged_positions.contains(&reference_position) {
            // Potentially contains a mutation, so we will write that, if applicable
            let variant = variant_map[&reference_position];
            match variant.genotype {
                Genotypes::Heterozygous => {
                    // sequence splicing isn't going to work because if I change the sequence, then it throws off my positions.
                    // 50/50 chance it gets applied
                    if rng.random()? < 0.5 {
                        base_to_write = variant.alternate.clone();
                    }
                },
                Genotypes::Homozygous => {
                    // Always applied
                    base_to_write = variant.alternate.clone()
                },
            }
        } else {
            // if no variant, test for an error
            let score = quality_scores[reference_position - quality_del_offset];
            let prob = sequencing_error_model.convert_score(score)?;
            if rng.random()? < prob {
                // This position contains a sequencing error
                debug!("Creating sequencing error");
                let error = sequencing_error_model
                    .generate_sequencing_error(reference_base, rng)?;
                match error {
                    SequencingErrorType::SnpError(base) => {
                        debug!("Snp error");
                        // write out new base, keep quality scores unchanged 
                        // We're converting this from our simplified representation to standard utf-8 character encoding to write to file
                        base_to_write = vec![base];
                    },
                    SequencingErrorType::DeletionError(length) => {
                        debug!("Deletion error");
                        // skip `length` bases, leave quality scores untouched (because it's close enough)
                        // First we check to make sure we haven't deleted too many bases this read
                        if reference_position + length >= sequence.len() {
                            // In this case, we would be deleting too many bases, so we'll skip
                        } else {
                            // for a deletion, we still writet the reference base, then we skip the next length bases
                            // base_to_write remains unchanged
                            quality_del_offset += length;
                            reference_position += length;
                        }
                    },
                    SequencingErrorType::InsertionError(vec) => {
                        debug!("Insertion error");
                        // Write out current base, then insertion vec with low quality scores, but don't increment index,
                        // so the next base is the one originally following the current base
                        let mut insert_len = vec.len();
                        let mut insertion: Vec<Nucleotide> = Vec::new();
                        insertion.push(reference_base);
                        // We have to be careful here not to write out more bases than a read length, for insertions
                        // at the ends of the reads. The +1 signifies the reference_base already included
                        if bases_written + insert_len >= read_length {
                            // If the insertion takes us past the read_length, we'll trim the insertion
                            insert_len = bases_written + insert_len - read_length;
                        };
                        for i in 0..insert_len {
                            insertion.push(vec[i])
                        }
                    }
                }
            }
        }
        for base in base_to_write {
            seq.push(base.into());
            bases_written += 1;
        }
        reference_position += 1;
    } 
    buffer.write(&seq.into_bytes())?;
    // Read has been written at this point
    
    // Write out the plus sign spacer
    buffer.write(b"\n+\n")?;
    // Convert quality scores to a character vec, write out `read_length` scores and a line break
    buffer.write(&quality_scores_to_char_vec(quality_scores)?)?;
    // Add a line break
    buffer.write(b"\n")?;
    Ok(())
}

pub fn quality_scores_to_char_vec(array: Vec<usize>) -> Result<Vec<u8>, FastqToolsError> {
    let mut score_vec = Vec::new();
    for score in array {
        score_vec.push(
            (score + 33) as u8
        )
    }
    Ok(score_vec)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::variants::{Variant, VariantType};
    use crate::structs::nucleotides::Nucleotide::*;
    use std::fs;

    #[test]
    fn test_find_map() {
        let sequence_block = PathBuf::from("chr1_001000_002000.json");
        let variant = Variant::new(
            VariantType::Deletion,
            1002,
            &vec![G, G, C],
            &vec![G],
            &mut vec![1,0]
        ).unwrap();
        let variant_vec = vec![variant];
        let mutated_map = MutatedMap::new(sequence_block.clone(), variant_vec).unwrap();
        assert_eq!(find_map((1050, 1550), &vec![mutated_map]).unwrap().sequence_block, sequence_block);
        
    }

    #[test]
    fn test_qual_score_to_write() {
        let qual_scores = vec![33, 25, 37, 28, 15, 33, 33, 37, 37, 25];
        let qual_string = "B:F=0BBFF:".as_bytes();
        assert_eq!(qual_string, quality_scores_to_char_vec(qual_scores).unwrap())
    }

    #[test]
    fn test_reverse_complement() {
        let read: Vec<Nucleotide> = vec![A, A, A, A, C, C, C, C, C];
        let revcomp: Vec<Nucleotide> = vec![G, G, G, G, G, T, T, T, T];
        assert_eq!(reverse_complement(read), revcomp);
    }

    #[test]
    fn test_apply_variants() {
        let sequence = vec![A, C, G, T, T, A, T, G];
        let variant1 = Variant::new(
                VariantType::SNP,
                1,
                &vec![T],
                &vec![C],
                &mut vec![1, 0],
        ).unwrap();
        let variant2 = Variant::new(
                VariantType::Deletion,
                3,
                &vec![T, T],
                &vec![T],
                &mut vec![0,1],
        ).unwrap();
        let variant_map = HashMap::from([
            (1, &variant1),
            (3, &variant2),
        ]);
        let flagged_positions = vec![1,3];
        let read_name_prefix = "neat_generated_";
        let file = PathBuf::from("fake.fq");
        let outfile = create_output_file(&file, true).unwrap();
        let mut buffer = GzEncoder::new(outfile, Compression::default());
        let qual_scores = vec![33, 25, 37, 28, 15, 33, 33, 37];
        let sequencing_error_model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let result = apply_variants_and_write_sequence(
            &sequence, 
            &flagged_positions, 
            &variant_map, 
            8, 
            read_name_prefix, 
            1, 
            1,
            &mut buffer, 
            qual_scores, 
            &sequencing_error_model, 
            &mut rng
        );

        match result {
            Ok(()) => assert!(true),
            Err(_) => assert!(false),
        }
        fs::remove_file(file).unwrap();
    }
}
