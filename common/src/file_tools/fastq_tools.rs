//! This library writes either single ended or paired-ended fastq files.
//! Need to update this method. We want to use the data structures and we want to make sure
//! this function is generic enough to work with the fragmented method we are implementing
//! This one needs a major overhaul, it is autogenerating quality scores etc. 
//! Will wait to get other things set up first
use std::io::Write;
use log::{debug, error, info};
use std::fs::File;
use std::path::PathBuf;
use std::collections::HashMap;
use simple_rng::{NeatRng, NeatRngError};
use flate2::Compression;
use flate2::write::GzEncoder;
use thiserror::Error;
use indicatif::ProgressBar;

use crate::structs::mutated_map::{MutatedMap, MutatedMapError};
use crate::structs::fasta_map::{FastaMapError, SequenceBlock};
use crate::models::quality_scores::QualityScoreModel;
use crate::structs::nucleotides::Nucleotide;
use crate::file_tools::file_io::create_output_file;
use crate::models::sequencing_error_model::{SeqModelError, SequencingErrorModel, SequencingErrorType};
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

/// write_fastq. This is freeform writing. I think each thread will get a sequence block
/// file is what I'm thinking, and it will write out a small fastq file of it's reads. 
/// After that I don't know. We might need to shuffle the reads? Not sure how that will
/// work exactly. Reading in fastqs and shuffling them would not be easy. That's the part
/// that's tripping me up. Also each file will have it's own order and numbering of reads
/// because I won't be able to synchronize that. Well, maybe I can. There is message passing
/// between threads. Should I write out temporary fastq files? Reading in all those
/// fastq fragments is going to be a huge amount of time. Maybe I want the temp fastqs
/// to be more machine readable than human.
/// 
/// Thought what if I use blockgzip via the gzp crate? I'm basically writing blocks.
/// Okay, so we create a buffer and we write the data to the buffer then I guess compress it?
/// I could set it to best compression, set to "best" Maybe we need to figure out an appropriate
/// buffer size first.
pub fn write_fastq(
    all_reads: &mut Vec<(String, usize, usize, Option<usize>, usize)>,
    mutated_maps: &HashMap<String, Vec<MutatedMap>>,
    read_length: usize,
    paired_ended: bool,
    output_files: (PathBuf, &mut Option<PathBuf>),
    quality_score_model: &QualityScoreModel,
    sequencing_error_model: &SequencingErrorModel,
    overwrite_output: bool,
    rng: &mut NeatRng,
) -> Result<(), FastqToolsError> {
    // Takes:
    // all_reads: All reads along with contig name.
    // mutated_maps: An vector of MutatedMap objects (one per block), indexed by contig name.
    //      This holds sequence and applicable mutations
    // contig_names: The keys to the other two maps.
    // read_length: The length of the reads. Any output reads will be trimmed to this length.
    // paired_ended: boolean to set paired ended mode on or off.
    // output_files: A tuple with the PathBuf objects to read1 and (if applicable) read2.
    // quality_score_model: The statistical model we will use to generate quality scores.
    //      We only need a reference to it, we don't want to copy this everywhere as these
    //      may get very large in practice
    // sequencing_error_model: The statistical model we'll need to generate read errors.
    //      Should work simply, will need to handle indels carefully. We don't want to 
    //      copy this one te every thread either.
    // rng: the random number generator for the run
    //
    // returns:
    // Error if there is a problem or else it successfully wrote the file(s)

    // TODO future improvement idea. Add in some randomly generated reads to simulate telomeres.
    let read_name_prefix = "neat_generated_".to_string();
    // shuffle reads
    rng.shuffle_in_place(all_reads)?;
    info!("Writing output fastq file(s)");
    let bar = ProgressBar::new(all_reads.len() as u64);
    // We'll start numbering reads at 1
    let mut read_index = 1;
    // Open our files for writing
    // Since Rust won't be able to write if there is a possibility the open file will not be there,
    // we'll have to use a dummy file for read 2, which will never get written to, then will be 
    // deleted at the end
    debug!("Output filename1: {:?}", &output_files.0);
    let outfile1 = create_output_file(&output_files.0, overwrite_output)?;
    debug!("Successfully created output files");
    let mut buffer_r1: GzEncoder<File>  = GzEncoder::new(
        outfile1, Compression::default()
    );
    debug!("We made it THIS far!");
    let mut files_written = Vec::new();
    for read in all_reads {
        bar.inc(1);
        let contig_name = &read.0;
        let contig_blocks = &mutated_maps[contig_name];
        let current_map = find_map((read.1, read.4), contig_blocks)?;
        let block = SequenceBlock::from(&current_map.sequence_block)?;
        let fragment: Vec<Nucleotide> = block.get_subseq(read.1, read.4)?;
        let quality_score_r1 = quality_score_model.generate_quality_scores(read_length, rng)?;
        let mut quality_score_r2: Option<Vec<usize>> = None;
        match paired_ended {
            true => quality_score_r2 = Some(quality_score_model.generate_quality_scores(read_length, rng)?),
            false => (),
        }
        
        // Find variants in the read and translate the coordinates to read coordinates
        let mut read_1_variants: HashMap<usize, &Variant> = HashMap::new();
        let mut read_1_flagged: Vec<usize> = Vec::new();
        for (pos, variant) in &current_map.variant_map {
            if (read.1..read.2).contains(&pos) {
                let start_pos = pos - read.1;
                read_1_variants.insert(start_pos, variant);
                read_1_flagged.push(start_pos)
            }
        }

        // This is for read1
        let result_r1 = apply_variants_and_write_sequence(
            fragment.clone(), 
            &read_1_flagged,
            &read_1_variants, 
            read_length, 
            &read_name_prefix, 
            read_index, 
            1,
            &mut buffer_r1, 
            quality_score_r1,
            sequencing_error_model,
            rng,
        );

        match result_r1 {
            Err(error) => return Err(error),
            _ => {
                read_index += 1;
                files_written.push(output_files.0.clone());
            },
        }

        if paired_ended {
            match output_files.1 {
                Some(file) => {
                    debug!("Output filename2: {:?}", &file);
                    let outfile2 = create_output_file(&file, overwrite_output)?;
                    let mut buffer_r2 = GzEncoder::new(outfile2, Compression::default());
                    // This is for read2, so we need to reverse the coordinates and do a reverse complement of the sequence
                    let pos1 = match read.3 {
                        Some(num) => num,
                        None => return Err(FastqToolsError::MalformedReadError("read.3".to_string()))
                    };
                    let pos2 = read.4;
                    let quality_scores = match quality_score_r2 {
                        Some(vec) => vec,
                        None => return Err(FastqToolsError::MissingQScores)
                    };

                    let r2_frag = reverse_complement(fragment);

                    // Find variants in the read and translate the coordinates to read coordinates
                    let mut read_2_variants: HashMap<usize, &Variant> = HashMap::new();
                    let reverse_index: Vec<usize> = (0..r2_frag.len()).rev().collect();

                    let mut r2_flagged_positions = Vec::new();
                    for (pos, variant) in &current_map.variant_map {
                        if (pos1..pos2).contains(&pos) {
                            let index = r2_frag.len() - pos1;
                            // The index is relative to the original fragment, now we need the reversed coordinate
                            let start_pos = reverse_index[index];
                            read_2_variants.insert(start_pos, variant);
                            r2_flagged_positions.push(start_pos);
                        }
                    }

                    let result_r2 = apply_variants_and_write_sequence(
                        r2_frag, 
                        &r2_flagged_positions,
                        &read_2_variants, 
                        read_length, 
                        &read_name_prefix, 
                        read_index, 
                        2,
                        &mut buffer_r2, 
                        quality_scores,
                        sequencing_error_model,
                        rng,
                    );

                    match result_r2 {
                        Err(error) => return Err(error),
                        _ => {
                            read_index += 1;
                            files_written.push(file.to_path_buf());
                        },
                    }
                },
                None => {},
            }
        }
    }
    bar.finish();
    info!("Fastq(s) successfully written: {:?}!", files_written);
    Ok(())
}

fn quality_scores_to_char_vec(array: Vec<usize>) -> Result<Vec<u8>, FastqToolsError> {
    let mut score_vec = Vec::new();
    for score in array {
        score_vec.push(
            (score + 33) as u8
        )
    }
    Ok(score_vec)
}

fn apply_variants_and_write_sequence (
    sequence: Vec<Nucleotide>,
    flagged_positions: &Vec<usize>,
    variant_map: &HashMap<usize, &Variant>,
    read_length: usize,
    read_name_prefix: &str,
    read_num: usize,
    read_strand: usize,
    buffer: &mut GzEncoder<File>,
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
    while (reference_position < sequence.len()) && (bases_written < read_length) {
        let reference_base = sequence[reference_position];
        let mut base_to_write: Vec<Nucleotide> = vec![reference_base];
        // Figure out what to actually write at this position
        if flagged_positions.contains(&reference_position) {
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
            let score = quality_scores[reference_position];
            let prob = sequencing_error_model.convert_score(score as usize)?;
            if rng.random()? < prob {
                // This position contains a sequencing error
                debug!("Creating sequencing error");
                let error = sequencing_error_model
                    .generate_sequencing_error(reference_base, rng)?;
                match error {
                    SequencingErrorType::SnpError(base) => {
                        // write out new base, keep quality scores unchanged 
                        // We're converting this from our simplified representation to standard utf-8 character encoding to write to file
                        base_to_write = vec![base];
                    },
                    SequencingErrorType::DeletionError(length) => {
                        // skip `length` bases, leave quality scores untouched (because it's close enough)
                        // First we check to make sure we haven't deleted too many bases this read
                        if reference_position + length >= sequence.len() {
                            // In this case, we would be deleting too many bases, so we'll skip
                        } else {
                            // for a deletion, we still writet the reference base, then we skip the next length bases
                            // base_to_write remains unchanged
                            reference_position += length - 1;
                        }
                    },
                    SequencingErrorType::InsertionError(vec) => {
                        // Write out current base, then insertion vec with low quality scores, but don't increment index,
                        // so the next base is the one originally following the current base
                        let mut insert_len = vec.len();
                        let mut insertion: Vec<Nucleotide> = Vec::new();
                        insertion.push(reference_base);
                        // We have to be careful here not to write out more bases than a read length, for insertions
                        // at the ends of the reads. The +1 signifies the reference_base already included
                        if bases_written + insert_len + 1 >= read_length {
                            // If the insertion takes us past the read_length, we'll trim the insertion
                            insert_len = bases_written + insert_len + 1 - read_length;
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
            sequence, 
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
