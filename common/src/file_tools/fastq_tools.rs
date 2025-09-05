//! This library writes either single ended or paired-ended fastq files.
//! Need to update this method. We want to use the data structures and we want to make sure
//! this function is generic enough to work with the fragmented method we are implementing
//! This one needs a major overhaul, it is autogenerating quality scores etc. 
//! Will wait to get other things set up first
use std::io::Write;
use std::fs;
use std::fs::File;
use std::path::PathBuf;
use std::iter::Extend;
use simple_rng::{NeatRng, NeatRngError};
use flate2::Compression;
use flate2::write::ZlibEncoder;
use tempfile::TempDir;
use thiserror::Error;

use crate::structs::fasta_map::{FastaMapError, SequenceBlock};
use crate::models::quality_scores::QualityScoreModel;
use crate::structs::nucleotides::Nucleotide;
use crate::file_tools::file_io::create_output_file;
use crate::models::sequencing_error_model::{SeqModelError, SequencingErrorModel, SequencingErrorType};
use crate::structs::variants::{Genotypes, Variant, VariantType};

#[derive(Error, Debug)]
pub enum BlockFastQError {
    #[error("Error writing bgzip fastq block {0}")]
    FastqWriteError(String),
    #[error("Error reading bgzip fastq block {0}")]
    FastqReadError(String),
    #[error("Mismatch between indexing and reads set for block {0}")]
    InvalidFastqBlock(String),
    #[error("Fastq Tools reported a FastaMap error: {0}")]
    FastaMapError(FastaMapError),
    #[error("Fastq tools reported a error model error: {0}")]
    ErrorModelError(SeqModelError),
    #[error("Fastq tools reported an IO error: {0}")]
    IoError(std::io::Error),
    #[error("Fastq tools reported an Rng Error: {0}")]
    RngError(NeatRngError),
}

impl From<NeatRngError> for BlockFastQError {
    fn from(error: NeatRngError) -> Self {
        BlockFastQError::RngError(error)
    }
}

impl From<std::io::Error> for BlockFastQError {
    fn from(error: std::io::Error) -> Self {
        BlockFastQError::IoError(error)
    }
}

impl From<FastaMapError> for BlockFastQError {
    fn from(error: FastaMapError) -> Self {
        BlockFastQError::FastaMapError(error)
    }
}

impl From<SeqModelError> for BlockFastQError {
    fn from(error: SeqModelError) -> Self {
        BlockFastQError::ErrorModelError(error)
    }
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
pub fn write_block_fastq_bgz(
    block_reads: &Vec<(usize, usize, usize, usize)>,
    read_number_vec: Vec<usize>,
    block_filename: &PathBuf,
    variants: Vec<&Variant>, 
    read_length: usize,
    paired_ended: bool,
    temp_dir: TempDir,
    quality_score_model: &QualityScoreModel,
    sequencing_error_model: &SequencingErrorModel,
    rng: &mut NeatRng,
) -> Result<(), BlockFastQError> {
    // Takes:
    // block reads: A list of tuples with coordinates of the reads for the block.
    // block_filename: A file name holding the sequence for this block. This should be an already mutated sequence.
    // variants: this should be a vetted list of variants to apply to this block.
    //      We assume here that all variants submitted will be processed.
    // read_length: The length of the reads. Any output reads will be trimmed to this length.
    // paired_ended: boolean to set paired ended mode on or off.
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
    //
    // Extract the filename from the block file name, this will hopefully give us 
    // unique filenames, if we've set things up right.
    let block_name = &block_filename
        .file_stem()
        .unwrap()
        .display()
        .to_string();
    // We don't want these files to have a valid extension, so they aren't accidentally 
    // read, but we need to distinguish them, so they are bgfq
    let filename1 = Some(temp_dir
        .path()
        .join(format!("{}_r1.bgfq", block_name)))
        .unwrap();
    let filename2 = Some(temp_dir
        .path()
        .join(format!("{}_r2.bgfq", block_name)))
        .unwrap();

    let read_name_prefix = "neat_generated_".to_string();
    
    // open the file and append lines
    let outfile1 = create_output_file(&filename1, false)
        .expect(&format!("Error creating output {:?}. If file alread exists, remove.", filename1));
    // open the second file and append lines
    let outfile2 = create_output_file(&filename2, false)
        .expect(&format!("Error opening output {:?}. If file alread exists, remove.", filename2));

    // length of read_number_vec and length of block_reads should be the same
    if read_number_vec.len() != block_reads.len() {
        return Err(BlockFastQError::InvalidFastqBlock(block_name.to_string()))
    }

    // load block into memory. Should only be 128kb long or so of sequence plus some strings. Easy peasy.
    let sequence_block = SequenceBlock::from(block_filename)?;
    // setting the buffer to 2*BUFSIZE which is approximately 0.5MB. idk, man.
    let mut buffer_r1  = ZlibEncoder::new(
        outfile1, Compression::default()
    );
    let mut buffer_r2 = ZlibEncoder::new(
        outfile2, Compression::default()
    );
    // Get the maximum number of deletions: We assume that all deletions
    // are included, the read is paired ended, and reads overlap.
    let mut max_num_deletions = 0;
    for variant in &variants {
        match variant.variant_type {
            VariantType::Deletion => max_num_deletions += variant.reference.len() - 1,
            _ => continue,
        }
    }
    // write out sequence by sequence, throwing in random errors and adding variants as appropriate
    // We combine the read1 and read2 reads so that we get a match on the naming. This could be done 
    // with a second loop as well, which would allow us to atomize the write sequences process, but 
    // the hope is that this is slightly more efficient. If this part is super slow, that may be
    // something to try
    for (index, read_num) in read_number_vec.iter().enumerate() {
        let (position1, position2, position3, position4) = block_reads[index];
        // Load actual sequence into a short lived variable
        // This is tricky, if the padding goes past the end of the read,
        // we'll get an index error.
        let write_padding = {
            if position4 - position1 < max_num_deletions {
                position4 - position1
            } else {
                max_num_deletions
            }
        };
        let mut current_read_1 = sequence_block
            .get_subseq(position1, position2 + write_padding)?;
        // Apply variants to read
        apply_variants_to_sequence(
            &mut current_read_1, 
            &variants, 
            position1, 
            position2, 
            rng,
            read_length, 
            write_padding, 
            &read_name_prefix, 
            *read_num, 
            &mut buffer_r1, 
            quality_score_model, 
            sequencing_error_model
        ).unwrap();
        
        if paired_ended {
            // Read 2 is read from the 3' end so we will write out the reverse complement
            let mut current_read_2 = reverse_complement(
                sequence_block
                    .get_subseq(position3 - write_padding, position4)?
            );

            // Apply variants to read 2
            apply_variants_to_sequence(
                &mut current_read_2, 
                &variants, 
                position3, 
                position4, 
                rng,
                read_length, 
                write_padding, 
                &read_name_prefix, 
                *read_num, 
                &mut buffer_r2, 
                quality_score_model, 
                sequencing_error_model
            ).unwrap();
        }
    }

    if !paired_ended {
        // We don't need this second file in this case
        fs::remove_file(filename2)?;
    }

    Ok(())
}

fn quality_scores_to_char_vec(array: Vec<usize>) -> Result<Vec<u8>, BlockFastQError> {
    let mut score_vec = Vec::new();
    for score in array {
        score_vec.push(
            (score + 33) as u8
        )
    }
    Ok(score_vec)
}

fn apply_variants_to_sequence (
    sequence: &mut Vec<Nucleotide>,
    variants: &Vec<&Variant>,
    position1: usize,
    position2: usize,
    rng: &mut NeatRng,
    read_length: usize,
    write_padding: usize,
    read_name_prefix: &str,
    read_num: usize,
    buffer: &mut ZlibEncoder<File>,
    quality_score_model: &QualityScoreModel,
    sequencing_error_model: &SequencingErrorModel,
) -> Result<(), BlockFastQError> {
    // This process will be needed for each read
    for variant in variants {
        if (position1..position2).contains(&variant.location) {
            // variant is in scope
            match variant.genotype {
                Genotypes::Heterozygous => {
                    // 50/50 chance it gets applied
                    if rng.random()? < 0.5 {
                        sequence_splicer(
                            sequence, 
                            variant,
                        )?;
                    }
                },
                Genotypes::Homozygous => {
                    // Always applied
                    sequence_splicer(
                        sequence, 
                        variant,
                    )?;
                },
            }
        }
    }
    // We build in some padding for deletions
    let qual_score_len = read_length + write_padding;
    // we may need to update quality scores
    let mut quality_scores_1 = quality_score_model
        .generate_quality_scores(qual_score_len, rng)?;

    let mut del_len_max_r1 = write_padding;

    // Write read name
    buffer.write(
        format!("@{} 1:{}\n", &read_name_prefix, &read_num
    ).as_bytes())?;

    let mut index = 0;
    let mut written = 0;

    // This section writes out the actual sequence. We stop when we reach the end of the read
    // or we write out 1 read_length worth of characters. Hopefully we catch all the edge cases.
    'outer: while (index < sequence.len()) && (written < read_length) {
        // This section will work like this:
        //     - for each index, check for an error with the quality score
        //     - if an error, handle according to type, including incrementing the index
        //     - Tricky part is, we only want to write out READ LENGTH number of bases,
        //           but there are indels of random length, so we have to be flexible
        //     - What if too many deletions accumulate in a read?
        let score = quality_scores_1[index];
        let prob = sequencing_error_model.convert_score(score as usize)?;
        
        // grab current base
        let current_base = sequence[index];
        
        // Lower score will translate to a higher probability of error
        if rng.random()? < prob {
            // sequencing error
            let error = sequencing_error_model
                .generate_sequencing_error(current_base, rng)?;
            match error {
                SequencingErrorType::SnpError(base) => {
                    // write out new base, keep quality scores unchanged 
                    // We're converting this from our simplified representation to standard utf-8 character encoding to write to file
                    buffer.write(&[base as u8])?;
                    written += 1;
                },
                SequencingErrorType::DeletionError(length) => {
                    // skip `length` bases, leave quality scores untouched (because it's close enough)
                    // First we check to make sure we haven't deleted too many bases this read
                    if length > del_len_max_r1 {
                        index += length;
                        del_len_max_r1 -= length;
                    } else {
                        // We added enough deletion errors to this read, skip
                    }
                },
                SequencingErrorType::InsertionError(vec) => {
                    // Write out current base, then insertion vec with low quality scores, but don't increment index,
                    // so the next base is the one originally following the current base
                    let insert_len = &vec.len();
                    let mut insertion: Vec<u8> = Vec::new();
                    insertion.push(current_base as u8);
                    written += 1;
                    // We have to be careful here not to write out more bases than a read length, for insertions
                    // at the ends of the reads
                    if written >= read_length { break 'outer }
                    for base in vec {
                        insertion.push(base as u8);
                        written += 1;
                        // We have to be careful here not to write out more bases than a read length, 
                        // for insertions at the ends of the reads
                        if written >= read_length { break 'outer }
                    }
                    buffer.write(&insertion)?;

                    let low_score = vec![5; *insert_len];
                    quality_scores_1.splice(
                        (index+1)..(index+1), low_score.iter().cloned()
                    );
                }
            }
        } else {
            buffer.write(&[current_base as u8])?;
            written += 1;
        }
        index = index + 1;
    }
    
    // Write out the plus sign spacer
    buffer.write(b"\n+\n")?;
    // Convert quality scores to a character vec, write out `read_length` scores and a line break
    let mut char_vec = quality_scores_to_char_vec(quality_scores_1)?;
    // Trim toinal length. If this throws an error, then we have a bug.
    char_vec = char_vec[..read_length].to_vec();
    // Add a line break
    char_vec.extend(b"\n");
    // and write
    buffer.write(&char_vec)?;
    Ok(())
}

fn sequence_splicer (
    sequence: &mut Vec<Nucleotide>,
    variant: &Variant,
) -> Result<(), BlockFastQError> {
    let reference_length = variant.reference.len();
    sequence.splice(
        variant.location..(variant.location + reference_length),
        variant.alternate.iter().cloned(),
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::variants::{Variant, VariantType};
    use crate::structs::nucleotides::Nucleotide::*;

    #[test]
    fn test_splicing_deletion() {
        let mut sequence = vec![C, G, G, G, C, A, A, G, T, C];
        let variant = Variant::new(
            VariantType::Deletion,
            2,
            &vec![G, G, C],
            &vec![G],
            &mut vec![1,0]
        ).unwrap();
        sequence_splicer(&mut sequence, &variant).unwrap();
        assert_eq!(sequence, vec![C, G, G, A, A, G, T, C,])
    }

    #[test]
    fn test_splicing_insertion() {
        let mut sequence = vec![C, G, G, G, C, A, A, G, T, C];
        let variant = Variant::new(
            VariantType::Insertion,
            2,
            &vec![G],
            &vec![G, T, T],
            &mut vec![1,1]
        ).unwrap();
        sequence_splicer(&mut sequence, &variant).unwrap();
        assert_eq!(sequence, vec![C, G, G, T, T, G, C, A, A, G, T, C,])
    }

    #[test]
    fn test_variables() {
        let v = 5;
        let mut x = v;
        for _ in 0..2 {
            x -= 1
        }
        assert_eq!(x, 3);
        assert_eq!(v, 5);
    }

    #[test]
    fn test_reverse_complement() {
        let read: Vec<Nucleotide> = vec![A, A, A, A, C, C, C, C, C];
        let revcomp: Vec<Nucleotide> = vec![G, G, G, G, G, T, T, T, T];
        assert_eq!(reverse_complement(read), revcomp);
    }
}
