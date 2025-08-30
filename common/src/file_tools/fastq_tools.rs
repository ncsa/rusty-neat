//! This library writes either single ended or paired-ended fastq files.
//! Need to update this method. We want to use the data structures and we want to make sure
//! this function is generic enough to work with the fragmented method we are implementing
//! This one needs a major overhaul, it is autogenerating quality scores etc. 
//! Will wait to get other things set up first
use std::io::Write;
use std::fs;
use std::path::Path;
use simple_rng::NeatRng;
use gzp::{bgzf, Compression};
use tempfile::TempDir;
use thiserror::Error;

use crate::file_tools::fasta_tools::BUFSIZE;
use crate::structs::fasta_map::{FastaMapError, SequenceBlock};
use crate::structs::nucleotides::sequence_array_to_string;
use crate::models::quality_scores::QualityScoreModel;
use crate::structs::nucleotides::complement;
use crate::file_tools::file_io::create_output_file;
use crate::models::sequencing_error_model::{SeqModelError, SequencingErrorModel};
use crate::structs::variants::Variant;

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
    IoError(std::io::Error)
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

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    // Returns the reverse complement of a vector of u8's representing a DNA sequence.
    let length = sequence.len();
    let mut rev_comp = Vec::new();
    for i in (0..length).rev() {
        rev_comp.push(complement(sequence[i]))
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
    block_variants: &Vec<Variant>,
    block_reads: &Vec<(usize, usize, usize, usize)>,
    read_number_vec: Vec<usize>,
    block_filename: &str,
    read_length: usize,
    paired_ended: bool,
    temp_dir: TempDir,
    quality_score_model: &QualityScoreModel,
    sequencing_error_model: &SequencingErrorModel,
    rng: &mut NeatRng,
) -> Result<(), BlockFastQError> {
    // Takes:
    // block_variants: A list of variants that are present on this sequence block.
    // block reads: A list of tuples with coordinates of the reads for the block.
    // block_filename: A file name holding the sequence for this block.
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
    let block_name = Path::new(block_filename)
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
    let mut outfile1 = create_output_file(&filename1, false)
        .expect(&format!("Error creating output {:?}. If file alread exists, remove.", filename1));
    // open the second file and append lines
    let mut outfile2 = create_output_file(&filename2, false)
        .expect(&format!("Error opening output {:?}. If file alread exists, remove.", filename2));

    // length of read_number_vec and length of block_reads should be the same
    if read_number_vec.len() != block_reads.len() {
        return Err(BlockFastQError::InvalidFastqBlock(block_name))
    }

    // load block into memory. Should only be 128kb long or so of sequence plus some strings. Easy peasy.
    let sequence_block = SequenceBlock::from(block_filename)?;
    // setting the buffer to 4*BUFSIZE which is approximately 1MB. idk, man.
    let buffer_r1  = bgzf::BgzfSyncWriter::with_capacity(
        outfile1, Compression::best(), 2*BUFSIZE
    );
    let buffer_r2 = bgzf::BgzfSyncWriter::with_capacity(
        outfile2, Compression::best(), 2*BUFSIZE
    );

    // write out sequence by sequence, throwing in random errors and adding variants as appropriate
    for (index, read_num) in read_number_vec.iter().enumerate() {
        let (position1, position2, position3, position4) = block_reads[index];
        let quality_scores = quality_score_model
            .generate_quality_scores(read_length as u32, &mut rng)?;

        // Load actual sequence into a short lived variable
        let current_read_1 = sequence_block
            .get_subseq(position1, position2)
            .unwrap();

        // Write read name
        buffer_r1.write(format!("@{} 1:{}", &read_name_prefix, read_num).as_bytes());
        let mut index = 0;
        let mut written = 0;
        while index < current_read_1.len() && written < read_length {
            // This section will work like this:
            //     - for each index, check for an error with the quality score
            //     - if an error, handle according to type, including incrementing the index
            //     - Tricky part is, we only want to write out READ LENGTH number of bases,
            //           but there are indels of random length, so we have to be flexible
            //     - What if too many deletions accumulate in a read?
            let score = quality_scores[index];
            let prob = sequencing_error_model.convert_score(score as usize);
            if rng.random()? < prob {

            }
            index = index + 1;
        }
        if paired_ended {
            buffer_r2.write(format!("@{} 2:{}", &read_name_prefix, read_num).as_bytes());
            for (base_index, base) in current_read_1.iter().enumerate() {
                todo!();
            }
        }
    }

    // for (contig, position1, position2) in reads_dataset.iter() {
       
    
        
    //     // Array as a string
    //     writeln!(
    //         &mut outfile1,
    //         // First read is the first read_length amount of the fragment
    //         "{}", sequence_array_to_string(&fragment.get(..read_length).unwrap())
    //     )?;
    //     // The stupid plus sign
    //     writeln!(&mut outfile1, "+")?;
    //     // Qual score of all F's for the whole thing.
    //     writeln!(&mut outfile1, "{}", quality_scores_to_str(quality_scores))?;
    //     index += 1;
    //     if paired_ended {
    //         // Need a quality score for this read as well
    //         let quality_scores = quality_score_model.generate_quality_scores(
    //             read_length, &mut rng
    //         );
    //         // The second read is the reverse complement at the end of the fragment
    //         let second_read = &fragment.get((fragment.len() - read_length)..).unwrap();
    //         // sequence name
    //         writeln!(
    //             &mut outfile2,
    //             "@{}/2",
    //             read_name_prefix.clone() + &index.to_string()
    //         )?;
    //         // Array as a string
    //         writeln!(
    //             &mut outfile2,
    //             "{}",
    //             sequence_array_to_string(&reverse_complement(second_read))
    //         )?;
    //         // The stupid plus sign
    //         writeln!(&mut outfile2, "+")?;
    //         // Qual score of all F's for the whole thing.
    //         writeln!(&mut outfile2, "{}", quality_scores_to_str(quality_scores))?;
    //     }
    // }
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

    #[test]
    fn test_reverse_complement() {
        let read: Vec<u8> = vec![0, 0, 0, 0, 1, 1, 1, 1];
        let revcomp: Vec<u8> = vec![2, 2, 2, 2, 3, 3, 3, 3];
        assert_eq!(reverse_complement(&read), revcomp);
    }

    // #[test]
    // fn test_write_fastq() {
    //     let mutated_reference_seq = vec![0, 2, 3, 0, 1, 3, 1, 0, 2, 3, 2, 3, 3, 1, 1, 3];
    //     let mutated_fasta_map = Box::new(HashMap::from([
    //         ("chr1".to_string(), mutated_reference_seq)
    //     ]));
    //     let mut reads_dataset = Box::new(vec![
    //         ("chr1".to_string(), 0, 4),
    //         ("chr1".to_string(), 6, 10),
    //     ]);
    //     let read_length = 4;
    //     let overwrite_output = true;
    //     let fastq_filename = "test_single";
    //     let paired_ended = false;
    //     let rng = NeatRng::new_from_seed(Vec::from([
    //         "Hello".to_string(),
    //         "cruel".to_string(),
    //         "world".to_string(),
    //     ]));
    //     let quality_score_model = QualityScoreModel::new();
    //     write_fastq(
    //         &mutated_fasta_map,
    //         &mut reads_dataset,
    //         read_length,
    //         overwrite_output,
    //         fastq_filename,
    //         paired_ended,
    //         quality_score_model.clone(),
    //         rng.clone()
    //     )
    //     .unwrap();
    //     let outfile1 = Path::new("test_single_r1.fastq");
    //     let outfile2 = Path::new("test_single_r2.fastq");
    //     assert!(outfile1.exists());
    //     assert!(!outfile2.exists());
    //     fs::remove_file(outfile1).unwrap();
    //
    //     let mut reads_dataset = Box::new(vec![
    //         ("chr1".to_string(), 0, 4),
    //         ("chr1".to_string(), 6, 14),
    //     ]);
    //
    //     let read_length = 4;
    //     let overwrite_output = false;
    //     let fastq_filename = "test_paired";
    //     let paired_ended = true;
    //
    //     write_fastq(
    //         &mutated_fasta_map,
    //         &mut reads_dataset,
    //         read_length,
    //         overwrite_output,
    //         fastq_filename,
    //         paired_ended,
    //         quality_score_model,
    //         rng.clone()
    //     )
    //         .unwrap();
    //
    //     let outfile1 = Path::new("test_paired_r1.fastq");
    //     let outfile2 = Path::new("test_paired_r2.fastq");
    //     assert!(outfile1.exists());
    //     assert!(outfile2.exists());
    //     fs::remove_file(outfile1).unwrap();
    //     fs::remove_file(outfile2).unwrap();
    // }

}
