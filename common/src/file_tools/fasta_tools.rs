//! Various tools to work with the fasta structs in this library, which will allow NEAT 
//! to read and write fasta data. Currently, only the fasta reader is implemented.
//! The reader takes an input fasta file, a length of overlap (which, in practice, will be based
//! on the fragment length the simulator is using) and a temporary directory TempDir reference.
//! 
//! The fasta file contains the sequences that will be broken up into smaller pieces
//! read_fasta writes out json files to this temp dir to act as a database to facilitate 
//! the simulation
//! 
//! Eventually, we will add a write fasta function to write out output fasta files, but there's 
//! some simulation details to work out and it is low-priority.

use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io;
use tempfile::TempDir;
use crate::file_tools::file_io::read_lines;
use crate::structs::fasta_map::{Contig, FastaMap, FastaMapError, RegionType, SequenceBlock};
use crate::structs::fasta_map::RegionType::{NRegion, NonNRegion};
use crate::structs::nucleotides::Nucleotide;
use log::error;

pub const BUFSIZE: usize = 131_072usize; // i.e., 128kb
pub const DICT_SIZE: usize = 32768usize; // i.e., 32kb

#[derive(Debug)]
pub enum FastaToolsError {
    ReadFastaError,
    IoError(io::Error),
    FastaMapInterfaceError(FastaMapError),
    FastaBlockNameError,
}

impl From<FastaMapError> for FastaToolsError {
    fn from(error: FastaMapError) -> Self {
        FastaToolsError::FastaMapInterfaceError(error)
    }
}

impl From<io::Error> for FastaToolsError {
    fn from(error: io::Error) -> Self {
        FastaToolsError::IoError(error)
    }
}

pub fn read_fasta(
    fasta_path: &str,
    overlap_len: usize,
    temp_dir: &TempDir
) -> Result<Box<FastaMap>, FastaToolsError> {
    // Here's the basic idea:
    //     1. Read in 128kb worth of sequence (or the entire sequence)
    //     2. Recode the sequence as a series of u8 digits [0, 4]
    //        representing A, C, G, T, N, respectively
    //     3. Write the sequence to a file as a single line string of digits
    //          a. If we have not reached the end of the sequence, save the tail end (up to
    //             fragment max) from the end of the sequence to use as the beginning of the
    //             next sequence.
    //          b. Save the file name, coordinates, and necessary info to retrieve the contents
    //             in an object that will get moved around.
    //     4. Read the next block, if it's part of the same sequence, then there will be some
    //        overlap for overlapping reads (overlap determined by max fragment size).
    //     5. Each sequence longer than 128kb will be composed of multiple blocks
    //     6. The system will need a lot of I/O to account for this, and space for what could be a
    //        metric boatload of temporary files. Maybe BufSize could be variable.
    //     7. Incidentally, this will also facilitate multi-threading, which is a long goal.
    //
    // Reads a fasta file and turns it into a HashMap and puts it in the heap
    if overlap_len >= BUFSIZE {
        error!("fragment max {overlap_len} must be less than {BUFSIZE}.");
        return Err(FastaToolsError::ReadFastaError)
    }
    // Contigs are entire sequences in a fasta. This will vector will essentially be a map of the
    // fasta file.
    let mut contigs: Vec<Contig> = Vec::new();
    // This is a vector of filenames from which we can reconstruct the sequence_block.
    let mut block_filenames: Vec<String> = Vec::new();

    let mut name_map: HashMap<String, String> = HashMap::new();
    let mut contig_order: VecDeque<String> = VecDeque::new();
    let mut current_key = String::new();
    // We're using X solely as a placeholder. Any X's should be ignored or passed over by the code
    let mut buffer = [Nucleotide::X; BUFSIZE];

    let lines = read_lines(fasta_path)?;
    // This will help us keep track of where we are in the read
    let mut block_start = 0;
    let mut buffer_index = 0;
    // counts up the bases in the contig
    let mut contig_length = 0;
    for line in lines {
        match line {
            Ok(l) => {
                // In FASTA format, header lines start with ">" and then contain alphanumeric strings
                // with data about the machine that produced them, a way to sort them relative to other 
                // reads, and if this is not the first loop, then this indicates
                // the end of the previous contig. We check that condition with "new_data" and if
                // there is, we dump the data to file. 
                if l.starts_with(">") {
                    // Checking if this is either the first loop or an empty read.
                    // After the first loop, then we need to log the previous contig
                    if contig_length == 0 && current_key.is_empty() {
                        // in this case, this is the first loop, skip this section
                    } else if contig_length - block_start == 0 && !current_key.is_empty() {
                        // in this case, we had two contig names in a row and no sequence, which is an invalid FASTA
                        let name = name_map[&current_key].clone();
                        error!("Read for contig {} is empty. Please check file format.", name);
                        return Err(FastaToolsError::ReadFastaError)
                    } else {
                        // the total number of valid bases of this buffer. We subtract our overall contig
                        // length, which gives us how many actual bases were read (skips overlaps), minus 
                        // the overall start point of this block and that tells us the actual size of the 
                        // buffer (since the tail will be filler, potentially, if the last block didn't fill the buffer)
                        let buffer_len = contig_length - block_start;

                        // Note that we only want to write out the buffer if there was anything added since the last dump
                        if buffer_len > 0 {
                            // This is a buffer buffer to hold the overlap in the files. This is the unique
                            // thing we are doing that makes this slightly different from standard compression
                            // but it maybe that the costs aren't worth it in terms of performance
                            // If this idea isn't working, then maybe simply block-gzipping the files and 
                            // creating an index off the bgzip header will be better. We could probably load
                            // 2 bgzip blocks at a time, and cover the overlap that way. But it will make
                            // variant generation more complicated.
                            let buffer_to_write = buffer[..buffer_len].to_vec(); // may need to add 1

                            let block_file = process_buffer_into_sequenceblock(
                                &buffer_to_write, 
                                &current_key, 
                                block_start, 
                                &temp_dir
                            )?;

                            // Store the filename
                            block_filenames.push(block_file);
                        };
                        // New data was written, add block to list and reset
                        // Add the contig to the list
                        let contig_map = build_contig_map(&block_filenames)?;
                        contigs.push(Contig::new(
                            current_key,
                            contig_length,
                            block_filenames,
                            contig_map,
                        )?);

                        // Reset buffers, no need to save since we're starting a new contig
                        buffer = [Nucleotide::X; BUFSIZE];
                        buffer_index = 0;

                        // reset variables to start building the next contig and sequence block
                        block_filenames = Vec::new();
                        contig_length = 0;
                        block_start = 0;

                    }
                    
                    // grab the name from the string, omitting the initial '>' character
                    let long_name = l[1..].to_string();
                    current_key = extract_key_name(&long_name);
                    name_map.entry(current_key.clone()).or_insert(long_name);
                    contig_order.push_back(current_key.clone());
                } else {
                    for char in l.chars() {
                        // If we have filled the buffer, time to write the file and clear the buffer
                        if buffer_index == BUFSIZE {
                            // grab the end of the buffer to create an overlap with the next section
                            let overlap_buffer: Vec<Nucleotide> = buffer[buffer_index-overlap_len..].to_vec();
                            // Encode the data and write to temporary file
                            let sequence_block_filename = process_buffer_into_sequenceblock(
                                &buffer, 
                                &current_key, 
                                block_start, 
                                temp_dir
                            )?;

                            // record sequence block filename for later retrieval
                            block_filenames.push(sequence_block_filename);

                            // Set start point for recording next block, offset to account for the overlap
                            block_start += BUFSIZE-overlap_len;
                            // Reset primary buffer and index
                            buffer = [Nucleotide::X; BUFSIZE];
                            // We checked above that overlap_len < BUFSIZE, so we are safe from overflowing the buffer
                            buffer_index = 0;
                            for i in 0..overlap_len {
                                buffer[i] = overlap_buffer[i];
                                // increment, which will leave buffer_index = overlap len < BUFSIZE so we are safe from 
                                // overflowing the buffer.
                                buffer_index += 1;
                                // Note that we are specifically NOT incrementing the contig length, as these bases
                                // have already been counted.
                            }
                        }
                        buffer[buffer_index] = Nucleotide::from(char);
                        // We've written data, so switch the flag
                        buffer_index += 1;
                        contig_length += 1;
                    }
                }
            },
            _ => { 
                error!("Unreadable line in fasta file");
                return Err(FastaToolsError::ReadFastaError) 
            },
        }
    }
    // End of file reached, write the last sequence and finish the final contig,
    // same as above, but with no need to reset anything.
    let buffer_len = contig_length - block_start;
    if buffer_len > 0 {
        // if we have data to write...
        let buffer_len = contig_length - block_start;
        let buffer_to_write = buffer[..buffer_len].to_vec();

        let block_file = process_buffer_into_sequenceblock(
            &buffer_to_write, 
            &current_key, 
            block_start, 
            temp_dir,
        )?;

        // Store the filename
        block_filenames.push(block_file);
    }

    // any new data was written, now we add the final block to the contig
    // Add the contig to the list
    contigs.push(Contig::new(
        current_key.clone(),
        contig_length.clone(),
        block_filenames.clone(),
        build_contig_map(&block_filenames)?.clone(),
    )?);

    // Return box object with final FastaMap object, which will be used for processing
    Ok(Box::new(
        FastaMap {
            contigs,
            name_map,
            contig_order,
        }
    ))
}

fn extract_key_name(full_name: &str) -> String {
    // Fasta names can have a variety of formats. Attempt to parse the format.
    // may need to add more delimiters, but these are the most common.
    // This helper function can be expanded as needed in the future
    let delimiters = ['|', ' '];
    let mut delim_index = full_name.len();
    for (index, char) in full_name.chars().enumerate() {
        if delimiters.contains(&char) {
            // Look for a delimiter
            delim_index = index;
            break;
        }
    };
    // Returns full name if no delimiter found
    full_name[..delim_index].to_string()
}

fn build_contig_map(block_files: &Vec<String>) -> Result<Vec<(usize, usize, RegionType)>, FastaToolsError> {
    // This is a 2D vector map representing the regions of the string of DNA.
    // example: If a DNA strand looks like this:
    //     NNNNNNNNAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNN
    // Then the map would be (0, 10, NRegion), (10, 103, NonNRegion), (103, 110, NRegion)
    // Note that the numbers are invented for illustration
    // The idea here is to help ensure we are getting reads from regions with actual data, not
    // masked regions, since those regions tend to be highly repetitive
    let mut map: Vec<(usize, usize, RegionType)> = Vec::new();

    // tracking variables to create the map
    let mut region_start = 0;
    let mut region_end = 0;
    let mut inside_n_region = false;
    for filename in block_files.iter() {
        // This will read the block data from the file
        let block: SequenceBlock = SequenceBlock::from(filename)?;
        let first_base = block.sequence[0];
        // Many sequences start with unknown bases (telomeres), usually composed
        // of highly repetitive regions that are difficult to distinguish for the
        // sequencer and difficult to align later. So we will try to skip those.
        // Also, failing to check this was resulting in a first map element of (0, 0, NonNRegion)
        match first_base {
            Nucleotide::N => {
                // Process the N region
                // set the flag correctly to avoid issues later
                inside_n_region = true;
                region_end += 1
            },
            // If we find ourselves with a real base, then we are safe to continue processing
            _ => region_end += 1,
        }
        // Now we look at the remainder
        let remainder: Vec<Nucleotide> = block.sequence[1..].to_vec();
        // Now we process the block
        for base in remainder {
            match base {
                Nucleotide::N => {
                    if inside_n_region {
                        // We are already inside an n region, so we just need to increment
                        region_end += 1
                    } else {
                        // Welcome to an n region
                        inside_n_region = true;
                        // write previous region, which was non-n
                        map.push((region_start, region_end, NonNRegion));

                        region_start = region_end.clone();
                        region_end = region_start.clone() + 1;
                    }
                },
                _ => {
                    if inside_n_region {
                        // Welcome to a non-n region
                        inside_n_region = false;
                        // write previous region, an n-region
                        map.push((region_start, region_end, NRegion));
                        region_start = region_end.clone();
                        region_end = region_start.clone() + 1;
                    } else {
                        // we were already in a non-n region, so just increment
                        region_end += 1;
                    }
                },
            }
        }
        // Need to push that last block
        let mut region = NonNRegion;
        if inside_n_region {
            region = NRegion
        }
        map.push((region_start, region_end, region))
    }
    Ok(map)
}

fn process_buffer_into_sequenceblock(
    buffer: &[Nucleotide], // The sequence in u8 format
    current_key: &str, // The contig short name
    block_start: usize, // Where, relative to the reference, this block starts
    temp_dir: &TempDir, // where to write files
) -> Result<String, FastaToolsError> {
    // This function creates a temporary sequence block to store data then serializes it into
    // json format and writes that to file. Later, SequenceBlock can reconstruct this struct instance
    // based on the json file. This way we can parallelize the simulation.
    let block_end = block_start+buffer.len();
    let mut sequence_block = SequenceBlock {
        contig: current_key.to_string(),
        ref_start: block_start,
        ref_end: block_end,
        sequence: buffer.to_vec(),
    };

    let sequence_file_path = temp_dir
        .path()
        .join(format!("{}_{:010}_{:010}.json", current_key.to_string(), block_start, block_end));

    if sequence_file_path.exists() {
        error!("BUG: Sequence file name collision!");
        return Err(FastaToolsError::FastaBlockNameError)
    };

    let tmp_file = File::create(&sequence_file_path)?;
    sequence_block.write_block(tmp_file)?;
    
    Ok(sequence_file_path.display().to_string())
}

// Not sure I need the code below, but it might be useful later, and it's pretty standard in terms
// of fasta writing format, so the method will be the same whatever the details
pub fn write_fasta() -> Result<(), FastaToolsError> {
    todo!();
    // // mutated_fasta: the mutated fasta object, storing information about variants and file locations of the sequence
    // // overwrite_output: Boolean that controls whether an existing file is overwritten.
    // // output_file: the prefix for the output file name
    // const LINE_LENGTH: usize = 70; // fasta is usually writes with 70 characters per line
    // let mut output_fasta = format!("{}.fasta", output_file);
    // let mut outfile = File::open(&mut output_fasta)
    //     .expect(&format!("Error opening {}", output_fasta));
    // let fasta_order = mutated_fasta.contig_order;
    // for contig in fasta_order {
    //     let sequence = &mutated_fasta.contigs[contig];
    //     let contig_length = sequence.len();
    //     // Write contig name
    //     writeln!(&mut outfile, ">{}", contig)?;
    //     // write sequences[ploid] to this_fasta
    //     let mut i = 0;
    //     while i < contig_length {
    //         let mut outlines = String::new();
    //         // Generate a few lines before writing
    //         'num_lines: for _ in 0..50 {
    //             for j in 0..LINE_LENGTH {
    
    //                 // Check to make sure we are not out of bounds, which will happen if we
    //                 // reach the end and there are not exactly 70 bases left (basically, every time)
    //                 if (i + j) >= contig_length {
    //                     // append the remaining bases
    //                     for base in sequence[i..].iter() {
    //                         outlines += &base.to_string();
    //                     }
    //                     outlines += "\n";
    //                     // ensure outer loop breaks after writing final line, break inner loop.
    //                     i += j;
    //                     // todo this is missing a line break after the last full 70 chars.
    //                     break 'num_lines
    //                 }
    //                 outlines += &sequence[i+j].to_string();
    //                 i += 1
    //             }
    //             outlines += "\n";
    //         }
    //         writeln!(&mut outfile, "{}", outlines)?;
    //     }
    // }
    // Ok(())
}

#[cfg(test)]
mod tests {
    use walkdir::WalkDir;
    use super::*;

    #[test]
    fn test_process_buffer_into_sequenceblock() {
        let test = [Nucleotide::A; BUFSIZE];
        let temp_dir = tempfile::tempdir().unwrap();
        let segment_file = process_buffer_into_sequenceblock(
            &test, 
            "chr1", 
            25, 
            &temp_dir
        ).unwrap();

        // below is the decode code from fasta_map
        let sequence_block = SequenceBlock::from(&segment_file).unwrap();
        assert_eq!(Vec::from(test), sequence_block.sequence);
    }

    #[test]
    fn test_read_fasta() {
        let test_fasta = "test_data/H1N1.fa";
        let temp_dir = tempfile::tempdir().unwrap();
        let test_map = read_fasta(test_fasta, 350, &temp_dir).unwrap();
        assert_eq!(test_map.contig_order[0], "H1N1_HA".to_string());
        let dir_list: Vec<String> = WalkDir::new(&temp_dir)
            .into_iter()
            .map(|s| s.unwrap().path().display().to_string())
            .collect();
        let mut contains_ha = false;
        let mut contains_np = false;
        for item in dir_list {
            if item.contains(&"H1N1_HA".to_string()) {
                contains_ha = true
            } else if item.contains(&"H1N1_NP".to_string()) {
                contains_np = true
            }
        }
        assert!(contains_ha);
        assert!(contains_np);
    }

    #[test]
    fn test_extract_key_name() {
        let name1 = "BBB|AAA|CCC";
        let name2 = "BBB AAA CCC";
        let name3 = "BBBAAACCC";
        assert_eq!(extract_key_name(name1), "BBB".to_string());
        assert_eq!(extract_key_name(name2), "BBB".to_string());
        assert_eq!(extract_key_name(name3), "BBBAAACCC".to_string());
    }

    #[test]
    fn test_build_contig_map() {
        let sequence = String::from("NNNNNNNNAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNN");
        let mut buffer = Vec::new();
        for char in sequence.chars() {
            buffer.push(
                Nucleotide::from(char)
            )
        }
        let current_key = "chr1";
        let block_start = 0;
        let temp_dir =  tempfile::tempdir().unwrap();
        let name = process_buffer_into_sequenceblock(
            &buffer, 
            current_key, 
            block_start, 
            &temp_dir,
        ).unwrap();
        let blocks = vec![name];
        let map: Vec<(usize, usize, RegionType)> = build_contig_map(&blocks).unwrap();
        assert_eq!(map.len(), 3);
        assert_eq!(map[0], (0, 8, RegionType::NRegion));
        assert_eq!(map[1], (8, 59, RegionType::NonNRegion));
        assert_eq!(map[2], (59, buffer.len(), RegionType::NRegion));
    }

}
