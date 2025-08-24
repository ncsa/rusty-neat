use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io;
use tempfile::TempDir;
use common::file_tools::read_lines;
use common::structs::fasta_map::{FastaMap, SequenceBlock, Contig, RegionType};
use common::structs::fasta_map::RegionType::{NRegion, NonNRegion};
use common::structs::nucleotides::char_to_base;

pub const BUFSIZE: usize = 131_072usize; // i.e., 128kb
pub const DICT_SIZE: usize = 32768usize; // i.e., 32kb

pub fn read_fasta(
    fasta_path: &str,
    overlap_len: usize,
    temp_dir: &TempDir
) -> Result<Box<FastaMap>, io::Error> {
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
        panic!("fragment max {overlap_len} must be less than {BUFSIZE}.");
    }
    // Contigs are entire sequences in a fasta. This will vector will essentially be a map of the
    // fasta file.
    let mut contigs: Vec<Contig> = Vec::new();
    // Sequence blocks are objects holding data on a block of sequence information,
    // including a file where the actual data (encoded as a string of simple u8s) can be retrieved.
    // the start and end points of the block, relative to the reference, and the length of the
    // sequence (usually 131_072, but could be less).
    let mut sequence_blocks: Vec<String> = Vec::new();

    let mut name_map: HashMap<String, String> = HashMap::new();
    let mut contig_order: VecDeque<String> = VecDeque::new();
    let mut current_key = String::new();
    let mut buffer = [4_u8; BUFSIZE];

    let lines = read_lines(fasta_path)?;
    // This will help us keep track of where we are in the read
    let mut block_start = 0;
    let mut buffer_index = 0;
    let mut block_index = 0;
    // counts up the bases in the contig
    let mut contig_length = 0;
    for line in lines {
        match line {
            Ok(l) => {
                // indicates a header line in fasta format
                // If this is not the first loop, then this indicates
                // the end of the previous contig.
                if l.starts_with(">") {
                    // Checking if this is either the first loop or an empty read.
                    // After the first loop, then we need to log the previous contig
                    if contig_length == 0 && current_key.is_empty() {
                        // in this case, this is the first loop, skip this section
                    } else if contig_length == 0 && !current_key.is_empty() {
                        // in this case, we had an empty read
                        let name = name_map[&current_key].clone();
                        panic!("Read for {} is empty. Please remove from fasta file and try again", name)
                    } else {
                        // the total number of valid bases of this buffer. We subtract our overall contig
                        // length, which gives us how many actual bases were read, minus the start point
                        // of this block and that tells us the actual size of the buffer (since the tail
                        // will be filler, potentially, if the last block didn't fill the buffer)
                        let buffer_len = contig_length - block_start;
                        let buffer_to_write = buffer[..buffer_len].to_vec();

                        let block_file = process_buffer_into_sequenceblock(
                            &buffer_to_write, &current_key, block_start, contig_length, temp_dir
                        ).unwrap();

                        // Store the filename
                        sequence_blocks.push(block_file);

                        // Reset buffers, no need to save since we're starting a new contig
                        buffer = [4_u8; BUFSIZE];
                        buffer_index = 0;

                        // Add the contig to the list
                        contigs.push(Contig::new(
                            current_key.clone(),
                            contig_length.clone(),
                            sequence_blocks.clone(),
                            build_contig_map(&sequence_blocks).clone(),
                        ).unwrap());

                        // reset variables to start building the next contig and sequence block
                        sequence_blocks = Vec::new();
                        contig_length = 0;
                        block_index = 0;
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
                        if buffer_index >= BUFSIZE {
                            // grab the end of the buffer to create an overlap with the next section
                            let overlap_buffer: Vec<u8> = buffer[BUFSIZE-overlap_len..].iter().cloned().collect();
                            // Since our buffer is full, we can just set block end to the full length
                            let block_end = block_start + buffer.len();

                            let block_file = process_buffer_into_sequenceblock(
                                &buffer, &current_key, block_start, block_end, temp_dir
                            ).unwrap();

                            // record sequence block filename for later retrieval
                            sequence_blocks.push(block_file);

                            // Set start point for recording next block, offset to account for the overlap
                            block_start += BUFSIZE-overlap_len;
                            // Reset primary buffer and index
                            buffer = [4_u8; BUFSIZE];
                            buffer_index = 0;
                            for i in 0..overlap_len {
                                buffer[i] = overlap_buffer[i];
                                // keep the buffer index sync'd with the buffer
                                buffer_index += 1;
                            }
                            // Increment block index
                            block_index += 1;
                        }
                        buffer[buffer_index] = char_to_base(char);
                        buffer_index += 1;
                        contig_length += 1;
                    }
                }
            },
            _ => { panic!("Unreadable line in fasta file"); },
        }
    }
    // End of file reached, write the last sequence and finish the final contig,
    // same as above, but with no need to reset anything.
    let buffer_len = contig_length - block_start;
    let buffer_to_write = buffer[..buffer_len].to_vec();

    let block_file = process_buffer_into_sequenceblock(
        &buffer_to_write, &current_key, block_start, contig_length, temp_dir
    ).unwrap();

    // Store the filename
    sequence_blocks.push(block_file);

    // Add the contig to the list
    contigs.push(Contig::new(
        current_key.clone(),
        contig_length.clone(),
        sequence_blocks.clone(),
        build_contig_map(&sequence_blocks).clone(),
    ).unwrap());

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
    let delimiters = ['|', ' '];
    let mut delim_index = 0;
    let mut flag = false;
    for (index, char) in full_name.chars().enumerate() {
        if delimiters.contains(&char) {
            delim_index = index;
            flag = true;
            break;
        }
    };
    if !flag {
        delim_index = full_name.len();
    };
    full_name[..delim_index].to_string()
}

fn build_contig_map(sequence_blocks: &Vec<String>) -> Vec<(usize, usize, RegionType)> {
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

    for block_filename in sequence_blocks {
        // This will read the block data from the file
        let block = SequenceBlock::from(block_filename.to_string());
        // Now we process the block
        let block_sequence = block.sequence;
        for base in block_sequence {
            match base {
                4 => {
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
    }
    map
}

fn process_buffer_into_sequenceblock(
    buffer: &[u8], // The sequence in u8 format
    current_key: &str, // The contig short name
    block_start: usize, // Where, relative to the reference, this block starts
    block_end: usize, // We take this explicitly in case we have a less than full buffer
    temp_dir: &TempDir, // where to write files
) -> io::Result<String> {
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
        panic!("BUG: Sequence file name collision!")
    };

    let tmp_file = File::create(&sequence_file_path)?;
    sequence_block.write_block(tmp_file)?;
    
    Ok(sequence_file_path.display().to_string())
}

// Not sure I need the code below, but it might be useful later
// pub fn write_fasta(
//     mutated_fasta: &FastaMap,
//     overwrite_output: bool,
//     output_file: &str,
// ) -> io::Result<()> {
//     // mutated_fasta: the mutated fasta object, storing information about variants and file locations of the sequence
//     // overwrite_output: Boolean that controls whether an existing file is overwritten.
//     // output_file: the prefix for the output file name
//     const LINE_LENGTH: usize = 70; // fasta is usually writes with 70 characters per line
//     let mut output_fasta = format!("{}.fasta", output_file);
//     let mut outfile = File::open(&mut output_fasta)
//         .expect(&format!("Error opening {}", output_fasta));
//     let fasta_order = mutated_fasta.contig_order;
//     for contig in fasta_order {
//         let sequence = &mutated_fasta.contigs[contig];
//         let contig_length = sequence.len();
//         // Write contig name
//         writeln!(&mut outfile, ">{}", contig)?;
//         // write sequences[ploid] to this_fasta
//         let mut i = 0;
//         while i < contig_length {
//             let mut outlines = String::new();
//             // Generate a few lines before writing
//             'num_lines: for _ in 0..50 {
//                 for j in 0..LINE_LENGTH {
    
//                     // Check to make sure we are not out of bounds, which will happen if we
//                     // reach the end and there are not exactly 70 bases left (basically, every time)
//                     if (i + j) >= contig_length {
//                         // append the remaining bases
//                         for base in sequence[i..].iter() {
//                             outlines += &base.to_string();
//                         }
//                         outlines += "\n";
//                         // ensure outer loop breaks after writing final line, break inner loop.
//                         i += j;
//                         // todo this is missing a line break after the last full 70 chars.
//                         break 'num_lines
//                     }
//                     outlines += &sequence[i+j].to_string();
//                     i += 1
//                 }
//                 outlines += "\n";
//             }
//             writeln!(&mut outfile, "{}", outlines)?;
//         }
//     }
//     Ok(())
// }

#[cfg(test)]
mod tests {
    extern crate walkdir;
    use walkdir::WalkDir;
    use std::fs::{read_to_string, File};
    use super::*;

    #[test]
    fn test_process_buffer_into_sequenceblock() {
        let test = [0_u8; BUFSIZE];
        let temp_dir = tempfile::tempdir().unwrap();
        let segment_file = process_buffer_into_sequenceblock(
            &test, "chr_1", 25, 25+BUFSIZE, &temp_dir
        ).unwrap();

        // below is the decode code from fasta_map
        let sequence_block = SequenceBlock::from(segment_file);
        assert_eq!(Vec::from(test), sequence_block.sequence);
    }

    #[test]
    fn test_read_fasta() {
        let test_fasta = "test_data/H1N1.fa";
        let temp_dir = tempfile::tempdir().unwrap();
        let test_map = read_fasta(test_fasta, 350, &temp_dir).unwrap();
        assert_eq!(test_map.contig_order[0], "H1N1_HA".to_string());
        let mut flag = false;
        let mut contents = String::new();
        for entry in WalkDir::new(&temp_dir) {
            match entry {
                Ok(entry) => {
                    // if flag == false {
                    //     contents = read_to_string(entry.path())
                    //         .expect("Should have been able to read the file");
                    //     flag = true;
                    // }
                    println!("{}", entry.path().display())
                },
                Err(e) => eprintln!("Error: {}", e),
                }
            }
        // println!("With text:\n{contents}");
        temp_dir.close().unwrap();
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
        let sequence = "NNNNNNNNAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNN";
    }

}
