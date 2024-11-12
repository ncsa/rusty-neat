use std::collections::{HashMap, VecDeque};
use std::ffi::OsString;
use std::fs::{read_to_string, File};
use std::io;
use std::hash::Hash;
use std::io::{BufRead, BufReader, ErrorKind, Write};
use common::file_tools::read_lines;
use common::structs::fasta_map::{FastaMap, SequenceBlock, Contig, RegionType};
use tempfile::TempDir;
use common::structs::fasta_map::RegionType::{NRegion, NonNRegion};
use common::structs::nucleotides::char_to_base;

pub const BUFSIZE: usize = 131_072usize; // i.e., 128kb
pub const DICT_SIZE: usize = 32768usize; // i.e., 32kb

pub fn read_fasta(
    fasta_path: &str,
    fragment_max: usize,
    temp_dir: TempDir
) -> Result<(Box<FastaMap>), io::Error> {
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
    if fragment_max >= BUFSIZE {
        panic!("fragment max {fragment_max} must be less than {BUFSIZE}. Possible bug.");
    }
    // Contigs are entire sequences in a fasta. This will vector will essentially be a map of the
    // fasta file.
    let mut contigs: Vec<Contig> = Vec::new();
    // Sequence blocks are objects holding data on a block of sequence information,
    // including a file where the actual data (encoded as a string of simple u8s) can be retrieved.
    // the start and end points of the block, relative to the reference, and the length of the
    // sequence (usually 131_072, but could be less).
    let mut sequence_blocks: Vec<SequenceBlock> = Vec::new();

    let mut name_map: HashMap<String, String> = HashMap::new();
    let mut contig_order: VecDeque<String> = VecDeque::new();
    let mut current_key = String::new();
    let mut buffer = [4_u8; BUFSIZE];
    let mut next_buffer = [4_u8; BUFSIZE];

    // Gather stats about the fasta file for sanity checking
    let stats = GatherStats::process_fasta(fasta_path);

    let lines = read_lines(fasta_path)?;
    // This will help us keep track of where we are in the read
    let mut block_start = 0;
    let mut local_index = 0;
    let mut block_index = 0;
    // to record the actual padding length
    let mut padding_length = 0;
    // counts up the bases in the contig
    let mut contig_length = 0;
    for line in lines {
        match line {
            Ok(l) => {
                // indicates a header line in fasta format
                // If this is not the first loop, then this indicates
                // the end of the previous contig.
                if l.starts_with(">") {
                    // if this is not the first loop, then we need to log the previous contig
                    if contigs.len() != 0 {
                        // the total number of valid bases of this buffer
                        let num_bases = contig_length - block_start;
                        // We still need to create the temp file and write the sequence to it.
                        // we could turn this section down to the for loop into a function
                        let (num_bases, file_path) = write_buffer_to_file(
                            &buffer,
                            contig_length,
                            block_start,
                            &temp_dir,
                            &current_key,
                            block_index,
                        )?;

                        // record sequence block info
                        sequence_blocks.push(SequenceBlock::new(
                            block_start.clone(),
                            num_bases.clone(),
                            file_path,
                        ).unwrap());

                        // Reset buffers
                        buffer = next_buffer.clone();
                        next_buffer = [4_u8; BUFSIZE];
                        // Reset local index to synchronize with new buffer on a new contig
                        local_index = 0;

                        let contig_map = build_contig_map(&sequence_blocks);

                        // Add the contig to the list
                        contigs.push(Contig::new(
                            current_key.clone(),
                            contig_length.clone(),
                            sequence_blocks.clone(),
                            contig_map.clone(),
                        ).unwrap());

                        // reset variables to start building the next contig and sequence block
                        sequence_blocks = Vec::new();
                        contig_length = 0;
                        block_index = 0;
                        block_start = 0;
                        // Fasta files always end on a line break or EOF, so
                        // we do not need to worry about swapping buffers or anything.
                        // Update buffer and next buffer
                        buffer = next_buffer.clone();
                        next_buffer = [4_u8; BUFSIZE];
                    }
                    // grab the name from the string, omitting the initial '>' character
                    let long_name = l[1..].to_string();
                    current_key = extract_key_name(&long_name);
                    name_map.entry(current_key.clone()).or_insert(long_name);
                } else {
                    for char in l.chars() {
                        // Catch any potential index out of range problems here
                        if local_index >= BUFSIZE {
                            let tail = buffer[BUFSIZE-fragment_max..];
                            // overwrite the beginning of the next buffer with the end of this one
                            // to create an overlap
                            padding_length = tail.len();
                            for i in 0..padding_length {
                                next_buffer[i] = tail[i];
                            }
                            let (num_bases, file_path) = write_buffer_to_file(
                                &buffer,
                                contig_length,
                                block_start,
                                &temp_dir,
                                &current_key,
                                block_index,
                            )?;
                            // record sequence block info
                            sequence_blocks.push(SequenceBlock::new(
                                block_start.clone(),
                                num_bases.clone(),
                                file_path.clone(),
                            ).unwrap());

                            // Set start point for recording next block
                            block_start += BUFSIZE-padding_length;
                            // Reset buffers
                            buffer = next_buffer.clone();
                            next_buffer = [4_u8; BUFSIZE];
                            // Reset local index to synchronize with new buffer (0) + padding_length
                            // which starts adding bases from the sequence after the padding
                            local_index = padding_length;
                            block_index += 1
                        }
                        *buffer[local_index] = *char_to_base(char);
                        local_index += 1;
                        contig_length += 1;
                    }
                }
            },
            _ => { panic!("AAAAAAAA"); },
        }
    }
    // End of file reached, write the last sequence and finish the final contig,
    // same as above, but with no need to reset anything.
    let (num_bases, file_path) = write_buffer_to_file(
        &buffer,
        contig_length,
        block_start,
        &temp_dir,
        &current_key,
        block_index,
    )?;

    // record sequence block info
    sequence_blocks.push(SequenceBlock::new(
        block_start.clone(),
        num_bases.clone(),
        file_path.clone(),
    ).unwrap());

    let contig_map = build_contig_map(&sequence_blocks);

    // Add the contig to the list
    contigs.push(Contig::new(
        current_key.clone(),
        contig_length.clone(),
        sequence_blocks.clone(),
        contig_map.clone(),
    ).unwrap());

    // Return box object with final FastaMap object, which will be used for processing
    Ok((Box::new(
        FastaMap {
            contigs,
            name_map,
            contig_order,
        }
    )))
}

fn extract_key_name(key: &str) -> String {
    // Fasta names can have a variety of formats. Attempt to parse the format.
    // may need to add more delimiters, but these are the most common.
    let delimiters = ['|', ' '];
    let mut delim_index = 0;
    let mut flag = false;
    for (index, char) in key.chars().enumerate() {
        if delimiters.contains(&char) {
            delim_index = index;
            flag = true;
            break;
        }
    };
    if !flag {
        delim_index = key.len();
    };
    key[..delim_index].to_string()

}

fn build_contig_map(sequence_blocks: &Vec<SequenceBlock>) -> Vec<(usize, usize, RegionType)> {
    // This is a 2D map representing the regions of the string of DNA.
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

    for mut block in sequence_blocks {
        let block_sequence = block.retrieve_all().unwrap();
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

fn write_buffer_to_file(
    buffer: &[u8; BUFSIZE],
    contig_len: usize,
    start: usize,
    temp_dir: &TempDir,
    key: &str,
    index: usize
) -> Result<(usize, OsString), io::Error> {
    let num_bases = contig_len - start;
    let sequence_file_path = temp_dir
        .path()
        .join(format!("{}_{:010}.seq", &key, index));
    if sequence_file_path.exists() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            "sequence file already exists (Unexpected!)"))
    };
    let mut tmp_file = File::create(&sequence_file_path)?;
    for i in 0..num_bases {
        tmp_file.write(buffer[i].as_bytes())?;
    };
    Ok((num_bases, sequence_file_path.into_os_string()))
}

struct GatherStats {
    // how many sequence lines the file contains
    line_count: usize,
    // Number of contigs, to verify
    number_of_contigs: usize,
    // Whole number times this genome will fill the buffer
    approximate_genome_size: usize,
}

impl GatherStats {
    fn process_fasta(filename: &str) -> Result<GatherStats, io::Error> {
        let mut line_count = 0;
        let mut number_of_contigs = 0;
        let mut line_width = 0;

        let file = read_lines(filename)?;
        for line in file {
            match line {
                Ok(l) => {
                    if l.starts_with(">") {
                        number_of_contigs += 1;
                    } else if line_width == 0 {
                        line_width = l.len();
                    }
                    line_count += 1
                },
                Some(Err(E)) => return Err(E),
                _ => panic!("Unknown fasta read error"),
            }
        };

        Ok(GatherStats {
            line_count,
            number_of_contigs,
            approximate_genome_size:(line_width * line_count)/BUFSIZE,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_fasta() {
        // let test_fasta = "test_data/H1N1.fa";
        // let (test_map, map_order) =
        //     read_fasta(test_fasta, 350).unwrap();
        // assert_eq!(map_order[0], "H1N1_HA".to_string());
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
