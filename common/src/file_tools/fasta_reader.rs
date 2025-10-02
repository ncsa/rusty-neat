//! Various tools to work with the fasta structs in this library, which will allow NEAT 
//! to read and write fasta data. Currently, only the fasta reader is implemented.
//! The reader takes an input fasta file, a length of overlap (which, in practice, will be based
//! on the fragment length the simulator is using) and a temporary directory TempDir reference.
//! 
//! The fasta file contains the sequences that will be broken up into smaller pieces
//! read_fasta writes out json files to this temp dir to act as a database to facilitate 
//! the simulationuse std::collections::HashMap;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use log::*;
use std::path::PathBuf;
use tempfile::TempDir;
use crate::file_tools::file_io::read_lines;
use crate::structs::{
    fasta_map::{
        Contig,
        FastaMap,
        FastaMapError,
        SequenceBlock,
        SequenceMap,
        RegionType::{
            NRegion,
            NonNRegion
        },
    },
    nucleotides::{
        Nucleotide,
        Nucleotide::{N, X, Maskeda, Maskedc, Maskedg, Maskedt},
        NucleotideSelector,
    },
};
use simple_rng::{NeatRng, NeatRngError};
use thiserror::Error;

pub const BUFSIZE: usize = 131_072usize; // i.e., 128kb
pub const DICT_SIZE: usize = 32768usize; // i.e., 32kb

#[derive(Debug, Error)]
pub enum FastaReaderError {
    #[error("Error reading fasta file")]
    ReadFastaError,
    #[error("IO error while reading fasta file")]
    IoError(#[from] io::Error),
    #[error("Error interfacing with FastaMap: {0}")]
    FastaMapInterfaceError(#[from] FastaMapError),
    #[error("Error in fasta block name")]
    FastaBlockNameError,
    #[error("Error in random number generator: {0}")]
    RngError(#[from] NeatRngError),
}

pub fn read_fasta(
    fasta_path: &PathBuf,
    nucleotide_selector: Option<&NucleotideSelector>,
    overlap_len: usize,
    temp_dir: &TempDir,
    mut rng: Option<&mut NeatRng>,
) -> Result<Box<FastaMap>, FastaReaderError> {
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
        return Err(FastaReaderError::ReadFastaError)
    }
    // Contigs are entire sequences in a fasta. This will vector will essentially be a map of the
    // fasta file.
    let mut contigs: Vec<Contig> = Vec::new();
    // This is a vector of filenames from which we can reconstruct the sequence_block.
    let mut block_filenames: Vec<PathBuf> = Vec::new();
    // This will keep our sequence maps stored at a higher level, so that we don't have to
    // read the entire sequence from file just to know where the N-regions are
    let mut block_sequence_maps: HashMap<PathBuf, Vec<SequenceMap>> = HashMap::new();
    // This hashmap maps the short name to a description.
    let mut name_map: HashMap<String, Option<String>> = HashMap::new();
    // This will give us the write order for the bam and vcf
    let mut contig_order: Vec<String> = Vec::new();
    // current contig
    let mut current_key = String::new();
    // We're using X solely as a placeholder. Any X's should be ignored or passed over by the code
    let mut buffer = [X; BUFSIZE];
    // Bufreader for the fasta file
    let lines = read_lines(fasta_path)?;
    // This will help us keep track of where we are in the read
    let mut block_start = 0;
    // And this helps us keep track of where we are in the buffer
    let mut buffer_index = 0;
    // counts up the bases in the contig
    let mut contig_length = 0;
    // Process the file
    for line in lines {
        match line {
            Ok(line_str) => {
                // In FASTA format, header lines start with ">" and then contain alphanumeric strings
                // with data about the machine that produced them, a way to sort them relative to other 
                // reads, and if this is not the first loop, then this indicates
                // the end of the previous contig. We check that condition with "new_data" and if
                // there is, we dump the data to file. These types of things always have a first loop/
                // last loop problem. Happy to take suggestions on how to make this less clunky.
                if line_str.starts_with(">") {
                    // Checking if this is either the first loop or an empty read.
                    // After the first loop, then we need to log the previous contig
                    if contig_length == 0 && current_key.is_empty() {
                        // in this case, this is the first loop, skip this section
                    } else if contig_length - block_start == 0 && !current_key.is_empty() {
                        // If skip_contig is true, then we don't throw an error
                        // in this case, we had two contig names in a row and no sequence, which is an invalid FASTA
                        let desc = name_map[&current_key].clone();
                        let mut name = current_key.to_owned();
                        match desc {
                            // reinsert space
                            Some(description) => name.push_str(&format!(" {}", &description)),
                            None => {
                                // name = current_key
                            },
                        }
                        error!("Read for contig {} is empty. Please check file for proper FASTA format.", name);
                        return Err(FastaReaderError::ReadFastaError)
                    } else {
                        // We have read some data and are ready to dump this buffer
                        //
                        // the total number of valid bases of this buffer. We subtract our overall contig
                        // length, which gives us how many actual bases were read (skips overlaps), minus 
                        // the overall start point of this block and that tells us the actual size of the 
                        // buffer (since the tail will be filler, potentially, if the last block didn't fill the buffer)
                        let buffer_len = contig_length - block_start;
                        // This is a buffer buffer to hold the overlap in the files. This is the unique
                        // thing we are doing that makes this slightly different from standard compression
                        // but it maybe that the costs aren't worth it in terms of performance
                        // If this idea isn't working, then maybe simply block-gzipping the files and 
                        // creating an index off the bgzip header will be better. We could probably load
                        // 2 bgzip blocks at a time, and cover the overlap that way. But it will make
                        // variant generation more complicated.
                        let buffer_to_write = buffer[..buffer_len].to_vec(); // may need to add 1

                        let (block_file, sequence_maps) = process_buffer_into_sequenceblock(
                            &buffer_to_write, 
                            &current_key, 
                            block_start, 
                            &temp_dir,
                            nucleotide_selector,
                            &mut rng,
                        )?;

                        // Store the filename
                        block_filenames.push(block_file.clone());

                        // Store the sequence_maps
                        block_sequence_maps.insert(block_file, sequence_maps);

                        let desc = name_map[&current_key].clone();
                        // Add the contig to the list
                        contigs.push(
                            Contig::new(
                                &current_key.clone(),
                                &desc,
                                contig_length,
                                block_filenames,
                                block_sequence_maps,
                            )?
                        );

                        // Reset buffers, no need to save since we're starting a new contig
                        buffer = [X; BUFSIZE];
                        buffer_index = 0;
                        // reset variables to start building the next contig and sequence block
                        block_filenames = Vec::new();
                        block_sequence_maps = HashMap::new();
                        contig_length = 0;
                        block_start = 0;
                    }
                    // grab the name from the string, omitting the initial '>' character
                    let long_name = line_str[1..].trim_end().to_string();
                    let (id, desc) = extract_key_name(&long_name);
                    current_key = id;
                    if name_map.contains_key(&current_key) {
                        // name collision
                        let desc_orig = name_map[&current_key].clone();
                        let mut orig_name = current_key.clone();
                        match desc_orig {
                            // reinsert space
                            Some(desc) => orig_name.push_str(&format!(" {}", desc)),
                            None => {
                                // Do nothing
                            },
                        }
                        error!("Name collision in fasta file, two entries with identical key names: {current_key} from {} {}", orig_name, long_name);
                        return Err(FastaReaderError::ReadFastaError)
                    }
                    name_map.insert(current_key.clone(), desc);
                    contig_order.push(current_key.clone());

                } else {
                    for char in line_str.chars() {
                        // If we have filled the buffer, time to write the file and clear the buffer
                        if buffer_index == BUFSIZE {
                            // grab the end of the buffer to create an overlap with the next section
                            let overlap_buffer: Vec<Nucleotide> = buffer[buffer_index-overlap_len..].to_vec();
                            // Encode the data and write to temporary file
                            let (sequence_block_filename, sequence_maps) = process_buffer_into_sequenceblock(
                                &buffer,
                                &current_key,
                                block_start,
                                temp_dir,
                                nucleotide_selector,
                                &mut rng,
                            )?;

                            // record sequence block filename for later retrieval
                            block_filenames.push(sequence_block_filename.clone());

                            // Store the sequence_maps
                            block_sequence_maps.insert(sequence_block_filename, sequence_maps);

                            // Set start point for recording next block, offset to account for the overlap
                            block_start += BUFSIZE-overlap_len;
                            // Reset primary buffer and index
                            buffer = [X; BUFSIZE];
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
                        buffer_index += 1;
                        contig_length += 1;
                    }
                }
            },
            _ => { 
                error!("Unreadable line in fasta file");
                return Err(FastaReaderError::ReadFastaError) 
            },
        }
    }
    if contig_order.is_empty() {
        error!("There were no contigs to read.");
        debug!("Check fasta for proper format.");
        return Err(FastaReaderError::ReadFastaError)
    }
    // End of file reached, write the last sequence and finish the final contig,
    // same as above, but with no need to reset anything.
    let buffer_len = contig_length - block_start;
    if buffer_len > 0 {
        // if we have data to write...
        let buffer_len = contig_length - block_start;
        let buffer_to_write = buffer[..buffer_len].to_vec();

        let (block_file, sequence_maps) = process_buffer_into_sequenceblock(
            &buffer_to_write, 
            &current_key, 
            block_start, 
            temp_dir,
            nucleotide_selector,
            &mut rng,
        )?;

        // Store the filename
        block_filenames.push(block_file.clone());

        // Store the sequence_maps
        block_sequence_maps.insert(block_file, sequence_maps);    
    }

    // any new data was written, now we add the final block to the contig
    // Add the contig to the list
    contigs.push(
        Contig::new(
            &current_key,
            &name_map[&current_key],
            contig_length,
            block_filenames,
            block_sequence_maps,
        )?
    );

    // Return box object with final FastaMap object, which will be used for processing
    Ok(Box::new(
        FastaMap {
            contigs,
            contig_order,
        }
    ))
}

fn extract_key_name(full_name: &str) -> (String, Option<String>) {
    // Fasta names can have a variety of formats. Our assumption here is that 
    // everything before the first space is a unique ID field, and everything after is 
    // a description. Oftentimes, the "common" name ends up in the description, which
    // can cause problems in the vcf. We may need to add a way to search both id
    // and description.
    let mut fields = full_name.splitn(2, char::is_whitespace);
    let id = fields.next().map(|s| s.to_owned()).unwrap();
    let desc = fields.next().map(|s| s.to_owned());
    (id, desc)
}

fn map_buffer(sequence: &[Nucleotide]) -> Result<Vec<SequenceMap>, FastaReaderError> {
    // This is a 2D vector map representing the regions of the string of DNA.
    // example: If a DNA strand looks like this:
    //     NNNNNNNNAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNN
    // Then the map would be (NRegion, 0, 10), (NonNRegion, 10, 103), (NRegion, 103, 110)
    // Note that the numbers are invented for illustration
    // The idea here is to help ensure we are getting reads from regions with actual data, not
    // masked regions, since those regions tend to be highly repetitive
    let mut map: Vec<SequenceMap> = Vec::new();

    // tracking variables to create the map
    let mut region_start = 0;
    let mut region_end = 1;
    let mut inside_n_region = false;
    let first_base = sequence[0];
    // Many sequences start with unknown bases (telomeres), usually composed
    // of highly repetitive regions that are difficult to distinguish for the
    // sequencer and difficult to align later. These can be represented by N or sometimes lower
    // case bases. So we will try to skip those.
    // Also, failing to check this was resulting in a first map element of (0, 0, NonNRegion)
    match first_base {
        N | X | Maskeda | Maskedc | Maskedg | Maskedt => {
            // Process the N region
            // set the flag correctly to avoid issues later
            inside_n_region = true;
            for base in &sequence[1..] {
                match base {
                    N | X | Maskeda | Maskedc | Maskedg | Maskedt => {
                        region_end += 1;
                    },
                    _ => {
                        inside_n_region = false;
                        map.push(SequenceMap::from(NRegion, 0, region_end));
                        region_start = region_end.clone();
                        region_end = region_start + 1;
                        break
                    }
                }
            }
            if inside_n_region {
                // This indicates our entire block was Ns
                map.push(SequenceMap::from(NRegion, 0, region_end));
                return Ok(map)
            }
        },
        // If we find ourselves with a real base, then we are safe to continue processing
        _ => {
            // Nothing to do in this case
        },
    }
    // Now we look at the remainder
    let remainder: Vec<Nucleotide> = sequence[region_end..].to_vec();
    // Now we process the block
    for base in remainder {
        match base {
            N | X | Maskeda | Maskedc | Maskedg | Maskedt => {
                if inside_n_region {
                    // We are already inside an n region, so we just need to increment
                    region_end += 1
                } else {
                    // Welcome to an n region
                    inside_n_region = true;
                    // write previous region, which was non-n
                    map.push(SequenceMap::from(NonNRegion, region_start, region_end));

                    region_start = region_end.clone();
                    region_end = region_start.clone() + 1;
                }
            },
            _ => {
                if inside_n_region {
                    // Welcome to a non-n region
                    inside_n_region = false;
                    // write previous region, an n-region
                    map.push(SequenceMap::from(NRegion, region_start, region_end));
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
    map.push(SequenceMap::from(region, region_start, region_end));
    Ok(map)
}

fn process_buffer_into_sequenceblock(
    buffer: &[Nucleotide],
    current_key: &str,
    block_start: usize,
    temp_dir: &TempDir,
    nucleotide_selector: Option<&NucleotideSelector>,
    rng: &mut Option<&mut NeatRng>,
) -> Result<(PathBuf, Vec<SequenceMap>), FastaReaderError> {
    // This function creates a temporary sequence block to store data then serializes it into
    // json format and writes that to file. Later, SequenceBlock can reconstruct this struct instance
    // based on the json file. This way we can parallelize the simulation.
    // params:
    //   buffer: The sequence to write in u8 format
    //   current_key: The contig short name
    //   block_start: Where, relative to the reference, this block starts
    //   temp_dir: where to write temp files
    //   nucleotide_selector: Optional for producing random bases
    //   rng: Optional for producing random bases
    let block_end = block_start+buffer.len();
    let sequence_map = map_buffer(buffer)?;
    let mut buffer_to_write: Vec<Nucleotide> = Vec::with_capacity(buffer.len());
    match nucleotide_selector {
        Some(model) => {
            match rng {
                Some(random_num) => {
                    for item in buffer {
                        match item {
                            N => {
                                buffer_to_write.push(model.sample_bases(random_num.random()?).get_masked());
                            },
                            X => {
                                //skip
                            }
                            _ => {
                                buffer_to_write.push(item.clone());
                            },
                        }
                    }
                },
                None => return Err(FastaReaderError::ReadFastaError),
            }
        },
        None => {
            for item in buffer {
                match item {
                    X => {
                        //skip
                    }
                    _ => {
                        buffer_to_write.push(item.clone());
                    },
                }
            }
        }
    }

    let mut sequence_block = SequenceBlock {
        contig: current_key.to_string(),
        ref_start: block_start,
        ref_end: block_end,
        sequence: buffer_to_write,
        sequence_map: sequence_map.clone(),
    };

    let sequence_file_path = temp_dir
        .path()
        .join(format!("{}_{:010}_{:010}.json", current_key.to_string(), block_start, block_end));

    if sequence_file_path.exists() {
        error!("BUG: Sequence file name collision!");
        return Err(FastaReaderError::FastaBlockNameError)
    };

    let tmp_file = File::create(&sequence_file_path)?;
    sequence_block.write_block(tmp_file)?;
    
    Ok((sequence_file_path, sequence_map))
}

pub fn count_fasta(filename: &PathBuf) -> Result<usize, FastaReaderError> {
    let open_file = read_lines(filename)?;
    let mut count = 0;
    for raw_line in open_file {
        match raw_line {
            Ok(l) => {
                if l.starts_with(">") {
                    count += 1;
                }
            },
            Err(error) => {
                error!("Error reading fasta file! {:?}", error);
                return Err(FastaReaderError::ReadFastaError)
            }
        }

    }
    Ok(count)
}

#[cfg(test)]
mod tests {
    use walkdir::WalkDir;
    use crate::structs::fasta_map::RegionType;

    use super::*;

    #[test]
    fn test_count_fasta() {
        // TODO add test
        assert!(true)
    }

    #[test]
    fn test_process_buffer_into_sequenceblock() {
        let test = [Nucleotide::A; BUFSIZE];
        let temp_dir = tempfile::tempdir().unwrap();
        let nucleotide_selector = NucleotideSelector::new();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let (segment_file, sequence_maps) = process_buffer_into_sequenceblock(
            &test, 
            "chr1", 
            25, 
            &temp_dir,
            Some(&nucleotide_selector),
            &mut Some(&mut rng),
        ).unwrap();

        // below is the decode code from fasta_map
        let sequence_block = SequenceBlock::from(&segment_file).unwrap();
        assert_eq!(sequence_block.sequence, Vec::from(test));
        assert_eq!(sequence_maps[0], SequenceMap::from(RegionType::NonNRegion, 0, BUFSIZE))
    }

    #[test]
    fn test_read_fasta() {
        let test_fasta = PathBuf::from("test_data/H1N1.fa");
        let temp_dir = tempfile::tempdir().unwrap();
        let nucleotide_selector = NucleotideSelector::new();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let test_map = read_fasta(
            &test_fasta,
            Some(&nucleotide_selector),
            350,
            &temp_dir,
            Some(&mut rng),
        ).unwrap();
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
        assert_eq!(extract_key_name(name1), ("BBB|AAA|CCC".to_string(), None));
        assert_eq!(extract_key_name(name2), ("BBB".to_string(), Some("AAA CCC".to_string())));
        assert_eq!(extract_key_name(name3), ("BBBAAACCC".to_string(), None));
    }

    #[test]
    fn test_build_contig_map_mixed() {
        let sequence = String::from("NNNAAAANNN");
        let mut buffer = Vec::new();
        for char in sequence.chars() {
            buffer.push(Nucleotide::from(char))
        }
        let map: Vec<SequenceMap> = map_buffer(buffer.as_slice()).unwrap();
        assert_eq!(map.len(), 3);
        assert_eq!(map[0], SequenceMap::from(NRegion, 0, 3));
        assert_eq!(map[1], SequenceMap::from(NonNRegion, 3, 7));
        assert_eq!(map[2], SequenceMap::from(NRegion, 7, buffer.len()));
    }

    #[test]
    fn test_build_contig_map_alln() {
        let sequence = String::from("NNNNNNNN");
        let mut buffer = Vec::new();
        for char in sequence.chars() {
            buffer.push(Nucleotide::from(char))
        }
        let map: Vec<SequenceMap> = map_buffer(buffer.as_slice()).unwrap();
        assert_eq!(map.len(), 1);
        assert_eq!(map[0], SequenceMap::from(NRegion, 0, 8));
    }

    #[test]
    fn test_build_contig_map_middle() {
        let sequence = String::from("AAAANNNAAAAAAAAAAAA");
        let mut buffer = Vec::new();
        for char in sequence.chars() {
            buffer.push(Nucleotide::from(char))
        }
        let map: Vec<SequenceMap> = map_buffer(buffer.as_slice()).unwrap();
        assert_eq!(map.len(), 3);
        assert_eq!(map[0], SequenceMap::from(NonNRegion, 0, 4));
        assert_eq!(map[1], SequenceMap::from(NRegion, 4, 7));
        assert_eq!(map[2], SequenceMap::from(NonNRegion, 7, buffer.len()));
    }
}
