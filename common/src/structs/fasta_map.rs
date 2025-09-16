//! This library contains low level structs designed to catalog and index the fasta file. 
//! Fasta files can be quite large and so this system is meant to speed up retreival when we
//! go to write output files, while keeping the burden of memory usage low, by writing the actual
//! sequences to file (encoded as Nucleotide) and retreiving them only as needed
//! We need to use the Variants struct to add variants to the contigs, making this a sort of second
//! level struct
use std::io;
use thiserror::Error;
use std::num::ParseIntError;
use std::{collections::HashMap, io::Write};
use serde::{Deserialize, Serialize};
use serde;
use serde_json;
use std::fs::{read_to_string, File};
use std::path::PathBuf;
use crate::structs::nucleotides::Nucleotide;
use log::error;

#[derive(Error, Debug)]
pub enum FastaMapError {
    #[error("FastaMap reported an IO Error: {0}")]
    IoError(#[from] io::Error),
    #[error("FastaMap reported bad coordinates.")]
    BadCoordinatesError,
    #[error("FastaMap reported an out of bounds error.")]
    OutOfBoundsError,
    #[error("FastaMap failed to validate the sequence map")]
    SeqMapValidationError,
    #[error("FastaMap reported a ContigNotFoundError")]
    ContigNotFoundError(&'static str),
    #[error("FastaMap reported a serde error: {0}")]
    SerdeValueError(#[from] serde::de::value::Error),
    #[error("FastaMap reported an error from serde_json: {0}")]
    SerdeJsonError(#[from] serde_json::Error),
    #[error("FastaMap reported a malformed filename error: {0}")]
    MalformedFilenameError(#[from] ParseIntError),
    #[error("Error retrieveing coordinates")]
    CoordinateRetrievalError,
}

pub fn decode_block_filename(path_string: &PathBuf) -> Result<(usize, usize), FastaMapError> {
    // This function parses the filename according to the established convention in NEAT
    // contigname_0001100000_0001200000.json
    // where the 2 groups of numbers are the start and end points of the sequence block
    // these coordinates are what we are trying to get from the name.
    let mut underscores: Vec<usize> = Vec::new();
    let mut index = 0;
    // We'll extract just the file stem part of the path and convert to string
    let filename = path_string.file_stem().unwrap().display().to_string();
    for char in filename.chars() {
        // Find all the underscores
        match char {
            '_' => {
                // Find the underscores
                underscores.push(index);
                index += 1; 
            },
            _ => index += 1,
        }
    }
    // the last two blocks: _000111000_000112000 contain the range
    let penultimate_underscore = underscores[underscores.len()-2];
    let ultimate_underscore = underscores[underscores.len()-1];
    // This should parse out to a usize position
    let start: usize = filename[(penultimate_underscore+1)..ultimate_underscore].parse()?;
    let end: usize = filename[(ultimate_underscore+1)..index].parse()?;
    Ok((start, end))
}

#[derive(Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct SequenceMap {
    // The purpose of this struct is to store where the N's were in the original fasta, since we
    // are overwriting those areas with random bases. This will help us place slightly more
    // meaningful mutations into the genome.
    pub region_type: RegionType,
    pub start: usize,
    pub end: usize,
}

impl SequenceMap {
    pub fn from(
        region_type: RegionType,
        start: usize,
        end: usize
    ) -> SequenceMap {
        SequenceMap {
            region_type,
            start,
            end,
        }
    }

    pub fn increment_end(&mut self, num: usize) {
        self.end += num;
    }
}

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct SequenceBlock {
    // This struct is for accessing the data stored in the sequence block file, a yaml file
    // The data file contains the following data:
    //   - contig: the key name from the contig list. This is the contig which this sequence is a part of
    //   - ref_start: the start point of this sequence relative to the full contig
    //   - ref_end: the end point of this sequence relative to the full contig
    //   - sequence: a list of Nucleotide data that encodes the sequence according to the pattern
    //         in the nucleotides package
    //   - sequence_map: A map of the region. This should allow us to quickly find areas to mutate.
    //        Note: this was moved from the contig level, because we were having to do a bunch of weird translations and
    //              remappings and so this ended up making more sense.
    pub contig: String,
    pub ref_start: usize,
    pub ref_end: usize,
    pub sequence: Vec<Nucleotide>,
    pub sequence_map: Vec<SequenceMap>,
}

impl SequenceBlock {
    pub fn from(mut filename: &PathBuf) -> Result<Self, FastaMapError> {
        // This is a check that the filename is properly encoded
        (_, _) = decode_block_filename(&filename)?;
        // Load the sequence block
        let json_data = read_to_string(&mut filename)?;
        let sequence_block: SequenceBlock = serde_json::from_str(&json_data)?;
        sequence_block.verify_regions()?;
        Ok(sequence_block)
    }

    pub fn get_start(&self) -> Result<usize, FastaMapError> {
        Ok(self.ref_start.clone())
    }
    
    pub fn get_end(&self) -> Result<usize, FastaMapError> {
        Ok(self.ref_end.clone())
    }

    pub fn contains(&self, (x, y): (usize, usize)) -> Result<bool, FastaMapError> {
        Ok((x >= self.get_start()?) && (y <= self.get_end()?))
    }

    pub fn get_non_n_regions(
        &self
    ) -> Result<Vec<&SequenceMap>, FastaMapError> {
        let mut return_regions = Vec::new();
        for region in &self.sequence_map {
            match region.region_type {
                RegionType::NonNRegion => return_regions.push(region),
                RegionType::NRegion => continue,
            }
        }
        Ok(return_regions)
    }

    fn verify_regions(&self) -> Result<(), FastaMapError> {
        // assume we have only N-regions and then look for valid regions.
        let mut all_ns = true;
        let mut prev_end = 0;
        for region in &self.sequence_map {
            match region.region_type {
                RegionType::NRegion => {
                    // We found a valid region
                    if all_ns { all_ns = false }
                    if prev_end > 0 {
                        // prev_end == 0 means we are in the first loop, where we skip this block
                        //
                        // If prev_end != segment.1 means are regions are either overlapping or not fully inclusive.
                        if prev_end != region.start {
                            return Err(FastaMapError::SeqMapValidationError)
                        }
                    };
                    // update prev_end for first and subsequent loops.
                    prev_end = region.end.clone();
                },
                RegionType::NonNRegion => {
                    prev_end = region.end.clone();
                },
            }
        };
        Ok(())
    }

    pub fn write_block(&mut self, mut file: File) -> io::Result<()> {
        // This will serialize the data to a byte vector and write it to the file
        let data = serde_json::to_vec(&self)?;
        file.write_all(&data)?;
        Ok(())
    }

    pub fn get_len(&self) -> usize {
        if self.ref_start >= self.ref_end {
            self.ref_start - self.ref_end
        } else {
            self.ref_end - self.ref_start
        }
    }

    pub fn get_seq_clone(&self) -> Result<Vec<Nucleotide>, FastaMapError> {
        Ok(self.sequence.clone())
    }

    pub fn get_subseq(&self, request_start: usize, request_end: usize) -> Result<Vec<Nucleotide>, FastaMapError> {
        // Fetches a subsequence of the original sequence
        if request_start >= request_end {
            error!("Bad coordinates for sequence block retrieval - start: {} >= end: {}", request_start, request_end);
            Err(FastaMapError::BadCoordinatesError)
        } else if (request_end <= self.ref_start) || 
                  (request_start > self.ref_end) || 
                  (request_end > self.ref_end) || 
                  (request_start < self.ref_start) {
            error!(
                "Requested coordinates out of bounds of sequence block -
                    request: ({}, {}), block ({}, {})",
                    request_start, request_end, self.ref_start, self.ref_end,
            );
            return Err(FastaMapError::OutOfBoundsError)
        } else {
            // start for this function is expected to be a coordinate from the overall contig,
            // so we modify by the known start and end point of this block to get the final
            // coordinates
            let block_start = request_start - self.ref_start;
            let block_end = request_end - self.ref_start;
            Ok(self.sequence[block_start..block_end].to_owned())
        }
    }
}

#[derive(Clone, Default, Debug, PartialEq, Serialize, Deserialize)]
pub enum RegionType {
    #[default]
    NRegion,
    NonNRegion,
}

#[derive(Debug, Clone)]
// Each Contig contains a vector of filenames, each containing a sequence block in json format.
// The contig level keeps things organized
pub struct Contig {
    // This is the short name of the contig
    pub name: String,
    // The overall length of the contig in bases
    pub contig_len: usize,
    // This is a vector of filenames. Each filename represents a sequence block and can be accessed via the contig block
    // The naming convention of the filenames will follow this format: contigname_0001100000_0001200000.jfa
    // Where contigname = the short name of the contig, the first string of digits is the start point of that file and the second block of digits is the end point 
    pub blocks: Vec<PathBuf>,
}

impl Contig {
    pub fn new(
        name: String,
        len: usize,
        blocks: Vec<PathBuf>,
    ) -> Result<Self, FastaMapError> {

        Ok(Contig {
            name,
            contig_len: len,
            blocks,
        })
    }
}

#[derive(Debug)]
// The FastaMap holds all the data for retrieving items from the fasta files and constructing outputs
pub struct FastaMap {
    // A vector of Contig objects, one for each contig in the fasta file
    pub contigs: Vec<Contig>,
    // This maps the short contig name, used throughout, to the full name from the fasta file
    pub name_map: HashMap<String, String>,
    // This is a Vec with the contig short names in the order they appeared in the fasta
    pub contig_order: Vec<String>,
}

impl FastaMap {
    pub fn new(
        contigs: Vec<Contig>, 
        name_map: HashMap<String, String>, 
        contig_order: Vec<String>
    ) -> Self {
        FastaMap { contigs, name_map, contig_order }
    }

    pub fn retrieve_contig(&self, request_name: String) -> Result<&Contig, FastaMapError> {
        for contig in self.contigs.iter() {
            if contig.name == request_name {
                return Ok(contig)
            }
        }
        Err(FastaMapError::ContigNotFoundError("Contig not found: {contig}"))
    }
}

#[cfg(test)]
mod tests {
    use std::fs::remove_file;
    use crate::structs::fasta_map::RegionType::NonNRegion;
    use crate::structs::nucleotides::Nucleotide::*;
    use super::*;

    #[test]
    fn test_sequence_block() {
        let mut seq_block = SequenceBlock {
            contig: "chrom1".to_string(),
            ref_start: 0,
            ref_end: 20,
            sequence: vec![A,A,A,A,A,A,A,A,A,A,A,A,T,T,T,A,A,A,T,A],
            sequence_map: vec![SequenceMap::from(NonNRegion, 0, 20)]
        };
        let filename = PathBuf::from("chrom1_000_020.json");
        // Test write works
        let file = File::create(&filename).unwrap();
        seq_block.write_block(file).unwrap();
        // let's read the file we just made and check the contents
        let seq_block_read = SequenceBlock::from(&filename).unwrap();
        assert_eq!(seq_block.contig, seq_block_read.contig);
        assert_eq!(seq_block.ref_start, seq_block_read.ref_start);
        assert_eq!(seq_block.ref_end, seq_block_read.ref_end);
        assert_eq!(seq_block.sequence, seq_block_read.sequence);
        remove_file(filename).unwrap();
    }

    #[test]
    fn test_filename_decoder() {
        let filename = PathBuf::from("chr1_001000_002000.json");
        let (a, b) = decode_block_filename(&filename).unwrap();
        assert_eq!(a, 1000);
        assert_eq!(b, 2000)
    }
    
    #[test]
    fn test_basic_map() {
        let filename = PathBuf::from("chrom1_0000_0020.json");
        let sequences = vec![filename];
        let contig = Contig {
            name: "chrom1".to_string(),
            contig_len: 20,
            blocks: sequences,
        };
        let name_map = HashMap::from([
            ("chrom1".to_string(), ">chrom1 foo bar\n".to_string())
        ]);
        let contigs = Vec::from([contig]);
        let contig_order = Vec::from(["chrom1".to_string()]);
        let fasta_map = FastaMap::new(
            contigs,
            name_map,
            contig_order,
        );

        let expected_output = "FastaMap { \
        contigs: [Contig { \
            name: \"chrom1\", \
            contig_len: 20, \
            blocks: [\
                \"chrom1_0000_0020.json\"] }], \
            name_map: {\"chrom1\": \">chrom1 foo bar\\n\"}, \
            contig_order: [\"chrom1\"] }";

        let test_out = format!("{:?}", fasta_map);

        assert_eq!(test_out, expected_output)
    }
}