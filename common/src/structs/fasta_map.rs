//! This library contains low level structs designed to catalog and index the fasta file. 
//! Fasta files can be quite large and so this system is meant to speed up retreival when we
//! go to write output files, while keeping the burden of memory usage low, by writing the actual
//! sequences to file (encoded as u8) and retreiving them only as needed
//! We need to use the Variants struct to add variants to the contigs, making this a sort of second
//! level struct
use std::io;
use thiserror::Error;
use std::num::ParseIntError;
use std::{collections::{HashMap, VecDeque}, io::Write};
use serde::{Deserialize, Serialize};
use serde;
use serde_json;
use std::fs::{read_to_string, File};
use std::path::Path;
use crate::structs::variants::Variant;
use log::{debug, error};

#[derive(Error, Debug)]
pub enum FastaMapError {
    #[error("FastaMap reported an IO Error: {0}")]
    IoError(#[from] io::Error),
    #[error("FastaMap reported bad coordinates.")]
    BadCoordinatesError,
    #[error("FastaMap reported an out of bounds error.")]
    OutOfBoundsError,
    #[error("FastaMap reported a ContigNotFoundError")]
    ContigNotFoundError(&'static str),
    #[error("FastaMap reported a serde error: {0}")]
    SerdeValueError(serde::de::value::Error),
    #[error("FastaMap reported an error from serde_json: {0}")]
    SerdeJsonError(serde_json::Error),
    #[error("FastaMap reported a malformed filename error: {0}")]
    MalformedFilenameError(ParseIntError),
}

impl From<serde_json::Error> for FastaMapError {
    fn from(error: serde_json::Error) -> Self {
        FastaMapError::SerdeJsonError(error)
    }
}

impl From<serde::de::value::Error> for FastaMapError {
    fn from(error: serde::de::value::Error) -> Self {
        FastaMapError::SerdeValueError(error)
    }
}

impl From<ParseIntError> for FastaMapError {
    fn from(error: ParseIntError) -> Self {
        FastaMapError::MalformedFilenameError(error)
    }
}

pub fn decode_file_name(path_string: &str) -> Result<(usize, usize), FastaMapError> {
    // This function parses the filename according to the established convention in NEAT
    // contigname_0001100000_0001200000.json
    // where the 2 groups of numbers are the start and end points of the sequence block
    // these coordinates are what we are trying to get from the name.
    let mut underscores: Vec<usize> = Vec::new();
    let mut index = 0;
    // We'll extract just the file stem part of the path and convert to string
    let filename = Path::new(path_string).file_stem().unwrap().display().to_string();
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

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct SequenceBlock {
    // This struct is for accessing the data stored in the sequence block file, a yaml file
    // The data file contains the following data:
    //   - contig: the key name from the contig list. This is the contig which this sequence is a part of
    //   - ref_start: the start point of this sequence relative to the full contig
    //   - ref_end: the end point of this sequence relative to the full contig
    //   - sequence: a list of u8 data that encodes the sequence
    pub contig: String,
    pub ref_start: usize,
    pub ref_end: usize,
    pub sequence: Vec<u8>,
}

impl SequenceBlock {
    pub fn from(mut filename: &str) -> Result<Self, FastaMapError> {
        // This is a check that the filename is properly encoded
        (_, _) = decode_file_name(&filename)?;

        let json_data = read_to_string(&mut filename)?;
        let sequence_block: SequenceBlock = serde_json::from_str(&json_data)?;
        Ok(sequence_block)
    }

    pub fn write_block(&mut self, mut file: File) -> std::io::Result<()> {
        // This will serialize the data to a byte vector and write it to the file
        let data = serde_json::to_vec(&self).unwrap();
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

    pub fn get_subseq(&self, request_start: usize, request_end: usize) -> Result<Vec<u8>, FastaMapError> {

        if request_start >= request_end {
            error!("Bad coordinates for sequence block retrieval - start: {} >= end: {}", request_start, request_end);
            return Err(FastaMapError::BadCoordinatesError)
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

#[derive(Clone, Debug)]
pub enum RegionType {
    NRegion,
    NonNRegion,
}

#[derive(Debug, Clone)]
// Each Contig contains a vector of filenames, each containing a sequence block in json format.
pub struct Contig {
    // This is the short name of the contig
    pub name: String,
    // The overall length of the contig in bases
    pub len: usize,
    // This is a vector of filenames. Each filename represents a sequence block and can be accessed via the contig block
    // The naming convention of the filenames will follow this format: contigname_0001100000_0001200000.jfa
    // Where contigname = the short name of the contig, the first string of digits is the start point of that file and the second block of digits is the end point 
    pub blocks: Vec<String>,
    // This is a map of the contig, giving the overall N-region versus NonN-regions of the contig
    pub contig_map: Vec<(usize, usize, RegionType)>,
}

impl Contig {
    pub fn new(
        name: String,
        len: usize,
        blocks: Vec<String>,
        map: Vec<(usize, usize, RegionType)>
    ) -> Result<Self, FastaMapError> {
        
        Ok(Contig {
            name,
            len,
            blocks,
            contig_map: map,
        })
    }

    pub fn get_sequence_block(&self, request_start: usize, request_end: usize) -> Result<Vec<String>, &'static str> {
        // This attempts to find a sequence block matching the first coordinate. Some coordinates may have two
        // files that have this location.
        // 
        // We may need a special case for when the location is in more than one block
        // Instead of a panic, we return a string error, in case the calling code is not necessarily expecting 
        // to find a result, it can deal with that itself.
        let mut found_blocks: Vec<String> = Vec::new();
        for block_name in &self.blocks {
            // Unwrapping here because we want a custom return error message
            let (start, end) = decode_file_name(block_name).unwrap();
            if ((request_start >= start) && (request_start < end)) && 
                ((request_end > start) && (request_end <= end)) {
                // In an attempt to keep features from spanning more than one sequence block, making
                // the simulaton more difficult, we're only looking for sequences that entirely contain
                // the feature. We could potentially improve on this idea in the future.
                //
                // block found! Note they may not be univque because of overlaps
                found_blocks.push(block_name.clone());
            }
        }
        if found_blocks.is_empty(){
            Err("No blocks found matching range")
        } else {
            Ok(found_blocks)
        }
    }

    pub fn map_variants(&self, variants_list: Vec<Variant>) -> Result<HashMap<Variant, Vec<String>>, FastaMapError> {
        // This returns a hashmap of the variants from the list as the key, with values being sequence blocks that the variant
        // intersects. This gives a complete map for the entire contig
        let mut variant_map: HashMap<Variant, Vec<String>> = HashMap::new();

        for variant in variants_list.iter() {
            let request_start = variant.location;
            // Guaranteed to work because a reference vector must have at least 1 element
            let request_end = request_start + variant.reference.len();
            let result = self.get_sequence_block(request_start, request_end);
            match result {
                Ok(result) => { variant_map.insert(variant.clone(), result); },
                // This is a little weird, but we'll debug later
                Err(e) => { debug!("{}", e) }
            }
        }
        Ok(variant_map)
    }
}

#[derive(Debug)]
// The FastaMap holds all the data for retrieving items from the fasta files
pub struct FastaMap {
    // A vector of Contig objects, one for each contig in the fasta file
    pub contigs: Vec<Contig>,
    // This maps the short contig name, used throughout, to the full name from the fasta file
    pub name_map: HashMap<String, String>,
    // This is a VecDeque with the contig short names in the order they appeared in the fasta
    // We need the correct order to output a valid file at the end.
    pub contig_order: VecDeque<String>,
}

impl FastaMap {
    pub fn new(
        contigs: Vec<Contig>, 
        name_map: HashMap<String, String>, 
        contig_order: VecDeque<String>
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
    use crate::structs::fasta_map::RegionType::{NRegion, NonNRegion};    use super::*;

    #[test]
    fn test_sequence_block() {
        let mut seq_block = SequenceBlock {
            contig: "chrom1".to_string(),
            ref_start: 0,
            ref_end: 20,
            sequence: vec![0,0,0,0,0,0,0,0,0,0,0,0,4,4,4,0,0,0,4,0],
        };
        let filename = "chrom1_000_020.json".to_string();
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
        let filename = "chr1_001000_002000.json".to_string();
        let (a, b) = decode_file_name(&filename).unwrap();
        assert_eq!(a, 1000);
        assert_eq!(b, 2000)
    }
    
    #[test]
    fn test_basic_map() {
        let filename = "chrom1_0000_0020.json".to_string();
        let sequences = vec![filename];
        let contig_map = Vec::from([
            (0, 12, NRegion),
            (12, 15, NonNRegion),
            (15, 18, NRegion),
            (18, 19, NonNRegion),
            (19, 20, NRegion),
        ]);
        let contig = Contig::new(
            "chrom1".to_string(),
            20,
            sequences,
            contig_map,
        ).unwrap();
        let name_map = HashMap::from([
            ("chrom1".to_string(), ">chrom1 foo bar\n".to_string())
        ]);
        let contigs = Vec::from([contig]);
        let contig_order = VecDeque::from(["chrom1".to_string()]);
        let fasta_map = FastaMap::new(
            contigs,
            name_map,
            contig_order,
        );

        let expected_output = "FastaMap { \
        contigs: [Contig { \
            name: \"chrom1\", \
            len: 20, \
            blocks: [\
                \"chrom1_0000_0020.json\"], \
            contig_map: [\
                (0, 12, NRegion), \
                (12, 15, NonNRegion), \
                (15, 18, NRegion), \
                (18, 19, NonNRegion), \
                (19, 20, NRegion)] }], \
            name_map: {\"chrom1\": \">chrom1 foo bar\\n\"}, \
            contig_order: [\"chrom1\"] }";

        let test_out = format!("{:?}", fasta_map);

        assert_eq!(test_out, expected_output)
    }
}