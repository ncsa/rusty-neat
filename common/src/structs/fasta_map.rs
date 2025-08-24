//! This library contains low level structs designed to catalog and index the fasta file. 
//! Fasta files can be quite large and so this system is meant to speed up retreival when we
//! go to write output files, while keeping the burden of memory usage low, by writing the actual
//! sequences to file (encoded as u8) and retreiving them only as needed

use std::{collections::{HashMap, VecDeque}, io::{Write, Read}};
use serde::{Deserialize, Serialize};
use serde::de::value::Error;
use serde_json;
use std::fs::{read_to_string, File};

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
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

impl From<String> for SequenceBlock {
    fn from(mut filename: String) -> Self {
        let json_data = read_to_string(&mut filename).unwrap();
        let sequence_block: SequenceBlock = serde_json::from_str(&json_data).unwrap();
        sequence_block
    }
}

impl SequenceBlock {

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

    pub fn get_subseq(&self, request_start: usize, request_end: usize) -> Result<Vec<u8>, Box<Error>> {

        if request_start >= request_end {
            panic!(
                "Bad coordinates for sequence block retrieval - start: {} >= end: {}", request_start, request_end
            );
        } else if (request_end <= self.ref_start) || 
                  (request_start > self.ref_end) || 
                  (request_end > self.ref_end) || 
                  (request_start < self.ref_start) {
            panic!(
                "Bad requested coordinates in sequence block retrieval -\
                 start: {}, end: {}, block start: {}, block end: {}",
                request_start, request_end, self.ref_start, self.ref_end,
            );
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
    pub contig_map: Vec<(usize, usize, RegionType)>
}

impl Contig {
    pub fn new(
        name: String,
        len: usize,
        blocks: Vec<String>,
        map: Vec<(usize, usize, RegionType)>
    ) -> Result<Self, &'static str> {

        Ok(Contig {
            name,
            len,
            blocks,
            contig_map: map,
        })
    }

    pub fn get_sequence_block(&self, location: usize) -> Result<SequenceBlock, &'static str> {
        let mut index = 0;
        let mut underscores = Vec::new();
        for block_name in &self.blocks {
            for char in block_name.chars() {
                // Find all the underscores
                match char {
                    '_' => {
                        // Find the underscores
                        underscores.push(index);
                        index += 1; 
                    },
                    '.' => {
                        // a period will denote the filename extension for our purposes
                        break
                    },
                    _ => index += 1,
                }
            }
            // the last two blocks: _000111000_000112000 contain the range
            let penultimate_underscore = underscores[underscores.len()-2];
            let ultimate_underscore = underscores[underscores.len()-2];
            // This should parse out to a usize position
            let start: usize = block_name[(penultimate_underscore+1)..ultimate_underscore].parse().unwrap();
            let end: usize = block_name[(ultimate_underscore+1)..index].parse().unwrap();
            if (location >= start) && (location < end) {
                // block found!
                return Ok(SequenceBlock::from(block_name.to_owned()));
            }
        }
        Err("Block not found")
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

impl From<(Vec<Contig>, HashMap<String, String>, VecDeque<String>)> for FastaMap {
    fn from(
        (contigs, name_map, contig_order): 
        (Vec<Contig>, HashMap<String, String>, VecDeque<String>)
    ) -> Self {
        Self { contigs, name_map, contig_order }
    }
}

impl FastaMap {
    pub fn retrieve_contig(&self, request_name: String) -> Result<Contig, &'static str> {
        for contig in self.contigs.iter() {
            if contig.name == request_name {
                return Ok(contig.clone())
            }
        }
        Err("Contig not found")
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
        let seq_block_read = SequenceBlock::from(filename.clone());
        assert_eq!(seq_block.contig, seq_block_read.contig);
        assert_eq!(seq_block.ref_start, seq_block_read.ref_start);
        assert_eq!(seq_block.ref_end, seq_block_read.ref_end);
        assert_eq!(seq_block.sequence, seq_block_read.sequence);
        remove_file(filename).unwrap();
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
        let fasta_map = FastaMap::from((
            contigs,
            name_map,
            contig_order,
        ));

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