use std::collections::{HashMap, VecDeque};
use std::ffi::OsString;
use std::fs::{read_to_string, File};
use std::path::PathBuf;
use std::str::FromStr;

// The SequenceBlock has coordinates that map back to the original reference. Local manipulations
// of the block can use local coordinates, but any variant locations etc will be modified by the
// ref start point
#[derive(Debug, Clone)]
pub struct SequenceBlock {
    // start is the index for when this feature begins on the reference sequence. Start is always
    // relative to the original sequence
    pub ref_start: usize,
    // end is where this feature ends, relative to the reference.
    pub ref_end: usize,
    pub length: usize,
    pub file: OsString,
}

impl SequenceBlock {
    pub fn new(ref_start: usize, length:usize, file: OsString) -> Result<Self, ()> {
        if length <= 0 {
            panic!(
                "Bad parameters for sequence block start: {} length: {}", ref_start, length
            )
        } else {
            Ok(SequenceBlock {
                ref_start,
                ref_end: ref_start + length,
                length,
                file,
            })
        }
    }

    pub fn len(&self) -> usize {
        self.length
    }

    pub fn retrieve_all(&self) -> Result<Vec<u8>, ()> {
        // There should only be one line, and since these are programmatically generated, there
        // shouldn't be an issue with security. If the format is wrong, then the map part should
        // throw an error. Had to use vec::from rather than collect for reasons unknown (can't
        // clone a u8?
        let mut result_vec = Vec::new();
        for line in read_to_string(&self.file).unwrap().lines() {
            for char in line.chars() {
                result_vec.push(char as u8);
            }
        };
        Ok(result_vec)
    }

    pub fn retrieve(&self, start: usize, end: usize) -> Result<Vec<u8>, &'static str> {
        if start >= end {
            panic!(
                "Bad coordinates for sequence block retrieval: {} >= end: {}", start, end
            );
        } else if end > self.ref_end || start < self.ref_start {
            panic!(
                "Bad coordinates in sequence block retrieval -\
                 start: {}, end: {}, block start: {}, block end: {}",
                start, end, self.ref_start, self.ref_end,
            );
        } else {
            // There should only be one line, and since these are programmatically generated, there
            // shouldn't be an issue with security. If the format is wrong, then the map part should
            // throw an error. Had to use vec::from rather than collect for reasons unknown (can't
            // clone a u8?
            let buffer: Vec<u8> = self.retrieve_all().unwrap();
            // start for this function is expected to be a coordinate from the overall contig,
            // so we modify by the known start and end point of this block to get the final
            // coordinates
            let block_start = start - self.ref_start;
            let block_end = end - self.ref_start;
            Ok(buffer[block_start..block_end].to_owned())
        }
    }
}

#[derive(Clone, Debug)]
pub enum RegionType {
    NRegion,
    NonNRegion,
}

#[derive(Debug)]
// Each Contig contains multiple sequence blocks, which map to files containing sequences that
// are 128Kb long.
pub struct Contig {
    pub name: String,
    pub len: usize,
    pub blocks: Vec<SequenceBlock>,
    pub contig_map: Vec<(usize, usize, RegionType)>
}

impl Contig {
    pub fn new(
        name: String,
        len: usize,
        blocks: Vec<SequenceBlock>,
        map: Vec<(usize, usize, RegionType)>
    ) -> Result<Self, &'static str> {

        Ok(Contig {
            name,
            len,
            blocks,
            contig_map: map,
        })
    }
}

#[derive(Debug)]
// The fasta map holds all the data for retrieving items from the fasta files
pub struct FastaMap {
    // A vector of Contig objects, one for each contig in the fasta file
    pub contigs: Vec<Contig>,
    // This maps the short contig name, used throughout, to the full name from the fasta file
    pub name_map: HashMap<String, String>,
    // This is a VecDeque with the contig short names in the order they appeared in the fasta
    // We need the correct order to output a valid file at the end.
    pub contig_order: VecDeque<String>,
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_map() {
        // use crate::structs::nucleotides::Nuc::*;
        //
        // let read = vec![
        //     Unk, Unk, Unk, Unk, Unk, Unk, Unk, Unk, Unk, Unk,
        //     Unk, Unk, Ada, Cyt, Thy, Unk, Unk, Unk, Ada, Unk
        // ];
        // let fasta_sequence = Box::new(HashMap::from([
        //     ("chr 1".to_string(), read)
        // ]));
        // let fasta_map = FastaMap::build_from_reference(&fasta_sequence).unwrap();
        //
        // let expected_output = "FastaMap { map: {\"chr 1\": \
        // (20, [\
        // FastaBlock { start: 0, block_type: NBlock }, \
        // FastaBlock { start: 12, block_type: Allowed }, \
        // FastaBlock { start: 15, block_type: NBlock }, \
        // FastaBlock { start: 18, block_type: Allowed }, \
        // FastaBlock { start: 19, block_type: NBlock }\
        // ])} }".to_string();
        //
        // let test_out = format!("{:?}", fasta_map);
        //
        // assert_eq!(test_out, expected_output)
    }
}