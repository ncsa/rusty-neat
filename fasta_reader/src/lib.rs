use std::collections::{HashMap, VecDeque};
use std::io;
use std::hash::Hash;
use common::file_tools::read_lines;

pub const BUFSIZE: usize = 131_072usize; // i.e., 128kb
pub const DICT_SIZE: usize = 32768usize; // i.e., 32kb

pub struct FastaBlock {

}

pub fn read_fasta(
    fasta_path: &str,
    fragment_max: usize,
) -> Result<(Box<HashMap<String, HashMap<u64, Vec<u8>>>>, VecDeque<String>), io::Error> {
    // Reads a fasta file and turns it into a HashMap and puts it in the heap
    let mut fasta_map: HashMap<String, HashMap<u64, Vec<u8>>> = HashMap::new();
    let mut fasta_order: VecDeque<String> = VecDeque::new();
    let mut current_key = String::new();
    let mut block_map: HashMap<u64, Vec<u8>>;
    let lines = read_lines(fasta_path)?;
    let mut temp_seq = String::new();
    // Index will be used to create hashmaps of readable areas and non-readable areas
    let mut index: u64 = 0;
    // For the purposes of NEAT, everything that isn't ACTG or actg is treated like an unknown
    // base (N). NEAT does not yet have the capability of sorting out some of the alternate bases
    // like R indicating a purine (A or G).

    // If we encounter N's or unknown bases, this will count them up for us.
    let mut n_run: u64 = 0;
    let mut within_n_block = false;
    let mut buffer = [4_u8; DICT_SIZE];
    let mut next_buffer = [4_u8; DICT_SIZE];
    let buffer_index = 0;
    for line in lines {
        match line {
            Ok(l) => {
                if l.starts_with(">") {
                    let (key, name) = extract_key_name(&l[1..]);
                }
            },
            _ => { panic!("AAAAAAAA"); },
        }
    }
    Ok((Box::new(fasta_map), fasta_order))
}

fn extract_key_name(key: &str) -> (String, String) {
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
    (key[..delim_index].to_string(), key.to_string())

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
        assert_eq!(extract_key_name(name1), ("BBB".to_string(), "BBB|AAA|CCC".to_string()));
        assert_eq!(extract_key_name(name2), ("BBB".to_string(), "BBB AAA CCC".to_string()));
        assert_eq!(extract_key_name(name3), ("BBBAAACCC".to_string(), "BBBAAACCC".to_string()));
    }
}
