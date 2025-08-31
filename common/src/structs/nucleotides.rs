//! Representing bases by a simple u8 num
//!     A = 0
//!     C = 1
//!     G = 2
//!     T = 3
//!     N (or other unknown chars) = 4
//! This is intended to make it easier to store them. Note that the bases are always in alphabetical order
//! with N tacked on at the end.
const ALLOWED_NUCS: [char; 4] = ['A', 'C', 'G', 'T'];

pub fn decode_base(b: u8) -> char {
    match b {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => 'N',
    }
}

pub fn complement(base: u8) -> u8 {
    // matches with the complement of each nucleotide.
    match base {
        0 => 3,
        1 => 2,
        2 => 1,
        3 => 0,
        _ => base,
    }
}


pub fn encode_base(c: char) -> u8 {
    // This defines the relationship between the 4 possible nucleotides in DNA and
    // the character from an input sequence. Everything that isn't a recognized base is
    // considered an N for Neat. Note that NEAT ignores soft masking.
    match c {
        'A' | 'a' => 0,
        'C' | 'c' => 1,
        'G' | 'g' => 2,
        'T' | 't' => 3,
        'U' | 'u' => 5, // Note this is prep work, not implemented
        _ => 4,
    }
}

#[allow(unused)]
pub fn encode_aa(c: char) -> u8 {
    // The purpose of this function is to encode standard amino acids as a u8
    match c {
        'A' => 255,
        'R' => 254,
        'N' => 253,
        'D' => 252,
        'C' => 251,
        'E' => 250,
        'Q' => 249,
        'G' => 248,
        'H' => 247,
        'I' => 246,
        'L' => 245,
        'K' => 244,
        'M' => 243,
        'F' => 242,
        'P' => 241,
        'S' => 240,
        'T' => 239,
        'W' => 238,
        'Y' => 237,
        'V' => 236,
        'X' => 235, // stop codon
        _ => panic!("Attempted to encode unknown amino acid: {c}")
    }
}

#[allow(unused)]
pub fn interpret_aa(name: &str) -> char {
    // Translates string representations of the amino acids into char
    match name.to_lowercase().as_str() {
        "alanine"       | "ala" => 'A',
        "arginine"      | "arg" => 'R',
        "asparagine"    | "asn" => 'N',
        "aspartic acid" | "asp" => 'D',
        "cysteine"      | "cys" => 'C',
        "glutamic acid" | "glu" => 'E',
        "glutamine"     | "gln" => 'Q',
        "glycine"       | "gly" => 'G',
        "histidine"     | "his" => 'H',        
        "isoleucine"    | "ile" => 'I',
        "leucine"       | "leu" => 'L',
        "lysine"        | "lys" => 'K',
        "methionine"    | "met" => 'M',
        "phenylaline"   | "phe" => 'F',
        "proline"       | "pro" => 'P',
        "serine"        | "ser" => 'S',
        "threonine"     | "thr" => 'T',
        "tryptophan"    | "trp" => 'W',
        "tyrosine"      | "tyr" => 'Y',
        "valine"        | "val" => 'V',
        _ => panic!("Unknown amino acid: {}", name),
    }
}

pub fn sequence_array_to_string(input_array: &Vec<u8>) -> String {
    // Converts a sequence vector into a string representing the DNA sequence
    let mut return_string = String::with_capacity(input_array.len());
    for nuc in input_array {
        return_string.push(decode_base(*nuc));
    }
    return_string
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_to_nuc() {
        let test_nuc = 'A';
        assert_eq!(encode_base(test_nuc), 0);
    }

    #[test]
    fn cast_type() {
        let num1: u8 = 0;
        assert_eq!(decode_base(num1), 'A')
    }

    #[test]
    fn test_complement() {
        let nuc1 = 0_u8;
        let nuc2 = 1_u8;
        let nuc3 = 2_u8;
        let nuc4 = 3_u8;
        let nuc5 = 4_u8;

        assert_eq!(complement(nuc1), nuc4);
        assert_eq!(complement(nuc2), nuc3);
        assert_eq!(complement(nuc3), nuc2);
        assert_eq!(complement(nuc4), nuc1);
        assert_eq!(complement(nuc5), nuc5);
    }
}
