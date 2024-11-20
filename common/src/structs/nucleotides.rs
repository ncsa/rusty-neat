// Representing bases by a simple u8 num
//     A = 0
//     C = 1
//     G = 2
//     T = 3
//     N (or other unknown chars) = 4
// This is intended to make it easier to store them. I think this will be an easier way to read
// the code, and I'm hoping it does not incur any extra computational burden

pub fn base_to_string(b: u8) -> String {
    base_to_char(b).to_string()
}

pub fn base_to_char(b: u8) -> char {
    match b {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        // Todo: replace the "N" with a randomly generated character from ACGT
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


pub fn char_to_base(c: char) -> u8 {
    // This defines the relationship between the 4 possible nucleotides in DNA and
    // the character from an input sequence. Everything that isn't a recognized base is
    // considered an N for Neat. Note that NEAT ignores soft masking.
    match c {
        'A' | 'a' => 0,
        'C' | 'c' => 1,
        'G' | 'g' => 2,
        'T' | 't' => 3,
        _ => 4,
    }
}

pub fn sequence_array_to_string(input_array: Vec<u8>) -> String {
    // Converts a sequence vector into a string representing the DNA sequence
    let mut return_string = String::with_capacity(input_array.len());
    for nuc in input_array {
        return_string.push(base_to_char(nuc));
    }
    return_string
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_to_nuc() {
        let test_nuc = 'A';
        assert_eq!(char_to_base(test_nuc), 0);
    }

    #[test]
    fn cast_type() {
        let num1: u8 = 0;
        assert_eq!(base_to_char(num1), 'A')
    }

    #[test]
    fn test_nuc_to_char() {
        let my_nuc = 1;
        assert_eq!(base_to_string(my_nuc), "C");
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
