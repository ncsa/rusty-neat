// Throughout this program, we are standardizing the use of a u8 representation of the nucleotides
//     A = 0
//     C = 1
//     G = 2
//     T = 3
//     N = 4
// This is intended to make it easier to store them. We thought about using the u8 representation
// of the character as built into Rust, but we'd then have to figure out the translations and keep
// track of extra numbers. So this is intended to simplify everything
enum Nuc {
    A,
    C,
    G,
    T,
}
impl Nuc {
    pub fn base_to_nuc(char_of_interest: char) -> u8 {
        // This defines the relationship between the 4 possible nucleotides in DNA and
        // a simple u8 numbering system. Everything that isn't a recognized base is a 4.
        // Note that NEAT ignores soft masking.
        return match char_of_interest {
            'A' | 'a' => Nuc::A,
            'C' | 'c' => 1,
            'G' | 'g' => 2,
            'T' | 't' => 3,
            _ => 4,
        }
    }

    pub fn u8_to_base(nuc_num: u8) -> char {
        // Canonical conversion from base u8 representation back into the character.
        // We're returning a string instead of a char to facilitate. No attempt to preserve or display
        // any soft masking.
        return match nuc_num {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'N',
        }
    }
}