// Throughout this program, we will build off the Nuc enum. Enums have implied numbering, so that:
//     A = 0
//     C = 1
//     G = 2
//     T = 3
//     N = 4
// This is intended to make it easier to store them. I think this will be an easier way to read
// the code, and I'm hoping it does not incur any extra computational burden
#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub enum Nuc {
    A,
    C,
    G,
    T,
    N,
}
impl Nuc {
    pub fn int_to_nuc(input_integer: usize) -> Nuc {
        // Similar to previous, but allows an integer representation, useful for sampling an index.
        return match input_integer {
            0 => Nuc::A,
            1 => Nuc::C,
            2 => Nuc::G,
            3 => Nuc::T,
            _ => Nuc::N,
        };
    }

    pub fn to_char(&self) -> char {
        // Canonical conversion from base u8 representation back into the character.
        // We're returning a string instead of a char to facilitate. No attempt to preserve or display
        // any soft masking.
        match self {
            Nuc::A => 'A',
            Nuc::C => 'C',
            Nuc::G => 'G',
            Nuc::T => 'T',
            Nuc::N => 'N',
        }
    }

    pub fn to_str(&self) -> String {
        self.to_char().to_string()
    }

    pub fn complement(&self) -> Nuc {
        // matches with the complement of each nucleotide.
        return match self {
            Nuc::A => Nuc::T,
            Nuc::C => Nuc::G,
            Nuc::G => Nuc::C,
            Nuc::T => Nuc::A,
            Nuc::N => Nuc::N,
        };
    }
}

pub fn base_to_nuc(char_of_interest: char) -> Nuc {
    // This defines the relationship between the 4 possible nucleotides in DNA and
    // the character from an input sequence. Everything that isn't a recognized base is
    // considered an N for Neat. Note that NEAT ignores soft masking.
    return match char_of_interest {
        'A' | 'a' => Nuc::A,
        'C' | 'c' => Nuc::C,
        'G' | 'g' => Nuc::G,
        'T' | 't' => Nuc::T,
        _ => Nuc::N,
    };
}

pub fn sequence_array_to_string(input_array: &[Nuc]) -> String {
    // Converts a sequence vector into a string representing the DNA sequence
    let mut return_string = String::new();
    for nuc in input_array {
        return_string += &nuc.to_str();
    }
    return_string
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_to_nuc() {
        let test_nuc = 'A';
        assert_eq!(base_to_nuc(test_nuc), Nuc::A);
    }

    #[test]
    fn cast_type() {
        let num1: u8 = 0;
        assert_eq!(Nuc::int_to_nuc(num1 as usize), Nuc::A)
    }

    #[test]
    fn test_num_to_nuc() {
        let num2: i8 = 1;
        let num3: i32 = 2;
        let num4: usize = 3;
        let num5: u64 = 9;
        assert_eq!(Nuc::int_to_nuc(num2 as usize), Nuc::C);
        assert_eq!(Nuc::int_to_nuc(num3 as usize), Nuc::G);
        assert_eq!(Nuc::int_to_nuc(num4), Nuc::T);
        assert_eq!(Nuc::int_to_nuc(num5 as usize), Nuc::N)
    }

    #[test]
    fn test_nuc_to_char() {
        let my_nuc = Nuc::C;
        assert_eq!(my_nuc.to_char(), 'C');
    }
    #[test]
    fn test_complement() {
        let nuc1 = Nuc::A;
        let nuc2 = Nuc::C;
        let nuc3 = Nuc::G;
        let nuc4 = Nuc::T;
        let nuc5 = Nuc::N;

        assert_eq!(nuc1.complement(), Nuc::T);
        assert_eq!(nuc2.complement(), Nuc::G);
        assert_eq!(nuc3.complement(), Nuc::C);
        assert_eq!(nuc4.complement(), Nuc::A);
        assert_eq!(nuc5.complement(), Nuc::N);
    }
}
