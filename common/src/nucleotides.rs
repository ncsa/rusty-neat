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
    pub fn base_to_nuc(char_of_interest: char) -> Nuc {
        // This defines the relationship between the 4 possible nucleotides in DNA and
        // a simple u8 numbering system. Everything that isn't a recognized base is a 4.
        // Note that NEAT ignores soft masking.
        return match char_of_interest {
            'A' | 'a' => Nuc::A,
            'C' | 'c' => Nuc::C,
            'G' | 'g' => Nuc::G,
            'T' | 't' => Nuc::T,
            _ => Nuc::N,
        }
    }

    pub fn nuc_to_char(&self) -> char {
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

    pub fn complement(&self) -> Nuc {
        // matches with the complement of each nucleotide.
        return match self {
            Nuc::A => Nuc::T,
            Nuc::C => Nuc::G,
            Nuc::G => Nuc::C,
            Nuc::T => Nuc::A,
            Nuc::N => Nuc::N,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_to_nuc() {
        let test_nuc = 'A';
        assert_eq!(Nuc::base_to_nuc(test_nuc), Nuc::A);
    }
    #[test]
    fn test_nuc_to_char() {
        let my_nuc = Nuc::C;
        assert_eq!(my_nuc.nuc_to_char(), 'C');
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