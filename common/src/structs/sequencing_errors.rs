use crate::structs::nucleotides::Nuc;

#[derive(PartialEq, Debug)]
pub struct SnpErr {
    nucleotide: Nuc,
}

impl SnpErr {
    pub fn new(nucleotide: Nuc) -> Self {
        SnpErr { nucleotide }
    }

    pub fn retrieve(&self) -> Nuc {
        self.nucleotide
    }
}

#[derive(Clone)]
pub struct IndelErr {
    // Both insertions and deletions have a length (positive for insertions, negative for deletions)
    length: i8,
    // Only insertions need to return the vector of added nucleotides.
    sequence: Option<Vec<Nuc>>,
}

impl IndelErr {
    pub fn new(length: i8, sequence: Option<Vec<Nuc>>) -> Self {
        IndelErr { length, sequence }
    }

    pub fn length(&self) -> i8 {
        self.length
    }

    pub fn sequence(&self) -> &Option<Vec<Nuc>> {
        &self.sequence
    }
}

pub enum SequencingError {
    // The only thing we store here are types of errors. Sequencing errors are generated on the fly,
    // and they are not stored for very long, so we only need a type to generate.
    SNPError(SnpErr),
    IndelError(IndelErr),
}

pub enum SequencingErrorType {
    SNPError,
    IndelError,
}

#[cfg(test)]
mod tests {
    use super::*;
    use structs::nucleotides::Nuc::*;

    #[test]
    fn test_snp() {
        let nucleotide = A;
        let snp = SnpErr { nucleotide };
        assert_eq!(snp, SnpErr::new(nucleotide));
        assert_eq!(snp.retrieve(), nucleotide);
    }

    #[test]
    fn test_insertion() {
        let length = 16;
        let sequence = vec![A, A, A, A, A, A, A, A, G, T, A, C, G, T, A, A];
        let insertion = IndelErr {
            length,
            sequence: Some(sequence.clone()),
        };
        assert_eq!(insertion.length(), 16);
        assert_eq!(insertion.sequence.unwrap(), sequence);
    }

    #[test]
    fn test_deletion() {
        let length = -16;
        let deletion = IndelErr {
            length,
            sequence: None,
        };
        assert_eq!(deletion.length(), -16);
        assert!(deletion.sequence.is_none());
    }
}
