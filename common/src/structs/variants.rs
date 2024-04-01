use crate::structs::nucleotides::Nuc;

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub enum VariantType {
    SNP,
    Indel,
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Variant {
    // This is a basic holder for a variant. There are several aspects of the variant that we
    // don't need to store, such as which ploid the variant appears on and the context sequence,
    // which we only need for the output files, we can generate at the time we need it, such as the
    // insertion sequence.
    //
    // This is the type of variant. Any new types will need to be added to VariantType
    // to be implemented here.
    pub variant_type: VariantType,
    // The reference allele of interest. This is either one base or several bases for a deletion.
    pub reference: &'static Vec<Nuc>,
    // The alternate allele of interest. This is either one base or several bases for an insertion.
    pub alternate: &'static Vec<Nuc>,
    // the genotype will tell us which ploid is mutated in vcf, fastq and bam.
    pub genotype: Vec<u8>,
    // Saved later for convenience
    pub is_homozygous: bool,
}

impl Variant {
    pub fn new(
        variant_type: VariantType,
        reference: &Vec<Nuc>,
        alternate: &Vec<Nuc>,
        genotype: Vec<u8>,
        is_homozygous: bool,
    ) -> Self {
        Variant {
            variant_type,
            reference,
            alternate,
            genotype,
            is_homozygous,
        }
    }

    pub fn is_insertion(&self) -> bool {
        self.reference.len() < self.alternate.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use structs::nucleotides::Nuc::*;

    #[test]
    fn test_get_len() {
        let my_variant = Variant::new(
            VariantType::Indel,
            &vec![A, A, A, T, C, G, T, T, T, A],
            &vec![T, G, A, A, T, G],
            vec![0, 1],
            true,
        );
        assert_eq!(my_variant.get_length(), 5)
    }

    #[test]
    fn test_contains() {
        let test_variant = Variant::new(
            VariantType::Indel,
            &vec![A, A, A, T, C, G, T, T, T, A],
            &vec![T, G, A, A, T, G],
            vec![0, 1],
            false,
        );
        let test_pos = 13;
        assert!(test_variant.contains(test_pos));
    }
}
