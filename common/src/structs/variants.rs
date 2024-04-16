use structs::nucleotides::Nuc;
use self::VariantType::*;

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub enum VariantType {
    SNP,
    Insertion,
    Deletion,
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
    pub reference: Vec<Nuc>,
    // The alternate allele of interest. This is either one base or several bases for an insertion.
    pub alternate: Vec<Nuc>,
    // the genotype will tell us which ploid is mutated in vcf, fastq and bam.
    pub genotype: Vec<u8>,
}

impl Variant {
    pub fn new(
        variant_type: VariantType,
        reference: &Vec<Nuc>,
        alternate: &Vec<Nuc>,
        genotype: Vec<u8>,
    ) -> Self {
        match variant_type {
            Insertion | Deletion => assert_ne!(reference.len(), alternate.len()),
            SNP => assert_eq!(reference.len(), alternate.len()),
        }

        Variant {
            variant_type,
            reference: reference.clone(),
            alternate: alternate.clone(),
            genotype: genotype.clone(),
        }
    }

    pub fn is_homozygous(&self) -> bool {
        return if self.genotype.contains(&0) { false } else { true }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use structs::nucleotides::Nuc::*;
    use structs::variants::VariantType::*;

    #[test]
    fn test_variant_creation() {
        let variant = Variant {
            // Not actually a valid variant, but just testing the constructor
            variant_type: SNP,
            reference: vec![A, C, T, G],
            alternate: vec![A, P(3)],
            genotype: vec![1, 1, 1],
        };
        assert_eq!(variant.variant_type, SNP);
        assert_eq!(variant.variant_type, SNP);
        assert!(variant.is_homozygous());
    }

    #[test]
    fn test_variant_new() {
        let variant = Variant::new(
            Insertion,
            &vec![A, P(3)],
            &vec![A, C, T, G],
            vec![0, 1],
        );
        assert!(!variant.is_homozygous())
    }

    #[test]
    fn test_deletion_new() {
        let variant = Variant::new(
            Deletion,
            &vec![A, C, T, G],
            &vec![A, P(3)],
            vec![0, 1],
        );
        assert!(!variant.is_homozygous())
    }

    #[test]
    #[should_panic]
    fn test_bad_variant_creation() {
        Variant::new(
            // with new it should catch this error
            SNP,
            &vec![A, C, T, G],
            &vec![A, P(3)],
            vec![1, 1, 1],
        );
    }
}
