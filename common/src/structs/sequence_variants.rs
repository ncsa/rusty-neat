//! This higher level struct will combine variant data and struct data, and will keep track
//! of which sequence blocks are affected
use crate::structs::{
    fasta_map::FastaMapError,
    variants::{Variant, VariantError}
};
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SeqBlockVariantError {
    #[error("ContigVariant reported a FastaMap error: {0}")]
    SeqBlockError(FastaMapError),
    #[error("ContigVariant reported a Variant error: {0}")]
    VariantError(VariantError),
    #[error("ContiVariant reported a mapping error")]
    MappingError,
}

#[derive(Debug)]
pub struct SequenceBlockVariants {
    // This struct will keep track of all variants on a contig.
    //
    // This is the filename for the seq_block
    seq_block: PathBuf,
    // The list of variants on the contig
    variants_list: Vec<Variant>,
}

impl SequenceBlockVariants {
    pub fn new(
        seq_block: PathBuf,
        variants: Vec<Variant>,
    ) -> Result<Self, SequenceBlockVariants> {
        // This creates a basic contig variants box without the mapping of variants to sequence blocks
        Ok(SequenceBlockVariants { seq_block, variants_list: variants })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::nucleotides::Nucleotide;
    use crate::structs::variants::{Variant, VariantType};
    use std::path::PathBuf;

    fn setup() -> Vec<Variant> {
        // create a variant
        let variant_type = VariantType::SNP;
        let location = 55;
        let reference: Vec<Nucleotide> = vec![Nucleotide::G];
        let alternate: Vec<Nucleotide> = vec![Nucleotide::T];
        let mut genotype: Vec<usize> = vec![1,0];
        let variant = Variant::new(variant_type, location, &reference, &alternate, &mut genotype).unwrap();
        let variants_list = vec![variant];
        variants_list
    }

    #[test]
    fn test_new_unmapped() {
        let vlist = setup();
        let fake_block = PathBuf::from("chr1.json");
        let sbv = SequenceBlockVariants::new(fake_block, vlist.to_owned()).unwrap();
        assert_eq!(sbv.seq_block, PathBuf::from("chr1.json"));
        assert_eq!(sbv.variants_list, vlist);
    }

}