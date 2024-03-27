use utils::nucleotides::Nuc;
#[derive(Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub enum VariantType {
    Indel,
    SNP,
}
pub struct Variant {
    // This is a basic holder for a variant. There are several aspects of the variant that we
    // don't need to store, such as which ploid the variant appears on and the context sequence,
    // which we only need for the output files, we can generate at the time we need it.
    //
    // chromosome: uniquely identifies the contig from the fasta file.
    // position: the position within the reference (0-indexed) where the first ref base is
    //      located.
    // reference: the vector of nucelotides on chromosome starting at position which are altered
    //      for this variant
    // alternate: the vector of nucleotides that take the place of reference.
    chromosome: String,
    position: usize,
    reference: Vec<Nuc>,
    alternate: Vec<Nuc>,
}