// This is a basic mutation with SNPs using a basic mutation model.
// mutate_fasta takes a fasta Hashmap and returns a mutated version and the locations of the
// mutations introduced
//
// mutate_sequence adds actual mutations to the fasta sequence

use std::collections::HashMap;
use common::structs::variants::{Variant, VariantType};

pub fn apply_mutations(
    sequence: &Vec<u8>,
    mutations_to_apply: &HashMap<usize, Variant>,
) -> Result<Vec<u8>, Err> {
    // Takes:
    // sequence: a vector of nucleotides representing the reference sequence.
    // mutations_to_apply: A hashmap of mutations to apply indexed by position.
    //
    // Returns:
    // A vector of nucleotides representing the mutated sequence
    //
    let mut mutated_sequence: Vec<u8> = Vec::new();
    let mut i = 0;
    while i < sequence.len() {
        if mutations_to_apply.contains_key(&i) {
            let variant_to_apply = mutations_to_apply.get(&i).unwrap();
            match variant_to_apply.variant_type {
                VariantType::Insertion => {
                    for base in &variant_to_apply.alternate {
                        mutated_sequence.push(base.clone());
                    }
                    i += 1;
                },
                VariantType::Deletion => {
                    mutated_sequence.push(sequence[i].clone());
                    i += variant_to_apply.reference.len();
                },
                VariantType::SNP => {
                    mutated_sequence.push(variant_to_apply.alternate[0].clone());
                    i += 1;
                },
            }
        } else {
            mutated_sequence.push(sequence[i].clone());
            i += 1;
        }
    }
    Ok(mutated_sequence)
}

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;
    use super::*;

    #[test]
    fn test_apply_mutations() {
        let seq1: Vec<u8> = vec![0, 0, 1, 3, 3, 2, 2, 1, 0, 3, 1, 1, 1];
        let mutations = HashMap::from([
            (3, Variant::new(
                VariantType::SNP,
                3,
                &vec![3],
                &vec![2],
                &vec![0, 1]
            )),
            (5, Variant::new(
                VariantType::Deletion,
                5,
                &vec![2, 2, 1],
                &vec![2],
                &vec![1, 1]
            )),
        ]);
        let mutant = apply_mutations(&seq1, &mutations).unwrap();
        assert_eq!(mutant.len(), seq1.len() - 2);
        assert_eq!(mutant[3], 2);
        assert_eq!(mutant[6], 0);
    }

    #[test]
    fn test_output_correct() {
        let mut fasta_map = HashMap::from([
            ("chr1".to_string(), vec![0, 1, 1, 2, 3, 2, 2, 1, 0, 2, 0]),
            ("chr2".to_string(), vec![3, 3, 0, 2, 1, 1, 3, 3, 0, 2, 2, 3, 3, 3, 0])
        ]);
        let chrom_order = VecDeque::from(
            ["chr1".to_string(), "chr2".to_string()]
        );
        let all_variants = HashMap::from([
            ("chr1".to_string(), HashMap::from([
                (3, Variant::new(
                    VariantType::SNP,
                    3,
                    &vec![2],
                    &vec![0],
                    &vec![0, 1]
                )),
                (5, Variant::new(
                    VariantType::Deletion, 5,
                    &vec![2, 2, 1],
                    &vec![2],
                    &vec![1, 1]
                ))
            ])),
            ("chr2".to_string(), HashMap::from([
                 (7, Variant::new(
                     VariantType::SNP,
                     7,
                     &vec![3],
                     &vec![2],
                     &vec![0, 1]
                 )),
                 (11, Variant::new(
                     VariantType::Insertion,
                     11,
                     &vec![3],
                     &vec![3, 3],
                     &vec![1, 1]
                 ))
            ])),
        ]);
        let output_map = HashMap::from([
            ("chr1", vec![0, 1, 1, 0, 3, 2, 0, 2, 0]),
            ("chr2", vec![3, 3, 0, 2, 1, 1, 3, 2, 0, 2, 2, 3, 3, 3, 3, 0])
        ]);
        for contig in chrom_order {
            let contig_variants = all_variants.get(&contig).unwrap().clone();
            let sequence_to_mutate = &fasta_map.get(&contig).unwrap().clone();
            let sequence_to_mutate = apply_mutations(
                &sequence_to_mutate,
                &contig_variants,
            ).expect("Error applyng mutations");

            fasta_map.insert(contig.clone(), sequence_to_mutate);
        }

        assert_eq!(fasta_map["chr1"], output_map["chr1"]);
        assert_eq!(fasta_map["chr2"], output_map["chr2"]);
    }
}