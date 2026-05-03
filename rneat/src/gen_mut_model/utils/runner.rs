use std::{collections::HashMap, path::PathBuf};
use itertools::Itertools;
use log::*;

use common::{
    models::{mutation_model::MutationModel, snp_trinuc_model::{TrinucFrame}}, 
    structs::{
        bed_record::BedRecord, 
        fasta_map::{FastaMap}, 
        nucleotides::{Nucleotide, ALLOWED_NUCS}, 
        variants::{Genotype, Variant, VariantType}
    }
};

use crate::gen_mut_model::errors::GenMutationModelError;

pub fn runner(
    fasta_map: Box<FastaMap>,
    filtered_mutations: HashMap<String, Vec<Variant>>,
    bed_table: HashMap<String, Vec<BedRecord>>,
    output_file: &PathBuf,
) -> Result<(), GenMutationModelError> {
    // This will generate a MutationModel and write it to the output file 
    // provided.
    // params:
    //   fasta_map: The database for the fasta file, with links back to the 
    //      original sequences
    //   filtered_mutations: The mutations from the original VCF, after any 
    //      appropriate bed filtering has been applied
    //   bed_table: The hash map representing the input bed
    //   output_file: Where to write the output model
    //
    // This will store the count of trinucleotide combo A transitioning to combo B
    // This is the raw data for the trinuc model.
    let mut trinuc_transition_count: HashMap<(TrinucFrame, TrinucFrame), usize> = HashMap::new();
    let mut snp_count = 0;
    let mut snp_transition_count: HashMap<(Nucleotide, Nucleotide), usize> = HashMap::new();
    let mut insertion_count: HashMap<usize, usize> = HashMap::new();
    let mut deletion_count: HashMap<usize, usize> = HashMap::new();
    let mut homozygous_count = 0;
    let (trinuc_count, bed_track_len) = count_trinculeotides (
        &fasta_map,
        &bed_table,
    )?;
    if trinuc_count.is_empty() {
        error!("Trinucleotide counts were empty");
        return Err(
            GenMutationModelError::TrinucCountError(
                "Trinuc counts are empty. Unkown error".to_string()
            )
        )
    }
    let mut total_reflen = 0;
    for contig in &fasta_map.contigs{
        let countable_regions = contig.get_non_n_regions()?;
        for region in countable_regions {
            // I think the coordinates are forward on all regions
            total_reflen += region.end - region.start;
        }
        let matching_variants = match filtered_mutations.get(&contig.name) {
            Some(v) => v,
            None => {
                debug!("No variants found for {}", contig.name);
                continue
            }
        };
        if matching_variants.is_empty() {
            debug!("No variants found for {}", contig.name);
            continue
        }
        for variant in matching_variants {
            match variant.variant_type {
                VariantType::SNP => {
                    snp_count += 1;
                    // VCF POS is 1-based; skip variants too close to contig edges
                    if variant.location < 2 {
                        debug!("Skipping edge variant at position {}", variant.location);
                        continue;
                    }
                    // Trinucleotide centered on the variant: [pos-1, pos, pos+1] (0-based)
                    // => get_block_subseq(pos-2, pos+1) with 1-based VCF pos
                    let trinuc_result = contig.get_block_subseq(variant.location-2, variant.location+1);
                    let trinuc_to_analyze = match trinuc_result {
                        Ok(t) => t,
                        Err(_) => {
                            debug!("Skipping edge variant at position {} (out of contig bounds)", variant.location);
                            continue;
                        }
                    };
                    let ref_array: [&Nucleotide; 3] = trinuc_to_analyze.iter().collect_array().unwrap();
                    let ref_frame = TrinucFrame::from(ref_array);
                    let alt_array = [ref_array[0], &variant.alternate[0], ref_array[2]];
                    let alt_frame = TrinucFrame::from(alt_array);
                    if *ref_array[1] != variant.reference[0] {
                        return Err(GenMutationModelError::BaseMismatch(format!("{:?}", variant)))
                    }
                    *trinuc_transition_count.entry((ref_frame, alt_frame)).or_default() += 1;
                    *snp_transition_count.entry((variant.reference[0], variant.alternate[0])).or_default() += 1;
                },
                VariantType::Insertion => {
                    let variant_len = variant.alternate.len() - variant.reference.len();
                    *insertion_count.entry(variant_len).or_default() += 1;
                },
                VariantType::Deletion => {
                    let variant_len = variant.reference.len() - variant.alternate.len();
                    *deletion_count.entry(variant_len).or_default() += 1;
                },
                _ => debug!("Unknown variant type, skipping for this analyis."),
            }
            // here we could get the variant population frequency from info extracted from the VCF
            // let var_pop_freq = find_caf(&variant.alternate);
            match variant.genotype {
                Genotype::Homozygous => {
                    homozygous_count += 1
                }
                Genotype::Heterozygous => {
                    // Do nothing
                }
            }
        }
    }
    // Compute probabilities
    let mut trinuc_mut_prob: HashMap<TrinucFrame, f64> = HashMap::new();
    let mut trinuc_trans_prob: HashMap<(TrinucFrame, TrinucFrame), f64> = HashMap::new();
    let mut snp_trans_frequency: HashMap<(Nucleotide, Nucleotide), f64> = HashMap::new();
    for (frame, count) in trinuc_count {
        // No need to run this loop if we never saw this frame in the first place
        if count == 0 {
            trinuc_mut_prob.insert(frame, 0.0);
            // nothing to look for
            continue;
        }
        // Thi number of times this frame appeared in any transition in the dataset.
        let mut frame_count = 0;
        for (key, value) in &trinuc_transition_count {
            if key.0 == frame {
                frame_count += *value;
            }
        }
        // The probability that this particular frame mutates.
        if count > 0 {
            trinuc_mut_prob.insert(frame, (frame_count as f64) / (count as f64));
        } else {
            // No dividsies by zeroesies. A count of zero means this particular trinucleotide did not appear in the dataset.
            trinuc_mut_prob.insert(frame, 0.0);
        }
        // So now we caluclate the transition, based on the probability that this trinuc
        // mutated, from this frame to each particular frame that it mutated to.
        if frame_count > 0 {
            // No need to run this loop if we never saw the trinuc.
            for (frame_tuple, count,) in &trinuc_transition_count{
                if frame_tuple.0 == frame {
                    // We already checked that frame_count > 0, so we're good to divide here.
                    trinuc_trans_prob.insert(
                        frame_tuple.clone(), 
                        (*count as f64)/(frame_count as f64),
                    );
                }
            }
        }
    }
    // Now we check just the middle of the frame, to see the base chance of that base transitioning
    for nuc1 in ALLOWED_NUCS.clone() {
        // this could probably be condensed into a mapping plus .iter().sum()
        let rolling_total = {
            let mut sub_tot = 0;
            for nuc2 in ALLOWED_NUCS.clone() {
                if snp_transition_count.contains_key(&(nuc1, nuc2)) {
                    sub_tot += snp_transition_count[&(nuc1, nuc2)]
                }
            };
            sub_tot
        };
        // Sucks we have to loop this twice. Maybe if we structure the map better?
        // but it's only O(4) so it's not bad.
        if rolling_total > 0 {
            // No divide by zero
            for nuc2 in ALLOWED_NUCS.clone() {
                if snp_transition_count.contains_key(&(nuc1, nuc2)) {
                    // probability that if nuc1 mutates into anything, it mutates into nuc2.
                    snp_trans_frequency.insert(
                        (nuc1, nuc2).clone(),
                        (snp_transition_count[&(nuc1, nuc2)] as f64) / (rolling_total as f64)
                    );
                }
            }
        }
    }
    // compute average snp and indel frequencies
    let total_insertions = {
        let mut subtotal = 0;
        for (_, v) in &insertion_count {
            subtotal += *v;
        }
        subtotal
    };
    let total_deletions = {
        let mut subtotal = 0;
        for (_, v) in &deletion_count {
            subtotal += *v;
        }
        subtotal
    };
    // There may have been unknown variants in our dataset
    let allowed_variant_count: f64 = (snp_count + total_insertions + total_deletions) as f64;
    let average_snp_frequency = (snp_count as f64) / allowed_variant_count;
    let average_deletion_frequency = (total_deletions as f64) / allowed_variant_count;
    let average_insertion_frequency = (total_insertions as f64) / allowed_variant_count;
    let variant_probs = vec![
        average_snp_frequency, 
        average_insertion_frequency, 
        average_deletion_frequency
    ];
    let homozygous_frequency = {
        if homozygous_count > 0 {
            (homozygous_count as f64) / allowed_variant_count
        } else {
            // A very small percentage allowed. Maybe this could be a parameter.
            0.001 / allowed_variant_count
            // This number comes from the original NEAT
        }
    };
    let average_mutation_rate = {
        if !bed_table.is_empty() {
            allowed_variant_count / (bed_track_len as f64)
        } else {
            allowed_variant_count / (total_reflen as f64)
        }
    };

    let ins_lengths: Vec<usize> = insertion_count.keys().cloned().collect();
    let ins_weights: Vec<f64> = insertion_count
        .values()
        .cloned()
        .map(|x| x as f64)
        .collect();
    let del_lengths: Vec<usize> = deletion_count.keys().cloned().collect();
    let del_weights: Vec<f64> = deletion_count
        .values()
        .cloned()
        .map(|x| x as f64)
        .collect();

    let result = MutationModel::from_raw_data(
        average_mutation_rate,
        homozygous_frequency,
        variant_probs,
        snp_trans_frequency,
        trinuc_mut_prob,
        trinuc_trans_prob,
        ins_lengths,
        ins_weights,
        del_lengths,
        del_weights,
    );

    match result {
        Ok(mut_model) => {
            mut_model.write_to_file(output_file)?;
            info!("Mutation model success! Wrote model to {:?}", output_file);
            Ok(())
        },
        Err(error) => {
            return Err(GenMutationModelError::MutModelError(error))
        }
    }
}

fn count_trinculeotides(
    fasta_map: &FastaMap,
    bed_table: &HashMap<String, Vec<BedRecord>>,
) -> Result<(HashMap<TrinucFrame, usize>, usize), GenMutationModelError> {
    // Counts the frequency of the various trinucleotide combinations in the dataset
    // fasta_map - Database for reading back the sequenced data.
    // bed_table - A table for the bed records

    // HashMap of trinuce to count, e.g., "AAA": 12
    let mut trinuc_count: HashMap<TrinucFrame, usize> = HashMap::new();
    let mut bed_track_len: usize = 0;
    if !bed_table.is_empty() {
        // There is actual bed input
        info!("Counting trinucleotide combinations in bed regions");
        for (chrom, bed_vec) in bed_table {
            let contig_block = fasta_map.retrieve_contig(chrom.to_string())
                .expect(&format!("Error retrieving contig {}", chrom));
            let regions_of_interest = bed_vec;
            for region in regions_of_interest {
                // Make sure our region is long enough. We need at least 3 bases to get any trinucleotides.
                if (region.end - region.start) >= 3 {
                    // Add to track length
                    bed_track_len += region.end - region.start;
                    // iterate over the regions, grab trinucleotides, profit
                    // The first trincu is one base in. The last is one base from the end
                    for i in (region.start+1)..(region.end-1) {
                        let trinuc = contig_block.get_block_subseq(i-1, i+2)?;
                        let frame = TrinucFrame::from((trinuc[0], trinuc[1], trinuc[2]));
                        *trinuc_count.entry(frame.to_owned()).or_default() += 1;
                    }
                } 
            }
        }
    } else {
        // No bed input
        info!("Counting trinucleotide combinations in bed regions");
        for contig in &fasta_map.contigs {
            let contig_block = fasta_map.retrieve_contig(contig.name.clone())
                .expect(&format!("Error retrieving contig {}", &contig.name));
            let regions_of_interest = contig_block.get_non_n_regions()
                .expect("Error retrieving non-n regions!");
            for region in regions_of_interest {
                // iterate over the regions, grab trinucleotides, profit
                // The first trincu is one base in. The last is one base from the end
                for i in (region.start+1)..(region.end-1) {
                    let trinuc = contig_block.get_block_subseq(i-1, i+2)?;
                    let frame = TrinucFrame::from((trinuc[0], trinuc[1], trinuc[2]));
                    *trinuc_count.entry(frame.to_owned()).or_default() += 1;
                }
            }
        }
    }
    Ok((trinuc_count, bed_track_len))
}

#[cfg(test)]
mod tests {
    use super::*;
    use common::{
        file_tools::{fasta_reader::read_fasta, vcf_tools::read_vcf},
        models::mutation_model::MutationModel,
    };
    use tempfile::tempdir;

    #[test]
    fn test_runner_with_snps() {
        // Exercises the full runner pipeline with SNP-only input, which also
        // exercises the IndelModel::default() fallback path in MutationModel::from_raw_data.
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let reference = PathBuf::from(format!("{}/test_data/references/H1N1.fa", manifest_dir));
        let vcf_path = PathBuf::from(format!("{}/test_data/vcfs/small_snps.vcf", manifest_dir));
        let temp = tempdir().unwrap();
        let fasta_map = read_fasta(&reference, None, 0, &temp, None).unwrap();
        let mutations = read_vcf(vcf_path).unwrap();
        let out_dir = tempdir().unwrap();
        let output_file = out_dir.path().join("test_model.json.gz");
        runner(fasta_map, mutations, HashMap::new(), &output_file).unwrap();
        assert!(output_file.exists());
        let model = MutationModel::from_file(&output_file).unwrap();
        assert!(model.mutation_rate > 0.0, "mutation_rate should be positive");
        // 1 of 3 SNPs has GT=1/1 (homozygous), so homozygous_frequency should be > 0
        assert!(model.homozygous_frequency > 0.0, "homozygous_frequency should be positive");
    }
}