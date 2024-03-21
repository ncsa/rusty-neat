use std::collections::HashSet;
use std::fmt::{Display, Formatter};
use log::info;
use rand::prelude::SliceRandom;
use rand::rngs::ThreadRng;
use utils::config::RunConfiguration;
use utils::fasta_tools::{read_fasta, write_fasta};
use utils::fastq_tools::write_fastq;
use utils::make_reads::generate_reads;
use utils::mutate::mutate_fasta;
use utils::vcf_tools::write_vcf;

#[derive(Debug)]
pub enum NeatErrors {
    FastaError,
    FastqError,
    VcfError,
}

impl Display for NeatErrors {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        todo!()
    }
}

pub fn run_neat(config: Box<RunConfiguration>, mut rng: &mut ThreadRng) -> Result<(), NeatErrors>{
    // Create the prefix of the files to write
    let output_file = format!("{}/{}", config.output_dir, config.output_prefix);

    // Reading the reference file into memory
    info!("Mapping reference fasta file: {}", &config.reference);
    let (fasta_map, fasta_order) = read_fasta(&config.reference);

    // Mutating the reference and recording the variant locations.
    info!("Mutating reference.");
    let (mutated_map, variant_locations) = mutate_fasta(
        &fasta_map,
        &mut rng
    );

    if config.produce_fasta {
        info!("Outputting fasta file");
        write_fasta(
            &mutated_map,
            &fasta_order,
            config.overwrite_output,
            &output_file,
        )
            .or(Err(NeatErrors::FastaError))
            .expect("Error writing fasta file");
    }

    if config.produce_vcf {
        info!("Writing vcf file");
        write_vcf(
            &variant_locations,
            &fasta_order,
            config.ploidy,
            &config.reference,
            config.overwrite_output,
            &output_file,
            &mut rng
        )
            .or(Err(NeatErrors::VcfError))
            .expect("Error writing VCF file");
    }

    let mut read_sets: HashSet<Vec<u8>> = HashSet::new();
    for (_name, sequence) in mutated_map.iter() {
        // defined as a set of read sequences that should cover
        // the mutated sequence `coverage` number of times
        let data_set = generate_reads(
            sequence,
            &config.read_len,
            &config.coverage,
            config.paired_ended,
            config.fragment_mean,
            config.fragment_st_dev,
            &mut rng
        );

        read_sets.extend(*data_set);
    }

    if config.produce_fastq {
        info!("Shuffling output fastq data");
        let mut outsets: Box<Vec<&Vec<u8>>> = Box::new(read_sets.iter().collect());
        outsets.shuffle(&mut rng);

        info!("Writing fastq");
        write_fastq(
            &output_file,
            config.overwrite_output,
            config.paired_ended,
            *outsets,
        )
            .or(Err(NeatErrors::FastqError))
            .expect("Problem writing fastq file");
        info!("Processing complete")
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs;
    use rand::thread_rng;
    use std::path::Path;
    use super::*;

    #[test]
    fn test_runner() -> Result<(), NeatErrors>{
        let config = RunConfiguration::build().build();
        run_neat(
            Box::new(config),
            &mut thread_rng(),
        ).unwrap();
        let fastq_file = Path::new("neat_out_r1.fastq").canonicalize().unwrap();
        fs::remove_file(fastq_file).unwrap();
        Ok(())
    }
}