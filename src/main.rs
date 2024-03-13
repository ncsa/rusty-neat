extern crate rand;
extern crate clap;
extern crate log;
extern crate simplelog;
extern crate serde_yaml;

mod utils;

use std::collections::{HashMap, HashSet};
use clap::{Parser};
use log::*;
use simplelog::*;
use rand::thread_rng;
use rand::prelude::*;

use utils::cli;
use utils::fasta_tools::{read_fasta, write_fasta};
use utils::config::{read_config_yaml, build_config_from_args};
use utils::mutate::mutate_fasta;
use utils::make_reads::generate_reads;
use utils::fastq_tools::write_fastq;
use utils::vcf_tools::write_vcf;

fn main() {

    TermLogger::init(
        LevelFilter::Trace,
        Config::default(),
        TerminalMode::Stdout,
        ColorChoice::Auto,
    ).unwrap();

    let mut rng = thread_rng();

    info!("Begin processing");
    // parse the arguments from the command line
    let args = cli::Cli::parse();

    // set up the config struct based on whether there was an input config. Input config
    // overrides any other inputs.
    let config = if args.config != "" {
        info!("Using Configuration file input: {}", &args.config);
        read_config_yaml(args.config)
    } else {
        info!("Using command line arguments.");
        debug!("Command line args: {:?}", &args);
        Ok(build_config_from_args(args).expect("Problem reading configuration yaml file"))
    }.unwrap();

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

    if config.produce_vcf {
        info!("Writing vcf file");
        write_vcf(
            &fasta_map,
            &variant_locations,
            &fasta_order,
            config.ploidy,
            &config.reference,
            &output_file,
            &mut rng).expect("Error writing vcf file")
    }

    if config.produce_fasta {
        info!("Outputting fasta file");
        write_fasta(
            &mutated_map,
            &fasta_order,
            &output_file,
        ).expect("Problem writing fasta file");
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
            &mut rng
        );

        // todo: we need to keep track of these reads better for paired endedness
        read_sets.extend(*data_set);
    }

    info!("Shuffling output fastq data");
    let mut outsets: Box<Vec<&Vec<u8>>> = Box::new(read_sets.iter().collect());
    outsets.shuffle(&mut rng);

    info!("Writing fastq");
    write_fastq(
        &output_file,
        *outsets,
    ).expect("Problem writing fastq file");
    info!("Processing complete")
}

