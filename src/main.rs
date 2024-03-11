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

fn main() {

    TermLogger::init(
        LevelFilter::Trace,
        Config::default(),
        TerminalMode::Stdout,
        ColorChoice::Auto,
    ).unwrap();

    let mut rng = thread_rng();

    info!("Begin processing");

    let args = cli::Cli::parse();

    let config = if args.config != "" {
        info!("Using Configuration file input: {}", &args.config);
        read_config_yaml(args.config)
    } else {
        info!("Using command line arguments.");
        debug!("Command line args: {:?}", &args);
        Ok(build_config_from_args(args).expect("Problem reading configuration yaml file"))
    }.unwrap();

    let output_file = format!("{}/{}", config.output_dir, config.output_prefix);

    info!("Mapping fasta file: {}", &config.reference);
    let (fasta_map, fasta_order) = read_fasta(&config.reference);
    // todo:
    // need to add this twice, produce two mutated fastas, or at least 2 separate mutation
    // datasets, each with half the mutation rate. Going to mean twice as much memory needed for
    // fasta creation, which isn't ideal
    info!("Mutating fasta");
    let mutated_map: Box<HashMap<String, Vec<Vec<u8>>>> = mutate_fasta(
        &fasta_map,
        config.ploidy,
        &mut rng
    );

    info!("Outputting fasta files");
    if config.produce_fasta == true {
        write_fasta(
            &mutated_map,
            &fasta_order,
            &output_file,
            config.ploidy
        ).expect("Problem writing fasta file");
    }

    let mut read_sets: HashSet<Vec<u8>> = HashSet::new();
    for (_name, sequences) in mutated_map.iter() {
        // defined as a set of read sequences that should cover the mutated sequence `coverage` number of times
        let data_set = generate_reads(
            &sequences,
            &config.read_len,
            &config.coverage,
            &mut rng
        );

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

