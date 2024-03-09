extern crate rand;
extern crate clap;
extern crate log;
extern crate simplelog;

mod utils;

use std::{env, process};
use std::collections::{HashMap, HashSet};
use clap::{Parser};
use log::*;
use simplelog::*;
use rand::thread_rng;
use rand::prelude::*;

use utils::cli;
use utils::fasta_tools::read_fasta;
use utils::config::{read_config_yaml, build_config_from_args, RunConfiguration};
use utils::mutate::mutate_fasta;
use utils::make_reads::generate_reads;
use utils::fastq_writer::write_fastq;

fn main() {

    TermLogger::init(
        LevelFilter::Trace,
        Config::default(),
        TerminalMode::Stdout,
        ColorChoice::Auto,
    )
        .unwrap();

    let mut rng = thread_rng();

    info!("Begin processing");

    let args = cli::Cli::parse();

    let config = if args.config != "" {
        info!("Using Configuration file input: {}", &args.config);
        read_config_yaml(args.config)
    } else {
        info!("Using command line arguments.");
        debug!("Command line args: {:?}", &args);
        build_config_from_args(args)
    }.unwrap();

    let output_file = format!("{}/{}.fastq", config.output_dir, config.output_prefix);
    info!(
        "Running neat on {} with {} bp read length and a coverage of {}.\n> Will output file: {}",
        config.reference, config.read_len, config.coverage, output_file
    );

    info!("Mapping fasta file: {}", &config.reference);
    let fasta_map = read_fasta(&config.reference);
    // todo:
    // need to add this twice, produce two mutated fastas, or at least 2 separate mutation
    // datasets, each with half the mutation rate. Going to mean twice as much memory needed for
    // fasta creation, which isn't ideal
    let mutated_map = mutate_fasta(&fasta_map, config.ploidy, rng);
    let mut rng = thread_rng();
    let strict_read_length: Option<bool> = Option::from(true);

    let mut read_sets: HashSet<Vec<u8>> = HashSet::new();
    for (name, sequences) in mutated_map.iter() {
        // defined as a set of read sequences that should cover the mutated sequence `coverage` number of times
        let data_set = generate_reads(
            &sequences,
            &config.read_len,
            &config.coverage,
            &mut rng,
            strict_read_length,
        );

        read_sets.extend(*data_set);
    }

    info!("Shuffling output fastq data");
    let mut outsets: Box<Vec<&Vec<u8>>> = Box::new(read_sets.iter().collect());
    outsets.shuffle(&mut rng);

    info!("Writing fastq: {}", output_file);
    write_fastq(
        &output_file,
        *outsets,
    ).unwrap();
    info!("Processing complete")
}

