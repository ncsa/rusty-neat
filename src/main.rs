mod utils;
use utils::cli;
use std::{env, process};
use clap::Parser;
use simplelog::*;
use rand::{thread_rng};
use rand::prelude::*;
use crate::utils::config::{read_config_yaml, build_config_from_args, Config};

fn main() {
    info!("Begin processing");

    let args = cli::Cli::parse();

    let config = if args.config != "" {
        info!("Using Configuration file input.");
        read_config_yaml(args.config)
    } else {
        info!("Using command line arguments.");
        build_config_from_args(args)
    }.unwrap();

    println!(
        "Running neat on {} with {} bp read length and a coverage of {}.\n> Will output file: {}",
        config.reference, config.read_len, config.coverage, format!("{}/{}", config.output_dir, config.output_prefix)
    )



    // let config = Config::build(&args).unwrap_or_else(|err| {
    //     println!("Problem parsing arguments: {err}");
    //     process::exit(1)
    // });
    //
    // println!("In file {}", config.fasta_path);
    //
    // if let Err(e) = rusty_neat::run(config) {
    //     println!("Application error: {e}");
    //     process::exit(1);
    // }
    //
    // let sequence = String::from("ACATCTACACGGAGCACGACGCACCNNNNNNccccNNNNCACTGCAGTGCATGACG");
    // logger::info!("Input sequence: {}", sequence);
    // let non_n_map = map_non_n_regions(&sequence);
    // if !non_n_map.contains(&true) {
    //     logger::info!("Input sequence contained no non-n values. Quitting.");
    //     std::process::exit(0)
    // }
    // logger::debug!("Non-n-map = {:?}", non_n_map);
    // let position = find_random_non_n(&non_n_map).unwrap();
    // logger::info!("Found a position! It is {}",position);
    // let code = &sequence.as_bytes()[position];
    // let nucleotide: char = *code as char;
    // logger::info!("Nucleotide to mutate: {:?}", nucleotide);
    // let new_nucleotide = generate_new_nucleotide(&nucleotide);
    // let result: String = sequence.chars().enumerate().map(|(i, c)| if i == position { new_nucleotide } else { c }).collect();
    // logger::info!("New sequence = {}", result);
    // logger::info!("End program");
}

