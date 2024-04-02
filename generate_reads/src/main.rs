extern crate clap;
extern crate itertools;
extern crate log;
extern crate rand;
extern crate rand_chacha;
extern crate rand_core;
extern crate rand_distr;
extern crate serde;
extern crate serde_json;
extern crate serde_yaml;
extern crate simplelog;

use common;
pub mod utils;

use clap::Parser;
use log::*;
use rand::thread_rng;
use rand::SeedableRng;
use rand_core::RngCore;
use simplelog::*;
use std::collections::HashMap;
use std::fs::File;
use rand_chacha::ChaCha20Rng;
use utils::cli;
use utils::config::{build_config_from_args, read_config_yaml};
use common::file_tools::check_parent;
use utils::runner::run_neat;

fn main() {
    info!("Begin processing");
    // parse the arguments from the command line
    let args = cli::Cli::parse();
    // log filter
    let level_filter = match args.log_level.to_lowercase().as_str() {
        "trace" => LevelFilter::Trace,
        "debug" => LevelFilter::Debug,
        "info" => LevelFilter::Info,
        "warn" => LevelFilter::Warn,
        "error" => LevelFilter::Error,
        "off" => LevelFilter::Off,
        _ => panic!(
            "Unknown log level, please set to one of \
            Trace, Debug, Info, Warn, Error, or Off (case insensitive)."
        ),
    };
    // Check that the parent dir exists
    let log_destination = check_parent(&args.log_dest).unwrap();
    // Set up the logger for the run
    CombinedLogger::init(vec![
        #[cfg(feature = "termcolor")]
        TermLogger::new(
            level_filter,
            Config::default(),
            TerminalMode::Stdout,
            ColorChoice::Always,
        ),
        #[cfg(not(feature = "termcolor"))]
        SimpleLogger::new(LevelFilter::Trace, Config::default()),
        WriteLogger::new(
            level_filter,
            Config::default(),
            File::create(log_destination).unwrap(),
        ),
    ])
    .unwrap();
    // set up the config struct based on whether there was an input config. Input config
    // overrides any other inputs.
    let config = if args.config != "" {
        info!("Using Configuration file input: {}", &args.config);
        read_config_yaml(args.config)
    } else {
        info!("Using command line arguments.");
        debug!("Command line args: {:?}", &args);
        build_config_from_args(args)
    };
    // Generate the RNG used for this run. If one was given in the config file, use that, or else
    // use thread_rng to generate a random seed, then seed using a SeedableRng based on StdRng
    let seed: u64;
    if !config.rng_seed.is_none() {
        seed = config.rng_seed.unwrap();
    } else {
        // We pick a random u64 to use as a seed
        let mut seed_rng = thread_rng();
        seed = seed_rng.next_u64();
    }
    info!("Generating random numbers using the seed: {}", seed);
    let rng: ChaCha20Rng = SeedableRng::seed_from_u64(seed);
    // run the generate reads main script
    run_neat(config, rng)
        .unwrap_or_else(|error| panic!("Neat encountered a problem: {:?}", error))
}
