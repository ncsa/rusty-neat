extern crate rand;
extern crate clap;
extern crate log;
extern crate simplelog;
extern crate serde_yaml;
extern crate rand_distr;
extern crate itertools;
extern crate rand_core;
extern crate rand_chacha;

mod utils;

use std::collections::HashMap;
use std::fs::File;
use clap::{Parser};
use log::*;
use simplelog::*;
use rand::SeedableRng;
use rand::thread_rng;
use rand_core::RngCore;
use utils::cli;
use utils::config::{read_config_yaml, build_config_from_args};
use utils::file_tools::check_parent;
use utils::runner::run_neat;
use utils::neat_rng::NeatRng;

fn main() {

    info!("Begin processing");
    // parse the arguments from the command line
    let args = cli::Cli::parse();

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
        )
    };

    // Check that the parent dir exists
    let log_destination = check_parent(&args.log_dest).unwrap();

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
        )
    ]).unwrap();

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
    let mut rng: NeatRng = SeedableRng::seed_from_u64(seed);

    run_neat(config, &mut rng).unwrap_or_else(|error| {
        panic!("Neat encountered a problem: {:?}", error)
    })
}

