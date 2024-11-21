extern crate clap;
extern crate itertools;
extern crate log;
extern crate serde;
extern crate serde_json;
extern crate serde_yaml;
extern crate simplelog;
extern crate fasta_reader;
extern crate simple_rng;

use common;
pub mod utils;
mod data;

use chrono::Utc;
use clap::Parser;
use log::*;
use simplelog::*;
use std::fs::File;
use utils::cli;
use utils::config::{build_config_from_args, read_config_yaml};
use common::file_tools::check_parent;
use utils::runner::run_neat;
use simple_rng::NeatRng;

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
    let log_destination = check_parent(&args.log_dest, true).unwrap();
    // Set up the logger for the run
    CombinedLogger::init(vec![
        #[cfg(feature = "termcolor")]
        TermLogger::new(
            level_filter,
            Config::default(),
            TerminalMode::Mixed,
            ColorChoice::Auto,
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
    let mut seed_vec: Vec<String> = Vec::new();
    let user_input = (&config.rng_seed).clone();
    let seed_vec = match user_input {
        Some(u) => {
            // For read it as a string, hopefully
            let raw_seed = u.to_owned();
            for seed_term in raw_seed.split_whitespace() {
                seed_vec.push(seed_term.to_string());
            }
            info!("Seed string to regenerate these exact results: {}", &raw_seed);
            seed_vec
        },
        _ => {
            // Since no seed was provided, we'll use a datetime stamp with nanoseconds
            // The seed can be any space separated or tab separated series of strings
            // e.g., "Every good boy does Fine"
            // seeds are case-sensitive
            let timestamp = Utc::now().format("%Y %m %d %H %M %S %f").to_string();
            for item in timestamp.split_whitespace() {
                seed_vec.push(item.to_string());
            }
            info!("Seed string to regenerate these exact results: {}", &timestamp);
            seed_vec
        },
    };
    let rng: NeatRng = NeatRng::new_from_seed(seed_vec);
    // run the generate reads main script
    run_neat(&config, rng)
        .unwrap_or_else(|error| panic!("Neat encountered a problem: {:?}", error))
}
