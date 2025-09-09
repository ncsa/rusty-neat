extern crate clap;
extern crate itertools;
extern crate log;
extern crate serde;
extern crate serde_json;
extern crate simplelog;
extern crate simple_rng;

use common;
pub mod utils;
pub mod errors;

use clap::Parser;
use log::*;
use simplelog::*;
use std::fs::File;
use utils::cli;
use crate::{errors::GenerateReadsErrors, utils::config::RunConfiguration};
use common::file_tools::folder_tools::check_parent;
use utils::runner::run_neat;
use simple_rng::NeatRng;

fn main() -> Result<(), GenerateReadsErrors> {
    info!("Begin processing");
    // parse the arguments from the command line
    let args = cli::Cli::parse();
    if args.config.is_empty() && args.reference.is_empty() {
        panic!("Must supply either a configuration file or a reference file to run NEAT!")
    }
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
        TermLogger::new(
            LevelFilter::Info,
            Config::default(),
            TerminalMode::Mixed,
            ColorChoice::Auto,
        ),
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
        RunConfiguration::from_yaml_file(args.config)?
    } else {
        info!("Using command line arguments.");
        debug!("Command line args: {:?}", &args);
        RunConfiguration::from_args(args)?
    };
    // Generate the RNG used for this run. If one was given in the config file, use that, or else
    // use thread_rng to generate a random seed, then seed using a SeedableRng based on StdRng
    let mut rng: NeatRng = NeatRng::new_from_seed(&config.seed_vec)
        .expect("Neat failed during rng creation!");
    // run the generate reads main script
    let result = run_neat(&Box::new(config), &mut rng);
    match result {
        Ok(()) => Ok(()),
        Err(error) => return Err(GenerateReadsErrors::MainError(Box::new(error))),
    }
}