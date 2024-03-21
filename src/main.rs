extern crate rand;
extern crate clap;
extern crate log;
extern crate simplelog;
extern crate serde_yaml;
extern crate rand_distr;
extern crate itertools;

mod utils;

use std::collections::HashMap;
use std::fs::File;
use clap::{Parser};
use log::*;
use simplelog::*;
use rand::thread_rng;

use utils::cli;
use utils::config::{read_config_yaml, build_config_from_args};
use utils::file_tools::check_parent;
use utils::runner::run_neat;

fn main() -> Result<(), std::fmt::Error> {

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

    let mut rng = thread_rng();

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

    run_neat(config, &mut rng).unwrap();
    Ok(())
}

