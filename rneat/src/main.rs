extern crate clap;
extern crate itertools;
extern crate log;
extern crate serde;
extern crate serde_json;
extern crate simplelog;
extern crate simple_rng;

pub mod filter_reads;
pub mod gen_reads;
pub mod gen_mut_model;

use std::env;
use std::time;
use common::{self, file_tools::file_io::create_output_file};

use common::file_tools::folder_tools::check_parent;
use simplelog::{
    ColorChoice, 
    CombinedLogger, 
    Config, 
    TermLogger, 
    TerminalMode, 
    WriteLogger
};
use thiserror::Error;
use log::*;
use clap::{value_parser, Arg, ArgAction, Command};
use std::path::PathBuf;

use crate::gen_mut_model::errors::GenMutationModelError;
use crate::{
    filter_reads::errors::FilterReadsError, 
    gen_reads::errors::GenerateReadsErrors,
};
/// This script parses arguments and checks them before submitting to the submodules, which currently 
/// include `gen-reads` and `filter-files`. As more are added, this can be expanded or refactored
/// argument parsing, timing, and other concerns will occur here.
#[derive(Error, Debug)]
pub enum NeatErrors {
	// Errors specific to each submodule
    #[error("Error while generating read dataset {0}")]
    GenerateReadsError(#[from] GenerateReadsErrors),
    #[error("Error while filtering reads with bed file {0}")]
    FilterReadsError(#[from] FilterReadsError),
    #[error("Error while generating mutation model {0}")]
    GenMutModelError(#[from] GenMutationModelError),
}

fn neat_commands() -> [Command; 3] {
	// These are the submodule commands. Any new commands added should go here.
    let configuration_arg = Arg::new("configuration_yaml")
        .long("configuration-yaml")
        .short('c')
        .help("Path to configuration file.")
        .action(ArgAction::Set)
        .exclusive(true)
        .default_missing_value("")
        .value_parser(value_parser!(PathBuf))
        .required(true);
    [
        Command::new("gen-reads")
            .about("Generates reads for an input dataset")
            .arg_required_else_help(true)
            .arg(
                &configuration_arg
            ),
        Command::new("filter-reads")
            .about("Filters the output of gen-reads")
            .arg_required_else_help(true)
            .arg(
                &configuration_arg
            ),
        Command::new("gen-mut-model")
            .about("Generate a mutation model based on an input mutations file")
            .arg_required_else_help(true)
            .arg(
                &configuration_arg
            )
    ]
}

fn main() -> Result<(), NeatErrors> {
    // parse the arguments from the command line
    let cmd = Command::new(env!("CARGO_CRATE_NAME"))
        .multicall(true)
        .subcommand(
            Command::new("rneat")
                .arg_required_else_help(true)
                .subcommand_value_name("SUB-COMMAND")
                .subcommand_help_heading("SUB-COMMANDS")
                .arg(
					// This arg controls the amount of stuff shown in the display log.
                    Arg::new("log_level")
                        .long("log-level")
                        .help("Sets the log level for the display log.")
                        .action(ArgAction::Set)
                        .num_args(0..=1)
                        .value_parser(["trace", "debug", "info", "warn", "error", "off"])
                        .default_value("trace")
                        .default_missing_value("trace")
                )
                .arg(
					// This arg controls where the log is written. The default is <current_working_dir>/.neat.log
                    Arg::new("log_dest")
                        .long("log-dest")
                        .help("Sets the log destination (full path with full filename) for the written log")
                        .action(ArgAction::Set)
                        .value_parser(value_parser!(PathBuf))
                        .default_missing_value("")
                )
                .subcommands(neat_commands())
        )
        .subcommands(neat_commands());
        
	// This parses the command arguments
    let matches = cmd.get_matches();
    let mut subcommand = matches.subcommand();
    let mut level_filter = "trace".to_string();
    let mut log_dest = env::current_dir().unwrap();
    log_dest.push(".neat.log");
    println!("log dest = {:?}", log_dest.display());
    if let Some(("rneat", cmd)) = subcommand {
        if cmd.contains_id("log_level") {
            level_filter = cmd
                .get_one::<String>("log_level")
                .expect("Must enter value for log-level")
                .to_string();
        } 
        if cmd.contains_id("log_dest") {
            log_dest = cmd
                .get_one::<PathBuf>("log_dest")
                .expect("Must provide a path with log-dest")
                .to_path_buf();
        }
        subcommand = cmd.subcommand();
    }

    // log filter
    let lev_filter_choice = {
        match level_filter.to_lowercase().as_str() {
            "trace" => LevelFilter::Trace,
            "debug" => LevelFilter::Debug,
            "info" => LevelFilter::Info,
            "warn" => LevelFilter::Warn,
            "error" => LevelFilter::Error,
            "off" => LevelFilter::Off,
            _ => unreachable!(
                "Parser should only allow one of the above values."
            ), 
        }
    };
    
    // Check that the parent dir exists
    let log_destination = check_parent(&log_dest, true).expect("Log destination parent path doesn't exist");
    let log_file = create_output_file(log_destination, true)
        .expect(&format!("Error creating log file: {:?}", log_destination));
    // Set up the logger for the run
    CombinedLogger::init(vec![
        TermLogger::new(
            LevelFilter::Info,
            Config::default(),
            TerminalMode::Mixed,
            ColorChoice::Auto,
        ),
        WriteLogger::new(
            lev_filter_choice,
            Config::default(),
            log_file,
        ),
    ]).expect("Error creating log file!");
    info!("////////////// Welcome to rneat! \\\\\\\\\\\\\\\\\\\\\\\\\\\\");
    let start = time::Instant::now();
    match subcommand {
        Some(("gen-reads", _)) => {
            if let Some(("gen-reads", cmd)) = subcommand {
                if cmd.contains_id("configuration_yaml") {
                    let file = cmd.get_one::<PathBuf>("configuration_yaml")
                            .expect("Must provide a path with configuration-yaml")
                            .to_path_buf();

                    if !file.is_file() {
                        return Err(
                            NeatErrors::FilterReadsError(
                                FilterReadsError::CliError(
                                    "Must supply either a configuration file or a reference file to run NEAT!".to_string()
                                )
                            )
                        )
                    }
                    info!("Running gen-reads to generate read data.");
                    let result = gen_reads::main(&file);
                    match result {
                        Err(error) => return Err(NeatErrors::GenerateReadsError(error)),
                        Ok(()) => info!("rneat gen-reads completed successfully"),
                    }
                } 
            }
        },
        Some(("filter-reads", _)) => {
            if let Some(("filter-reads", cmd)) = subcommand {
				info!("Running rneat filter-reads");
                // extract the flags and values
                if cmd.contains_id("configuration_yaml") {
                    let file = cmd.get_one::<PathBuf>("configuration_yaml")
                            .expect("Must provide a path with configuration-yaml")
                            .to_path_buf();

                    if !file.is_file() {
                        return Err(
                            NeatErrors::FilterReadsError(
                                FilterReadsError::CliError(
                                    "Must supply either a configuration file or a reference file to run NEAT!".to_string()
                                )
                            )
                        )
                    }
                    info!("Running filter-reads to filter rneat-generated read data.");
                    let result = filter_reads::main(&file);
                    match result {
                        Err(error) => return Err(NeatErrors::FilterReadsError(error)),
                        Ok(()) => info!("rneat filter-reads completed succesfully"),
                    }
                }
            }
        },
        Some(("gen-mut-model", _)) => {
            if let Some(("gen-mut-model", cmd)) = subcommand {
				info!("Running rneat filter-reads");
                // extract the flags and values
                if cmd.contains_id("configuration_yaml") {
                    let file = cmd.get_one::<PathBuf>("configuration_yaml")
                            .expect("Must provide a path with configuration-yaml")
                            .to_path_buf();

                    if !file.is_file() {
                        return Err(
                            NeatErrors::FilterReadsError(
                                FilterReadsError::CliError(
                                    "Must supply either a configuration file or a reference file to run NEAT!".to_string()
                                )
                            )
                        )
                    }
                    info!("Running gen-mut-model to generate a model for use in rneat.");
                    let result = gen_mut_model::main(&file);
                    match result {
                        Err(error) => return Err(NeatErrors::GenMutModelError(error)),
                        Ok(()) => info!("rneat gen-mut-model completed succesfully"),
                    }
                }
            }
        },
        _ => unreachable!("Parser should ensure no other subcommand choice is made.")
    }

    let elapsed_time = time::Instant::now() - start;
    info!("Processing finished in {} milliseconds", elapsed_time.as_millis());
    if elapsed_time.as_millis() < 1000 {
        info!("Processing finished in {} milliseconds", elapsed_time.as_millis());
    } else if elapsed_time.as_secs() > 300 {
        info!("Processing finished in {} minutes", ((elapsed_time.as_secs() as f64)/60.0))
    } else {
        info!("Processing finished in {} seconds", elapsed_time.as_secs());
    }
    Ok(())
}
