extern crate clap;
extern crate itertools;
extern crate log;
extern crate serde;
extern crate serde_json;
extern crate simplelog;
extern crate simple_rng;

pub mod filter_reads;
pub mod gen_reads;

use std::env;
use same_file::is_same_file;
use common::{self, file_tools::file_io::create_output_file};

use common::file_tools::folder_tools::{check_create_dir, check_parent};
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

use crate::{
    filter_reads::{
        FilterReadsError,
    }, 
    gen_reads::{
        errors::GenerateReadsErrors,
    },
};

#[derive(Error, Debug)]
pub enum NeatErrors {
    #[error("Error while generating read dataset")]
    GenerateReadsError(#[from] GenerateReadsErrors),
    #[error("Error while filtering reads with bed file")]
    FilterReadsError(#[from] FilterReadsError),
}

fn neat_commands() -> [Command; 2] {
    [
        Command::new("gen-reads")
            .about("Generates reads for an input dataset")
            .arg_required_else_help(true)
            .arg(
                Arg::new("configuration_yaml")
                    .long("configuration-yaml")
                    .short('c')
                    .help("Path to configuration file.")
                    .exclusive(true)
                    .action(ArgAction::Set)
                    .default_missing_value("")
                    .value_parser(value_parser!(PathBuf))
            ),
        Command::new("filter-reads")
            .about("Filters the output of gen-reads")
            .arg_required_else_help(true)
            .arg(
                Arg::new("bed_file")
                    .long("bed-file")
                    .short('b')
                    .help("Path to bed file containing desired regions.")
                    .action(ArgAction::Set)
                    .default_missing_value("")
                    .value_parser(value_parser!(PathBuf))
                    .required(true)
            )
            .arg(
                Arg::new("file_to_filter")
                    .long("file-to-filter")
                    .short('f')
                    .help("Path to fastq_r1 file containing reads to filter.")
                    .action(ArgAction::Set)
                    .default_missing_value("")
                    .value_parser(value_parser!(String))
            )
            .arg(
                Arg::new("output_file")
                    .long("output-file")
                    .short('o')
                    .help("File (including path) to output file where to write files.")
                    .action(ArgAction::Set)
                    .default_missing_value("")
                    .value_parser(value_parser!(PathBuf))
                    .required(true)
            ),
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
    info!("////////////// Welcome to rusty-neat! \\\\\\\\\\\\\\\\\\\\\\\\");

    match subcommand {
        Some(("gen-reads", _)) => {
            if let Some(("gen-reads", cmd)) = subcommand {
                if cmd.contains_id("configuration_yaml") {
                    let file = cmd.get_one::<PathBuf>("configuration_yaml")
                            .expect("Must provide a path with configuration-yaml")
                            .to_path_buf();
                    info!("Using configuration file: {:?}", file);
                    info!("Running Generate Reads.");
                    let result = gen_reads::main(&file);
                    match result {
                        Err(error) => return Err(NeatErrors::GenerateReadsError(error)),
                        _ => info!("rneat gen-reads completed successfully"),
                    }
                } 
            }
        },
        Some(("filter-reads", _)) => {
            if let Some(("filter-reads", matches)) = subcommand {
                // extract the flags and values
                let mut is_fastq = true;
                let mut is_gzip = false;
                let bed_path: PathBuf = {
                    if let Some(bed_file) = matches.get_one::<PathBuf>("bed_file") {
                        info!("filter-reads bed file: {:?}", bed_file);
                        if !bed_file.is_file() {
                            return Err(NeatErrors::FilterReadsError(FilterReadsError::FileNotFound(format!("{:?}", bed_file))))
                        }
                        bed_file.to_path_buf()
                    } else {
                        return Err(NeatErrors::FilterReadsError(FilterReadsError::FileNotFound("No bed file received.".to_string())))
                    }
                };
                let file_to_filter: PathBuf = {
                    if let Some(input_file) = matches.get_one::<String>("file_to_filter") {
                        info!("filter-reads file to filter: {:?}", input_file);
                        let file_path = PathBuf::from(input_file);
                        let extension = file_path.extension();
                        match extension {
                            Some(ext) => {
                                if ext != "gz" {
                                    is_gzip = true;
                                    let shorter_fn = file_path.with_extension("");
                                    let extension2 = shorter_fn.extension();
                                    match extension2 {
                                        Some(ext2) => {
                                            if ext2 == "vcf" {
                                                is_fastq = false;
                                            }
                                        },
                                        None => {
                                            panic!("No file extension for gzipped file. Filetype not recognized")
                                        },
                                    }
                                } else if ext == "vcf" {
                                    is_fastq = false;
                                }
                            },
                            None => return Err(
                                NeatErrors::FilterReadsError(
                                    FilterReadsError::MalformedFileName(
                                        "Input file has no valid extensions.".to_string()
                                    )
                                )
                            ),
                        }
                        PathBuf::from(input_file)
                    } else {
                        return Err(
                            NeatErrors::FilterReadsError(
                                FilterReadsError::FileNotFound(
                                    "Input file not received.".to_string()
                                )
                            )
                        )
                    }
                };
                let output_file: PathBuf = {
                    if let Some(output_opt) = matches.get_one::<PathBuf>("output_file") {
                        info!("filter-reads output dir: {:?}", output_opt);
                        check_create_dir(output_opt);
                        output_opt.to_path_buf()
                    } else {
                        return Err(
                            NeatErrors::FilterReadsError(
                                FilterReadsError::FileNotFound(
                                    "No output file received.".to_string()
                                )
                            )
                        )
                    }
                };
                if is_same_file(&file_to_filter, &output_file).unwrap() {
                    return Err(
                        NeatErrors::FilterReadsError(
                            FilterReadsError::FileNotFound(
                                "Output and input files are the same. Cannot read and write the same file".to_string()
                            )
                        )
                    )
                }
                info!("Running filter reads!");
                let result: Result<(), FilterReadsError> = filter_reads::main(
                    bed_path,
                    file_to_filter,
                    is_gzip,
                    is_fastq,
                    output_file,
                );
                match result {
                    Err(error) => return Err(NeatErrors::FilterReadsError(error)),
                    _ => {
                        info!("rneat filter-reads completed successfully!");
                    },
                }
            }
        },
        _ => unreachable!("Parser should ensure no other subcommand choice is made.")
    }

    Ok(())
}

