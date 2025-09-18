extern crate clap;
extern crate itertools;
extern crate log;
extern crate serde;
extern crate serde_json;
extern crate simplelog;
extern crate simple_rng;

pub mod filter_reads;
pub mod gen_reads;

use std::{env, fs::File};
use common;

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
                    .exclusive(true)
                    .action(ArgAction::Set)
                    .default_missing_value("")
                    .value_parser(value_parser!(PathBuf))
                    .required(true)
            )
            .arg(
                Arg::new("fastq1_file")
                    .long("fastq1-file")
                    .short('1')
                    .help("Path to fastq_r1 file containing reads to filter.")
                    .exclusive(true)
                    .action(ArgAction::Set)
                    .default_missing_value("")
                    .value_parser(value_parser!(PathBuf))
            )
            .arg(
                Arg::new("fastq2_file")
                    .long("fastq2-file")
                    .short('2')
                    .help("Path to fastq_r2 file containing reads to filter.")
                    .exclusive(true)
                    .action(ArgAction::Set)
                    .default_missing_value("")
                    .value_parser(value_parser!(PathBuf))
            )
            .arg(
                Arg::new("vcf_file")
                    .long("vcf-file")
                    .short('v')
                    .help("Path to vcf file containing variants to filter.")
                    .exclusive(true)
                    .action(ArgAction::Set)
                    .default_missing_value("")
                    .value_parser(value_parser!(PathBuf))
            )
            .arg(
                Arg::new("output_dir")
                    .long("output-dir")
                    .short('o')
                    .help("Path to output dir where to write files.")
                    .exclusive(true)
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
                        .exclusive(true)
                        .action(ArgAction::Set)
                        .num_args(0..=1)
                        .require_equals(true)
                        .value_parser(["trace", "debug", "info", "warn", "error", "off"])
                        .default_value("trace")
                        .default_missing_value("trace")
                )
                .arg(
                    Arg::new("log_dest")
                        .long("log-dest")
                        .help("Sets the log destination (full path with full filename) for the written log")
                        .exclusive(true)
                        .action(ArgAction::Set)
                        .default_missing_value("")
                )
                .subcommands(neat_commands())
        )
        .subcommands(neat_commands());
        

    let matches = cmd.get_matches();
    let mut subcommand = matches.subcommand();
    let mut level_filter = String::new();
    let mut log_dest = env::current_dir().unwrap().with_file_name(".neat.log");
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
    let log_destination = check_parent(&log_dest, true).unwrap();
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
            File::create(log_destination).unwrap(),
        ),
    ])
    .unwrap();

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
                        _ => {},
                    }
                } 
            }
        },
        Some(("filter-reads", _)) => {
            if let Some(("filter_reads", matches)) = subcommand {
                // extract the flags and values
                let bed_path: PathBuf = {
                    if let Some(bed_file) = matches.get_one::<String>("bed_file") {
                        info!("filter-reads bed file: {bed_file}");
                        PathBuf::from(bed_file)
                    } else {
                        return Err(NeatErrors::FilterReadsError(FilterReadsError::FileNotFound("No bed file received.".to_string())))
                    }
                };
                let fastq1_file: Option<PathBuf> = {
                    if let Some(fastq_file) = matches.get_one::<String>("fastq1_file") {
                        info!("filter-reads fastq_r1 file: {fastq_file}");
                        Some(PathBuf::from(fastq_file))
                    } else {
                        None
                    }
                };
                let fastq2_file: Option<PathBuf> = {
                    if let Some(fastq_file) = matches.get_one::<String>("fastq2_file") {
                        info!("filter-reads fastq_r2 file: {fastq_file}");
                        Some(PathBuf::from(fastq_file))
                    } else {
                        None
                    }
                };
                let vcf_file: Option<PathBuf> = {
                    if let Some(vcf_file) = matches.get_one::<String>("vcf_file") {
                        info!("filter-reads fastq_r2 file: {vcf_file}");
                        Some(PathBuf::from(vcf_file))
                    } else {
                        None
                    }
                };
                if fastq1_file.is_none() && fastq2_file.is_none() && vcf_file.is_none() {
                    warn!("Received no input files to filter. Exiting.");
                    debug!("If this exit is unexpected, check input flags and paths for files. Use full paths when possible.");
                    // Hard to say if this is an error or bad input.
                    return Ok(())
                }
                let output_dir: PathBuf = {
                    if let Some(otpt_dir) = matches.get_one::<String>("output_dir") {
                        info!("filter-reads output dir: {otpt_dir}");
                        PathBuf::from(otpt_dir)
                    } else {
                        return Err(NeatErrors::FilterReadsError(FilterReadsError::FileNotFound("No output dir received.".to_string())))
                    }
                };
                info!("Running filter reads!");
                let result: Result<(), FilterReadsError> = filter_reads::main(
                    bed_path,
                    fastq1_file,
                    fastq2_file,
                    vcf_file,
                    output_dir,
                );
                match result {
                    Err(error) => return Err(NeatErrors::FilterReadsError(error)),
                    _ => {},
                }
            }
        },
        _ => unreachable!("Parser should ensure no other subcommand choice is made.")
    }

    Ok(())
}

