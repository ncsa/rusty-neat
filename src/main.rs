#![allow(dead_code)]
mod logger;
use rusty_neat::{Config};
use std::{env, process};



fn main() {
    logger::init();
    logger::info!("Begin processing");

    // Argument parser
    let args: Vec<String> = env::args().collect();
    let config = Config::build(&args).unwrap_or_else(|err| {
        println!("Problem parsing arguments: {err}");
        process::exit(1)
    });

    println!("In file {}", config.fasta_path);

    if let Err(e) = rusty_neat::run(config) {
        println!("Application error: {e}");
        process::exit(1);
    }
}

