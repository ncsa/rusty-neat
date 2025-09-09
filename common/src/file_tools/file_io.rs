//! This contains only one function at the moment, a general opener and reader
use std::io::{BufReader, BufRead, Result, Lines};
use std::fs::File;
use std::path::Path;

pub fn read_lines<T: AsRef<Path> + ?Sized>(filename: &T) -> Result<Lines<BufReader<File>>> {
    // This opens file and creates a buffer to read lines
    let file = File::open(filename)?;
    Ok(BufReader::new(file).lines())
}

pub fn create_output_file<T: AsRef<Path>> (filename: &T, overwrite_file: bool) -> Result<File> {
    if Path::new(filename.as_ref()).exists() && !overwrite_file {
        // The file already exists and we're in non overwrite mode
        panic!("Attempting to overwrite an existing file: {}", filename.as_ref().display())
    } else {
        File::options()
            .create(true)
            .write(true)
            .open(&filename)
    }
}

pub fn append_to_file(filename: &str) -> Result<File> {
    // if the file doesn't exist, we'll create it. If it does, we'll append to it
    if Path::new(&filename).exists() {
        File::options()
            .write(true)
            .append(true)
            .open(&filename)
    } else {
        File::options()
            .create(true)
            .write(true)
            .open(&filename)
    }
}