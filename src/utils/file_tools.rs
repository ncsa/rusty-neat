use std::fs::{File, OpenOptions};
use std::io;
use std::io::{BufRead, Error};
use std::path::Path;
use log::error;

pub fn read_lines(filename: &str) -> io::Result<io::Lines<io::BufReader<File>>> {
    // This creates a buffer to read lines
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn open_file(mut filename: &mut str, overwrite_file: bool) -> File {
    if overwrite_file {
        File::options().create(true).write(true).open(&mut filename).unwrap()
    } else {
        File::options().create_new(true).append(true).open(&mut filename).unwrap()
    }
}

pub fn check_parent_and_create(filename: &str) -> io::Result<&Path> {
    // checks that the parent dir exists and then if so creates the file object open
    // and ready to write
    let file_path = Path::new(filename);
    if !file_path.parent().unwrap().exists() {
        error!("Path to log file not found!");
        assert!(file_path.parent().unwrap().exists())
    }
    Ok(file_path)
}