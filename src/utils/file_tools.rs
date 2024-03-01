use std::fs::File;
use std::io;
use std::io::BufRead;

pub fn read_lines(filename: &str) -> io::Result<io::Lines<io::BufReader<File>>> {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}