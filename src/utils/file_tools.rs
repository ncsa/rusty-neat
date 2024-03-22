// Various file tools needed throughout the code.

use std::fs::File;
use std::io;
use std::io::{BufRead, Error};
use std::path::Path;

pub fn read_lines(filename: &str) -> io::Result<io::Lines<io::BufReader<File>>> {
    // This creates a buffer to read lines
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn open_file(mut filename: &mut str, overwrite_file: bool) -> Result<File, Error> {
    if overwrite_file {
        File::options().create(true).write(true).open(&mut filename)
    } else {
        File::options().create_new(true).append(true).open(&mut filename)
    }
}

pub fn check_parent(filename: &str) -> io::Result<&Path> {
    // checks that the parent dir exists and then if so creates the Path object open
    // and ready to write
    let file_path = Path::new(filename);
    assert!(file_path.parent().unwrap().exists());
    Ok(file_path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_parent() {
        let filename = "data/H1N1.fa";
        check_parent(filename).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_check_parent_fail() {
        let filename = "fake/test.fa";
        check_parent(filename).unwrap();
    }
}