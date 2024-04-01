// Various file tools needed throughout the code.
use log::warn;
use std::fs::File;
use std::io::{BufRead, Error};
use std::path::Path;
use std::{fs, io};

pub fn read_lines(filename: &str) -> io::Result<io::Lines<io::BufReader<File>>> {
    // This creates a buffer to read lines
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn open_file(mut filename: &mut str, overwrite_file: bool) -> Result<File, Error> {
    if overwrite_file && Path::new(filename).exists() {
        File::options().create(true).write(true).open(&mut filename)
    } else {
        File::options()
            .create_new(true)
            .append(true)
            .open(&mut filename)
    }
}

pub fn check_parent(filename: &str) -> io::Result<&Path> {
    // checks that the parent dir exists and then if so creates the Path object open
    // and ready to write
    let file_path = Path::new(filename);
    let parent = file_path.parent().unwrap();
    if !parent.exists() {
        check_create_dir(parent);
    };
    Ok(file_path)
}

pub fn check_create_dir(dir_to_check: &Path) {
    if !dir_to_check.is_dir() {
        warn!("Directory not found, creating: {:?}", dir_to_check);
        fs::create_dir(dir_to_check).expect("Error creating the directory");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_parent() {
        let filename = "test_data/H1N1.fa";
        check_parent(filename).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_check_parent_fail() {
        let filename = "fake/test.fa";
        check_parent(filename).unwrap();
    }
}
