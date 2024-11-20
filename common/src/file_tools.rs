// Various file tools needed throughout the code.
use log::warn;
use std::fs::File;
use std::io::{BufReader, BufRead, Error, Lines};
use std::path::Path;
use std::{env, fs, io};
use std::ffi::OsString;

pub fn read_lines(filename: &str) -> io::Result<Lines<BufReader<File>>> {
    // This creates a buffer to read lines
    let file = File::open(filename)?;
    Ok(BufReader::new(file).lines())
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

pub fn check_parent(filename: &str, create: bool) -> io::Result<&Path> {
    // checks that the parent dir exists and then if so creates the Path object open
    // and ready to write
    let file_path = Path::new(filename);
    let parent = file_path.parent().unwrap();
    if !parent.exists() && create {
        check_create_dir(parent);
    } else if !parent.exists() {
        println!("{}", env::current_dir().unwrap().to_str().unwrap());
        panic!("Directory {} does not exist!", parent.to_str().unwrap());
    }
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
        check_parent(filename, false).unwrap();
    }

    #[test]
    fn test_check_parent_fail() {
        let filename = "fake/test.fa";
        assert!(!Path::new("fake").is_dir());
        check_parent(filename, true).unwrap();
        assert!(Path::new("fake").is_dir());
        fs::remove_dir("fake").unwrap()
    }
}
