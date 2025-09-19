//! Various file tools needed throughout the code.
use log::*;
use std::path::PathBuf;
use std::{fs, io};

pub fn check_parent(filename: &PathBuf, create: bool) -> io::Result<&PathBuf> {
    // checks that the parent dir exists and then if so creates the Path object open
    // and ready to write
    let file_path = filename;
    let parent = file_path
        .parent()
        .unwrap()
        .to_path_buf();
    if !parent.exists() && create {
        check_create_dir(&parent);
    } else if !parent.exists() {
        error!("Directory {} does not exist!", parent.to_str().unwrap());
        panic!("Error creating file log!");
    }
    Ok(file_path)
}

pub fn check_create_dir(dir_to_check: &PathBuf) {
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
        let filename = PathBuf::from("test_data/H1N1.fa");
        check_parent(&filename, false).unwrap();
    }

    #[test]
    fn test_check_parent_fail() {
        let filename = PathBuf::from("fake/test.fa");
        assert!(!(PathBuf::from("fake").is_dir()));
        check_parent(&filename, true).unwrap();
        assert!(PathBuf::from("fake").is_dir());
        fs::remove_dir("fake").unwrap()
    }
}
