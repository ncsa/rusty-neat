//! IO for NEAT
use std::fs::{File, OpenOptions, Path, PathBuf};
use std::io::prelude::*;
use std::env;

const BYTES_PER_LINE: usize = 16;

fn is_compressed(file: &PathBuf) -> bool {
    todo!()
}

pub fn read_file() {

    // Code to open a file for appending
    // let f = OpenOptions::new()
    //     .read(true)
    //     .write(true)
    //     .create(true)
    //     .append(true)
    //     .open(path)?;

    let hello = PathBuf::from("/tmp/hello.txt");

    let arg1 = env::args().nth(1);
    let fname = arg1.expect("usage: fview FILENAME");

    let mut f = File::open(&fname).expect("Unable to open file.");
    let mut pos = 0;
    let mut buffer = [0; BYTES_PER_LINE];

    while let Ok(_) = f.read_exact(&mut buffer) {
        print!("[0x{:08x}] ", pos);
        for byte in &buffer {
            match *byte {
                0x00 => print!(".  "),
                0xff => print!("## "),
                _ => print!("{:02x} ", byte),
            }
        }

        println!();
        pos += BYTES_PER_LINE
    }
}

