//! This contains only one function at the moment, a general opener and reader
use std::io::{BufReader, BufRead, Result, Lines, self, Write};
use std::fs::File;
use std::path::PathBuf;

use flate2::read::GzDecoder;

pub fn read_lines(filename: &PathBuf) -> Result<Lines<BufReader<File>>> {
    // This opens file and creates a buffer to read lines
    let file = File::open(filename)?;
    Ok(BufReader::new(file).lines())
}

pub fn read_gzip_lines(filename: &PathBuf) -> Result<Lines<BufReader<GzDecoder<File>>>> {
    // Reads a file from gzipped format into lines
    let file = File::open(filename)?;
    Ok(BufReader::new(GzDecoder::new(file)).lines())
}

pub fn open_safe(filename: &PathBuf) -> Result<BufReader<File>> {
    // This opens a file for reading without compression
    let file = File::open(filename)?;
    Ok(BufReader::new(file))
}

pub fn create_output_file (filename: &PathBuf, overwrite_file: bool) -> Result<File> {
    if filename.is_file() && !overwrite_file {
        // The file already exists and we're in non overwrite mode
        panic!("Attempting to overwrite an existing file: {}", filename.display())
    } else {
        File::create(&filename)
    }
}

pub fn append_to_file(filename: &PathBuf) -> Result<File> {
    // if the file doesn't exist, we'll create it. If it does, we'll append to it
    if filename.is_file() {
        File::options()
            .write(true)
            .append(true)
            .open(&filename)
    } else {
        File::create(&filename)
    }
}

pub struct VectorBuffer {
    // The purpose of this struct is to mock a buffer so we can run single-ended read generation
    buffer: Vec<u8>,
}

impl VectorBuffer {
    pub fn new() -> VectorBuffer {
        VectorBuffer { buffer: Vec::new() }
    }
}

impl Write for VectorBuffer {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.buffer.extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        // In this case, flushing does nothing because
        // VectorBuffer only simulates buffering to stdout or other medium
        Ok(())
    }
}
