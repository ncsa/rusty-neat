//! This contains only one function at the moment, a general opener and reader
use std::io::{BufReader, BufRead, Result, Lines, self, Write};
use std::fs::File;
use std::path::Path;

pub fn read_lines<T: AsRef<Path> + ?Sized>(filename: &T) -> Result<Lines<BufReader<File>>> {
    // This opens file and creates a buffer to read lines
    let file = File::open(filename)?;
    Ok(BufReader::new(file).lines())
}

pub fn create_output_file<T: AsRef<Path>> (filename: &T, overwrite_file: bool) -> Result<File> {
    if Path::new(filename.as_ref()).is_file() && !overwrite_file {
        // The file already exists and we're in non overwrite mode
        panic!("Attempting to overwrite an existing file: {}", filename.as_ref().display())
    } else {
        File::create(&filename)
    }
}

pub fn append_to_file<T: AsRef<Path>>(filename: &T) -> Result<File> {
    // if the file doesn't exist, we'll create it. If it does, we'll append to it
    if Path::new(&filename.as_ref()).is_file() {
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
