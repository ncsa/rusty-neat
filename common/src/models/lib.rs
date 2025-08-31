//! This library contains some data and functions that are useful in the other models
//! 
use std;
use std::fs::File;
use std::io::{BufReader, Write};
use gzp::{bgzf, Compression};
use serde::Deserialize;
use crate::file_tools::file_io::create_output_file;

pub fn model_writer(data: String, filename: &str) -> std::io::Result<()> {
    // This will take any serializable model and write it to file.
    let fileout = create_output_file(&filename, false)
        .expect(&format!("Error creating output {}. If file alread exists, remove.", &filename));
    let mut buffer = bgzf::BgzfSyncWriter::new(
        fileout, Compression::fast()
    );
    buffer.write(data.as_bytes())?;
    Ok(())
}

pub fn model_reader<T>(filename: &str) -> Result<T, std::io::Error>
where
    T: for<'de> Deserialize<'de>,
{
    // This will take any serializable model and write it to file.
    let filein = File::open(filename)
        .expect(&format!("Error opening file {}", &filename));
    let reader = BufReader::new(filein);
    let buffer = bgzf::BgzfSyncReader::new(reader);
    let value: T = serde_json::from_reader(buffer).unwrap();
    Ok(value)
}