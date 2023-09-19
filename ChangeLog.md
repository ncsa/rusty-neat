# RustyNeat v0.1

9/5/2023 
- Initialized Cargo package with hello world default script.
- Added function to select a random integer from a list (for selecting positions)

9/9/2012
- Created a first pass program that mutates a hardcoded DNA string at one random non-N position.
- Added a logging function to print out time-stamped log items both to stdout and to file.
- Need to improve error handling all around, taking advantage of Rusts error handling system.
- Design decision: After toying with some options to create a standardized RNG to use throughout the program, it seems the best option is to use an on-the-fly rng per function for now. This eliminates reproducibility between runs, but should enable parallelism of processing.
- Redesign from the ground up. Trying to use KISS as much as possible, and only introduce complexity as it becomes necessary, or as I learn Rust tricks that may help this endeavor.

<<<<<<< HEAD
# rusty_neat v0.2

Goal: To release a version of NEAT in rust capable of reading an arbitrary fasta file, 
and outputting a version of that file with 1% of the bases mutated as SNPs.
=======
9/11/2023
- Found a sqlite crate for rust. This might be a good way to store the reference information for our genome. We can have one row per chromosome. Name, info, sequence. Read the data in once, store it for long term retrieval. Website: https://docs.rs/sqlite/latest/sqlite/
>>>>>>> 853f8db214fcbc9fa6b4b44840b9795cce94d6ac
