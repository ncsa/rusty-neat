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

# rusty-neat v0.1

Goal: To release a version of NEAT in rust capable of reading an arbitrary fasta file, 
and outputting a version of that file with 1% of the bases mutated as SNPs.
=======
9/11/2023
- Found a sqlite crate for rust. This might be a good way to store the reference information for our genome. We can have one row per chromosome. Name, info, sequence. Read the data in once, store it for long term retrieval. Website: https://docs.rs/sqlite/latest/sqlite/

# rusty-neat v0.1.0

Moved the code from a brach under NEAT (https://github.com/ncsa/neat) to its own repo https://github.com/ncsa/rusty-neat. This project will be considered a separate, related project. We will use the same underlying design ideas as NEAT, but it will essentially be different and won't be guaranteed in any way to be able to reproduce what the original NEAT did (though we aim to get it as close as possible, the differences in Python and Rust would make byte-by-byte compatability very difficult). I have implemented an idea I prototyped in a private repo (https://github.com/joshfactorial/basic_neat), of reading in the file in what I hope is a straightforward way. I was looking for stuff like running out of memory and speed. It can read in the `E. coli` fasta file in about 30 seconds on my laptop, which is pretty fast, compared to Python. What's more is it can store the whole thing in memory, due to the efficient way Rust uses memory. Instead of copying the entire data set every time in memory, quickly eating it up, Rust forces you to decide when to pass a read-only version, copy the dataset, or use up the dataset once and lose it forever (all valid choices depending on context). So I'm feeling hopeful that this will work. I am publishing an Alpha called v0.1.0 for people to test, if they feel so inclined. See the README for more details.