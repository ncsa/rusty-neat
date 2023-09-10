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