# RustyNeat v1.0.0

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
========
9/9/2025

# rusty-neat v1.0.0
- It's been a couple years since I had an update to rusty-neat. This new vrison is version 1.0.0! It is not yet fully functional, as it is missing the BAM creation, though it does work on some test data: H1N1.fa and ecoli.fa. It's very fast for H1N1, but that's not too surprising. It's slow right now on fastq creation. I already know what I need to do to fix that part, I put an [issue](https://github.com/ncsa/rusty-neat/issues/65) up here that I am working on in develop. But for now, rusty-neat functions and can be executed on some arbitrary data. Will update the readme as well. Once I get the fastq generation idea to work (will require changes to fastq_tools and runner scripts), next step will be parallelizing the main contig loop in runner. That should give us some speedup. Also, we went from using way too much memory to using almost no memory. I'm not sure which is better. It would be nice to find a balance. I still think that a database could work. I'm wondering, maybe we index on position. Just use like a json-based database: {"contig": "chr1", "sequence": {"1":A, "2": C, etc}} We'd want something that had fast retrieval for writing. Something to investigate, maybe.
- The code is significantly different now, I read in all the old models from the original NEAT2.0, and stored them as default models. I used a lot of the same algorithms that they used in NEAT 2.0, figuring it was closer to rust than NEAT 3.X. 
- rusty-neat writes out a ton of files. It should clean them up afterward, but I'm not sure it does yet, may need to investigate if there is an added command required. But it stores them in temp and linux cleans that up. Basically, we are splitting up the input fasta into smaller chunks (tbd on if there is an optimal size) then writing those to file for later retrival. Need to test various chunk sizes. I am basically doing blocked gzip-type things, but with some twists. If someone really understood bgzip format, it could probably be used directly. I'm not there yet.
- My runtimes running on a WSL2 Ubuntu setup with the latest rust installed: 311 milliseconds for H1N1.fa, 13 minutes for ecoli. It's all the IO during that step. I knew it was ugly as I was writing it but it's a rough draft. Obviously all the repeat code in fastq_tools is not ideal.
- The output files, so far have looked good. I realized I haven't tested paired-ended yet.
- I decided to map out the N-regions then overwrite them with random valid bases. The main reason being that it kept dipping into N regions and trying to use them for snp creation and that isn't going to work. NEAT 2.0 had a data structure that kept track of how long N-regions were, then would skip them if they were longer than a read length or overwrite them if they were short. That might work. Potential room for improvement, but it's fast so far.
- Eliminated bam creation for now. The option is in there, but it warns you that it doesn't work yet. Once fastqs are sorted and maybe after parallelism, it'll be next up.

=========
9/11/2025

# rusty-neat v1.1.0
- We have refactored the fastq code. Now a temporary fastq is written to file per block read in by the fasta reader. This should allow us to parallelize the temporary fastqs. Then it's simply a matter of concatenating them together. Works great! I got through single-ended E coli depth 10 in 36 seconds. I got through paired-ended E coli in 78 seconds. 
- Caveats: the fastqs are sorted by read start. Included in each name is the exact location where the sequence was drawn from. Now, when you run an alignment, the reads will be named in such a way that will have <chromosome name>_start-pos_length in the name, so you can shuffle the reads and run an alignment, and easily check to see if your alignment is in the correct contig and location.