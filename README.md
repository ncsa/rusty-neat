# The rusty-neat project
Rusty-neat is a Rust port of NEAT (https://github.com/ncsa/neat), a genetic simulation program that creates next-gen-looking fastq files along with an accompanying golden vcf file showing all of the variants inserted by NEAT. It also generates "noise" in the form of sequencing errors as it is writing out files. These features can help you hone in your alignment software to your data. 

As of now, rusty-neat is still under development. This is our 1.0.0 version, but it is still in alpha. There are several features yet to be completed. One of which is the ability to read data from real data, to hone NEAT's statistics. There are other features we hope to add too. Check the changelog for more information. 

# How to use rusty-neat

Download the executable in the release (current version 1.0.0).

```
./rusty-neat -h
```

displays help

Use the help menu to see the available options and leave an issue if you find something bad happening. This data is not currently considered usable for anything requiring real rigor, but this is the first iteration toward a final product. 

To compile and run rusty-neat yourself, you will need the Rust environment (https://www.rust-lang.org/tools/install), with cargo. You will also need git installed for your operating system. You will then need to git clone and cd into the repo directory. From your home directory in Linux the process might look something like:

```angular2html
~/$ git clone git@github.com:ncsa/rusty-neat.git
~/$ cd rusty-neat
~/rusty-neat/$
```

For Windows users, you're probably going to want to use WSL, but go ahead and test it if there's a Rust for windows. (Post help requests in the Issues tab if you need assistance or have suggestions).

Once in the repo, you can build the program either in debug (default) or release mode. The main difference is how much info it gives you if there is an error. Release mode also has some optimizations to run it faster.

```angular2html
~/rusty-neat/$ cargo build --release
```
If you prefer to run the package directly without using the binary, you can also use
```angular2html
~/rusty-neat/$ cargo run
```
Rust will download any required packages. Compiling Rust code is the slowest part of the process. The final binary will be built and the program run immediately after in the second case. To run the program manually, from the repo main dir, run

```angular2html
~/rusty-neat/$ ./target/debug/rusty-neat -h
```

You can either run rusty-neat as a cli or using a configuration file. Most of these features should be active, but there is currently no way to generate your own data, so stick with default models for now. Also, there is not yet a way to create bams. A command line execution might look like this:

```angular2html
~/rusty-neat/$ ./target/debug/rusty-neat --log-level Info -r /path/to/my_bacterium.fa 
```
or for a configuration file:
```angular2html
~/rusty-neat/$ ./target/debug/rusty-neat -C /path/to/filled/in/config.yml
```
If you record the output in the logs of Seed string to regenerate these exact results: XXXXXXX, you should be able to use that string as input with rng_seed and reproduce your results (untested as of yet).

##Fastq Output
============
the fastq output will have a key name that identifies the block where teh read was drawn from. This should allow you to asseses how well the aligner funcution. The name will have the format <contig short name>_<block start>_<block_end>

```angular2html
@neat_generated_Chromosome_0000000000_0000131072_1:1
CTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTTAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGC
+
>AC7<GDEGGGGEFGA<GFCGG;GGGGGF>GEGGEGGGFGFEFGCEGGGGGGCG:AEFGFFGG>FG;GDGA9$GGAGF=GFG=EFFCGGGFGGGGGC$BFGEFFAGGG9F7E@>?GFGGGG>EBFGFDG)DGC6DEDFA2EG:EGG%FFB
```
For example, comes from the chromosome in the reference (which was simply named "Chromosome") between 0 and 131072, which is the size 
of a block read by the fasta reader. An improvement may be to give the exact coordinates of the read. We will continue to think about bam creation, but this may be a good method as well to see how well the alignment worked. Perhaps we could write a processing script that took the output of GATK alignment and checked to see how well the reads were aligned relative to where they came from. That might be easier than writing bam tools.

That's it so far! Test out the features and let us know what you think!