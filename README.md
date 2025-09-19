# The rusty-neat project
Weclome to `rneat`, a Rust port of NEAT (https://github.com/ncsa/neat), a genetic simulation program that creates next-gen-looking fastq files along with an accompanying golden vcf file showing all of the variants inserted by NEAT. It also generates "noise" in the form of sequencing errors as it is writing out files. These features can help you hone in your alignment software to your data. 

As of now, `rneat` is still under development. This is our 1.3.0 version, but it is still in beta. There are several features yet to be completed. One of which is the ability to read data from real data, to hone NEAT's statistics. There are other features we hope to add too. Check the changelog for more information. 

# How to use `rneat`

Download the executable in the release (current version 1.3.0).

```
$ rneat --help
Usage: rneat [OPTIONS] [SUB-COMMAND]

SUB-COMMANDS:
  gen-reads     Generates reads for an input dataset
  filter-reads  Filters the output of gen-reads
  help          Print this message or the help of the given subcommand(s)

Options:
      --log-level [<log_level>]  Sets the log level for the display log. [default: trace] [possible values: trace, debug, info, warn, error, off]
      --log-dest <log_dest>      Sets the log destination (full path with full filename) for the written log
  -h, --help                     Print help
```
To check options for a subcommand:
```
$ rneat gen-reads --help
Generates reads for an input dataset

Usage: rneat gen-reads [OPTIONS]

Options:
  -c, --configuration-yaml <configuration_yaml>  Path to configuration file.
  -h, --help                                     Print help
```
```
$ rneat filter-reads --help
Filters the output of gen-reads

Usage: rneat filter-reads [OPTIONS] --bed-file <bed_file> --output-file <output_file>

Options:
  -b, --bed-file <bed_file>              Path to bed file containing desired regions.
  -f, --file-to-filter <file_to_filter>  Path to fastq_r1 file containing reads to filter.
  -o, --output-file <output_file>        File (including path) to output file where to write files.
  -h, --help                             Print help
```

Use the help menu to see the available options and leave an issue if you find something bad happening. This data is not currently considered usable for anything requiring real rigor, but this is the first iteration toward a final product. 

To compile and run `rneat` yourself, you will need the Rust environment (https://www.rust-lang.org/tools/install), with cargo. You will also need git installed for your operating system. You will then need to git clone and cd into the repo directory. From your home directory in Linux the process might look something like:

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
~/rusty-neat/$ ./target/release/rneat -h
```

`rneat` uses a configuration file to read values it needs for the run. Most of these features should be active, but there is currently no way to generate your own data, so stick with default models for now. Also, there is not yet a way to create bams. A command line execution might look like this:

```angular2html
~/rusty-neat/$ ./target/release/rneat -C /path/to/filled/in/config.yml
```
If you record the output in the logs of Seed string to regenerate these exact results: XXXXXXX, you should be able to use that string as input with rng_seed and reproduce your results (untested as of yet).

Fastq Output
============
the fastq output will have a key name that identifies the block where teh read was drawn from. This should allow you to asseses how well the aligner funcution. The name will have the format `<contig short name>_<fragment_start>_<fragment_end>`

```angular2html
@neat_generated_Chromosome_0000000000_0000131072_1:1
CTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTTAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGC
+
>AC7<GDEGGGGEFGA<GFCGG;GGGGGF>GEGGEGGGFGFEFGCEGGGGGGCG:AEFGFFGG>FG;GDGA9$GGAGF=GFG=EFFCGGGFGGGGGC$BFGEFFAGGG9F7E@>?GFGGGG>EBFGFDG)DGC6DEDFA2EG:EGG%FFB
```
For example, comes from the chromosome in the reference (which was simply named "Chromosome") between 0 and 131072, which is the size 
of a block read by the fasta reader. An improvement may be to give the exact coordinates of the read. We will continue to think about bam creation, but this may be a good method as well to see how well the alignment worked. Perhaps we could write a processing script that took the output of GATK alignment and checked to see how well the reads were aligned relative to where they came from. That might be easier than writing bam tools.

Reading Bed Data
================
`rneat` can now read bed data and filter the reads based on regions specified. This comes with some caveats. First is that rust can't handle any header lines or non-header rows in the bed file. Second, the names must match what is in the fasta. Third, `rneat` can only filter, and does not use the rest of the bed. If this is an issue please raise an issue and we will see about adding more features.

Bed data will have some challenges. For example, if the contig names in the bed file don't match the assumed contig name from the fasta file, how will rust be able to know which is what? To make things work, we are going to assume, for now, that the bed file contig names match the names as derived by `rneat`. To determine this:

First, in a linux terminal, use your input fasta to find the names of the contigs like this:

```angular2html
$ grep "^>" input.fasta
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
>NC_001134.8 Saccharomyces cerevisiae S288C chromosome II, complete sequence
>NC_001135.5 Saccharomyces cerevisiae S288C chromosome III, complete sequence
>etc
```
This will give you the fill fasta name. `rneat` determines the short name of the input fasta by skipping the initial '>' character, then taking the string up to the first delimiter (space or '|'). In the above example, the short names are "NC_001133.9", "NC_001134.8", etc. 

Next, we check the bed. This example awk command will look for unique values in the first column of your bed file and print out what it finds. 

```angular2html
$ awk '!seen[$1]++ {print $1}' file.bed
1
2
3
etc
```
In this case, the bed file uses a simple numbering scheme to number the chromosomes. We can see a mapping with the roman numeral chromosomes in the name, but it's tricky to program. So, the following command will allow you, the user, to map these chromosomes and thus instruct NEAT. Yours will vary based on the inputs.

```angular2html
awk -F'\t' '{$1=($1=="1"?"NC_001133.9":$1); $1=($1=="2"?"NC_001134.8":$1); print}' OFS='\t' file.bed > file_renamed.bed
```
You will have to tailor this command to your dataset, but this will allow you to map the names from the bed to the fasta in a way `rneat` will understand.

Filtering Your Data
====================

There are two ways to filter with `rneat`. The first is to add a path to a bed in the config with the keyword `filter_output`

In configuration:
```
filter-output: /path/to/my.bed
```

Now, calling `rneat gen-reads`:
```
$ rneat gen-reads -c simple.yml
```

Will produce normal fastq and vcfs, alongside filtered fastq and vcfs, all in gzip format.

The second way is directly on existing files:
```
$ rneat gen-reads -b /path/to/my.bed -f myfastq_r1.fastq.gz -o myfastq_r1_filtered.fastq.gz
```
Will filter myfastq_r1.fastq into myfastq_r1_filtered.fastq.gz. This should work seamlessly for both fastq and vcf, whether it is gzipped or not. 

IMPORTANT: This feature only works on fastq files generated by `rneat`. It is untested on generic VCF files, but should work. The reason it only works with `rneat`-generated fastq files is that `rneat` uses a naming scheme that the parser can use to identify where the read is located (enabling filtering), whereas an average fastq will not have that info.

Try it out and let us know if you run into any issues!

We will continue to improve this process in the future. If you have suggestions for parsing names to facilitate this, let us know in the Issues section!

That's it so far! Test out the features and let us know what you think!