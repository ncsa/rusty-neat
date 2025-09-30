# The rusty-neat project
Weclome to `rneat`, a Rust port of NEAT (https://github.com/ncsa/neat), a genetic simulation program that creates next-gen-looking fastq files along with an accompanying golden vcf file showing all of the variants inserted by NEAT. It also generates "noise" in the form of sequencing errors as it is writing out files. These features can help you hone in your alignment software to your data. 

As of now, `rneat` is still under development. We are working toward the goal of the original NEAT of being able to read user data and construct a model based on that data. We do however, have functional reads being generated, with variants tracked in a vcf file. This is our 1.2.0 version.

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

To run filter reads, check the help menu.

```
$ rneat filter-reads --help
Filters the output of gen-reads

Usage: rneat filter-reads --configuration-yaml <configuration_yaml>

Options:
  -c, --configuration-yaml <configuration_yaml>  Path to configuration file.
  -h, --help                                     Print help
```

Use the help menu to see the available options and leave an issue if you find something bad happening.

To compile and run `rneat` yourself, you will need the Rust environment (https://www.rust-lang.org/tools/install), with cargo. You will also need git installed for your operating system. You will then need to git clone and cd into the repo directory. From your home directory in Linux the process might look something like:

```
~/$ git clone git@github.com:ncsa/rusty-neat.git
~/$ cd rusty-neat
~/rusty-neat/$
```

For Windows users, you're probably going to want to use WSL, but go ahead and test it if there's a Rust for windows. (Post help requests in the Issues tab if you need assistance or have suggestions).

Once in the repo, you can build the program either in debug (default) or release mode. The main difference is how much info it gives you if there is an error. Release mode also has some optimizations to run it faster.

```
~/rusty-neat/$ cargo build --release
```

If you prefer to run the package directly without using the binary, you can also use
```
~/rusty-neat/$ cargo run
```

Rust will download any required packages. Compiling Rust code is the slowest part of the process. The final binary will be built and the program run immediately after in the second case. To run the program manually, from the repo main dir, run

```
~/rusty-neat/$ ./target/release/rneat -h
```

`rneat` uses a configuration file to read values it needs for the run. Most of these features should be active, but there is currently no way to generate your own data, so stick with default models for now. Also, there is not yet a way to create bams. A command line execution might look like this:

```
~/rusty-neat/$ ./target/release/rneat -C /path/to/filled/in/config.yml
```

If you record the output in the logs of Seed string to regenerate these exact results: XXXXXXX, you should be able to use that string as input with rng_seed and reproduce your results (needs more thorough testing).

Fastq Output
============
The fastq output will have a key name that identifies the fragment where teh read was drawn from. This should allow you to asseses how well the aligner funcution. The name will have the format `<contig short name>_<fragment_start>_<fragment_end>`

```angular2html
@neat_generated_Chromosome_0000000000_0000353_1:1
CTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTTAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGC
+
>AC7<GDEGGGGEFGA<GFCGG;GGGGGF>GEGGEGGGFGFEFGCEGGGGGGCG:AEFGFFGG>FG;GDGA9$GGAGF=GFG=EFFCGGGFGGGGGC$BFGEFFAGGG9F7E@>?GFGGGG>EBFGFDG)DGC6DEDFA2EG:EGG%FFB
```
For example, comes from the chromosome in the reference (which was simply named "Chromosome") between 0 and 353, which is the size 
of a fragment in the simulated DNA. If there is a pair with this read, it will have the same coordinates, though it started at index 352 instead.

Reading Bed Data
================
`rneat` can now read bed data and filter the reads based on regions specified. This comes with some caveats. First is that rust can't handle any header lines or non-header rows in the bed file. Second, the names must match what is in the fasta. Third, `rneat` can only filter, and does not use the rest of the bed. If this is an issue please raise an issue and we will see about adding more features.

Bed data will have some challenges. For example, if the contig names in the bed file don't match the assumed contig name from the fasta file, how will rust be able to know which is what? To make things work, we are going to assume, for now, that the bed file contig names match the names as derived by `rneat`. To determine this:

First, in a linux terminal, use your input fasta to find the names of the contigs like this:

```
$ grep "^>" input.fasta
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
>NC_001134.8 Saccharomyces cerevisiae S288C chromosome II, complete sequence
>NC_001135.5 Saccharomyces cerevisiae S288C chromosome III, complete sequence
>etc
```

This will give you the fill fasta name. `rneat` determines the short name of the input fasta by skipping the initial '>' character, then taking the string up to the first delimiter (space or '|'). In the above example, the short names are "NC_001133.9", "NC_001134.8", etc. 

Next, we check the bed. This example awk command will look for unique values in the first column of your bed file and print out what it finds. 

```
$ awk '!seen[$1]++ {print $1}' file.bed
1
2
3
etc
```

In this case, the bed file uses a simple numbering scheme to number the chromosomes. We can see a mapping with the roman numeral chromosomes in the name, but it's tricky to program. So, the following command will allow you, the user, to map these chromosomes and thus instruct NEAT. Yours will vary based on the inputs.

```
awk -F'\t' '{$1=($1=="1"?"NC_001133.9":$1); $1=($1=="2"?"NC_001134.8":$1); print}' OFS='\t' file.bed > file_renamed.bed
```

You will have to tailor this command to your dataset, but this will allow you to map the names from the bed to the fasta in a way `rneat` will understand.

Filtering Your Data
====================

To run filter-reads, you must first copy the filter_reads_template.yml file from rusty-neat/template_config/ to a directory of your choosing. Then you can edit the  file in your favorite editor. The configuration only has three fields:

```
bed_file: # required
```
A filename (full path is best) for a bed file with regions to target for filtering. It should be tab-separated in standard bed format.

In configuration:
```
files_to_filter: [
  # required, in list form
]
```
Filenames (full path is best) for files to filter by the bed file. For example, if you ran paired-ended fastqs with the vcf, your list would look like:

```
files_to_filter: [
   /my/home/output/bacteria_r1.fastq.gz,
   /my/home/output/bacteria_r2.fastq.gz,
   /my/home/output/bacteria.vcf.gz,
]
```
Note the full extension is important, but filter-reads should be able to handle both gzipped and unzipped files, if you decide to unzip them.

```
filter_key: .
```
This is the key appended to the filename, before extensions. For example, if your filename is "miscanthus_r1.fast.gz", and you set the key as "_only_genes", the output file would be "miscanthus_r1_only_genes.fastq.gz. The default is "_filtered".

Once your files are entered and your config is saved, you can run rneat:

```
$ rneat filter-reads -c my_config.yml
```
Your output will give you the filenames and success status.

IMPORTANT: This feature only works on fastq files generated by `rneat`. It is untested on generic VCF files, but should work. The reason it only works with `rneat`-generated fastq files is that `rneat` uses a naming scheme that the parser can use to identify where the read is located (enabling filtering), whereas an average fastq will not have that info.

Generating a Mutation Model
====================
`rneat` now includes the ability to read in real data and learn the parameters to reproduce it in a simulation run. So far, this is limited to the mutation model, usinge the `rneat gen-mut-model` subcommand. 

```
$ rneat gen-mut-model -c gen_mut_model_config.yml
```
The inputs for this are an arbitrary VCF file with a list of mutations on a single individual and the original fasta file that the vcf was built from. `rneat` will compute statistics for indels and snps and build the trinucleotide model of SNP generation, as in the original NEAT. Please submit any feature requests you have for expanding this feature.

Caveats: Only one sample can be read at this point. Currently, the high-mutation regions and common variants features from NEAT are not yet implemented. Please submit or comment on a feature request if you desire this feature.

Try it out and let us know if you run into any issues!

We will continue to improve this process in the future. If you have suggestions for parsing names to facilitate this, let us know in the Issues section!

That's it so far! Test out the features and let us know what you think!