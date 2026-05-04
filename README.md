# The rusty-neat project
Weclome to `rneat`, a Rust port of NEAT (https://github.com/ncsa/neat), a genetic simulation program that creates next-gen-looking fastq files along with an accompanying golden vcf file showing all of the variants inserted by NEAT. It also generates "noise" in the form of sequencing errors as it is writing out files. These features can help you hone in your alignment software to your data. 

As of now, `rneat` is still under development. We are working toward the goal of the original NEAT of being able to read user data and construct a model based on that data. We do however, have functional reads being generated, with variants tracked in a vcf file. This is our 1.2.0 version.

# How to use `rneat`

## Prerequisites
One of the packages requires cmake to use. For Debian/Ubuntu this should be a simple `sudo apt install cmake` and for RHEL/Rocky type distros this should be `sudo dnf install cmake`

Download the executable in the release (current version 1.2.1).

```bash
$ rneat --help
Usage: rneat [OPTIONS] [SUB-COMMAND]

SUB-COMMANDS:
  gen-reads            Generates reads for an input dataset
  filter-reads         Filters the output of gen-reads
  gen-mut-model        Generates a mutation model from real VCF data
  gen-seq-error-model  Generates a sequencing error model from real FASTQ data
  help                 Print this message or the help of the given subcommand(s)

Options:
      --log-level [<log_level>]  Sets the log level for the display log. [default: trace] [possible values: trace, debug, info, warn, error, off]
      --log-dest <log_dest>      Sets the log destination (full path with full filename) for the written log
  -h, --help                     Print help
```

To check options for a subcommand:

```bash
$ rneat gen-reads --help
Generates reads for an input dataset

Usage: rneat gen-reads [OPTIONS]

Options:
  -c, --configuration-yaml <configuration_yaml>  Path to configuration file.
  -h, --help                                     Print help
```

To run filter reads, check the help menu.

```bash
$ rneat filter-reads --help
Filters the output of gen-reads

Usage: rneat filter-reads --configuration-yaml <configuration_yaml>

Options:
  -c, --configuration-yaml <configuration_yaml>  Path to configuration file.
  -h, --help                                     Print help
```

To run gen-mut-model:

```bash
$ rneat gen-mut-model --help
Generates a mutation model from input VCF data

Usage: rneat gen-mut-model [OPTIONS]

Options:
  -c, --configuration-yaml <configuration_yaml>  Path to configuration file.
  -h, --help                                     Print help
```

Use the help menu to see the available options and leave an issue if you find something bad happening.

To compile and run `rneat` yourself, you will need the Rust environment (https://www.rust-lang.org/tools/install), with cargo. You will also need git installed for your operating system. You will then need to git clone and cd into the repo directory. From your home directory in Linux the process might look something like:

```bash
~/$ git clone git@github.com:ncsa/rusty-neat.git
~/$ cd rusty-neat
~/rusty-neat/$
```

For Windows users, you're probably going to want to use WSL, but go ahead and test it if there's a Rust for windows. (Post help requests in the Issues tab if you need assistance or have suggestions).

Once in the repo, you can build the program either in debug (default) or release mode. The main difference is how much info it gives you if there is an error. Release mode also has some optimizations to run it faster.

```bash
~/rusty-neat/$ cargo build --release
```

If you prefer to run the package directly without using the binary, you can also use
```bash
~/rusty-neat/$ cargo run
```

Rust will download any required packages. Compiling Rust code is the slowest part of the process. The final binary will be built and the program run immediately after in the second case. To run the program manually, from the repo main dir, run

```bash
~/rusty-neat/$ ./target/release/rneat -h
```

`rneat` uses a configuration file to read values it needs for the run. Most of these features should be active, but there is currently no way to generate your own data, so stick with default models for now. Also, there is not yet a way to create bams. A command line execution might look like this:

```bash
~/rusty-neat/$ ./target/release/rneat -C /path/to/filled/in/config.yml
```
If you record the output in the logs of Seed string to regenerate these exact results: XXXXXXX, you should be able to use that string as input with rng_seed and reproduce your results (untested as of yet).

Fastq Output
============
the fastq output will have a key name that identifies the block where teh read was drawn from. This should allow you to asseses how well the aligner funcution. The name will have the format `<contig short name>_<fragment_start>_<fragment_end>`

```bash
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

```bash
$ grep "^>" input.fasta
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
>NC_001134.8 Saccharomyces cerevisiae S288C chromosome II, complete sequence
>NC_001135.5 Saccharomyces cerevisiae S288C chromosome III, complete sequence
>etc
```
This will give you the fill fasta name. `rneat` determines the short name of the input fasta by skipping the initial '>' character, then taking the string up to the first delimiter (space or '|'). In the above example, the short names are "NC_001133.9", "NC_001134.8", etc. 

Next, we check the bed. This example awk command will look for unique values in the first column of your bed file and print out what it finds. 

```bash
$ awk '!seen[$1]++ {print $1}' file.bed
1
2
3
etc
```
In this case, the bed file uses a simple numbering scheme to number the chromosomes. We can see a mapping with the roman numeral chromosomes in the name, but it's tricky to program. So, the following command will allow you, the user, to map these chromosomes and thus instruct NEAT. Yours will vary based on the inputs.

```bash
awk -F'\t' '{$1=($1=="1"?"NC_001133.9":$1); $1=($1=="2"?"NC_001134.8":$1); print}' OFS='\t' file.bed > file_renamed.bed
```
You will have to tailor this command to your dataset, but this will allow you to map the names from the bed to the fasta in a way `rneat` will understand.

Filtering Your Data
====================

To run filter-reads, you must first copy the filter_reads_template.yml file from rusty-neat/template_config/ to a directory of your choosing. Then you can edit the  file in your favorite editor. The configuration only has three fields:

```bash
bed-file: /path/to/my.bed # required
```
A filename (full path is best) for a bed file with regions to target for filtering. It should be tab-separated in standard bed format.

In configuration:
```bash
files_to_filter: [
  # required, in list form
]
```
Filenames (full path is best) for files to filter by the bed file. For example, if you ran paired-ended fastqs with the vcf, your list would look like:

```bash
files_to_filter: [
   /my/home/output/bacteria_r1.fastq.gz,
   /my/home/output/bacteria_r2.fastq.gz,
   /my/home/output/bacteria.vcf.gz,
]
```
Note the full extension is important, but filter-reads should be able to handle both gzipped and unzipped files, if you decide to unzip them.

```bash
filter_key: .
```
This is the key appended to the filename, before extensions. For example, if your filename is "miscanthus_r1.fast.gz", and you set the key as "_only_genes", the output file would be "miscanthus_r1_only_genes.fastq.gz. The default is "_filtered".

Once your files are entered and your config is saved, you can run rneat:

```bash
$ rneat filter-reads -c my_config.yml
```
Your output will give you the filenames and success status.

IMPORTANT: This feature only works on fastq files generated by `rneat`. It is untested on generic VCF files, but should work. The reason it only works with `rneat`-generated fastq files is that `rneat` uses a naming scheme that the parser can use to identify where the read is located (enabling filtering), whereas an average fastq will not have that info.

Generating a Mutation Model
====================
`rneat` now includes the ability to read in real data and learn the parameters to reproduce it in a simulation run. So far, this is limited to the mutation model, using the `rneat gen-mut-model` subcommand.

```bash
$ rneat gen-mut-model -c gen_mut_model_config.yml
```

The inputs are a single-sample VCF and the reference FASTA the VCF was called against. `rneat` computes statistics for indels and SNPs and builds the trinucleotide model of SNP generation, as in the original NEAT. The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `mutation_model` config key.

Copy `template_config/gen_mut_model_template.yml` to a directory of your choosing and fill in the fields:

```yaml
# required: path to the reference FASTA used to call the VCF
reference: /path/to/reference.fa

# required: path to a single-sample VCF file
# The VCF must include FORMAT and SAMPLE columns, and GT must be present in FORMAT.
# QUAL=. is accepted and treated as 0.
vcf_file: /path/to/variants.vcf

# required: output path for the generated model; should end in .json.gz
output_file: /path/to/output_model.json.gz

# optional: restrict model learning to regions in this BED file
# set to . (dot) to use the entire reference (default)
bed_file: .

# optional: set to true to overwrite an existing output file (default: false)
overwrite_output: false

# optional: custom 4x4 SNP transition matrix TSV (rows/columns: A C G T).
# A single header line is allowed. Diagonal values are zeroed automatically.
# Overrides the transition matrix inferred from VCF SNP data.
transition_matrix_file: /path/to/matrix.tsv
```

**VCF requirements:**
- Single sample only — multi-sample VCFs are not yet supported
- Each variant record must have `GT` in the `FORMAT` column; `rneat` hard-errors if GT is missing
- `QUAL=.` is accepted and treated as quality score 0

Caveats: Only one sample can be read at this point. Currently, high-mutation regions and common variants features from NEAT are not yet implemented. Please submit or comment on a feature request if you desire this feature.

Try it out and let us know if you run into any issues!

We will continue to improve this process in the future. If you have suggestions for parsing names to facilitate this, let us know in the Issues section!

Generating a Sequencing Error Model
====================
`rneat` can also learn a sequencing error model from real FASTQ data using the `rneat gen-seq-error-model` subcommand. This reads per-base quality scores from the FASTQ to build a Markov quality score model and computes an average base error rate. Optionally, you can supply a BAM file or a custom TSV to set the SNP transition matrix (which base errors are most likely to become which other bases).

```bash
$ rneat gen-seq-error-model -c gen_seq_error_model_config.yml
```

The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `seq_error_model` config key.

Copy `template_config/gen_seq_error_model_template.yml` to a directory of your choosing and fill in the fields:

```yaml
# required: path to input FASTQ file (.fastq or .fastq.gz)
fastq_file: /path/to/reads.fastq.gz

# required: output path for the generated model; should end in .json.gz
output_file: /path/to/seq_error_model.json.gz

# optional: set to true to overwrite an existing output file (default: false)
overwrite_output: false

# optional: maximum number of reads to use for model learning; 0 = unlimited (default: 0)
max_reads: 0

# optional: quality score ASCII offset; 33 for Illumina 1.8+/Sanger, 64 for older Illumina (default: 33)
qual_offset: 33

# optional: aligned BAM file to infer the SNP transition matrix from read-vs-reference mismatches.
# The BAM must have MD tags. If your aligner did not add them, run:
#   samtools calmd -b aligned.bam reference.fa > aligned_with_md.bam
# Ignored if transition_matrix_file is also set.
bam_file: /path/to/aligned.bam

# optional: custom 4x4 SNP transition matrix TSV (rows/columns: A C G T).
# A single header line is allowed. Diagonal values are ignored.
# Takes precedence over bam_file.
transition_matrix_file: /path/to/matrix.tsv
```

**SNP transition matrix priority:**
1. `transition_matrix_file` (explicit TSV) — highest priority
2. `bam_file` (inferred from MD-tagged BAM mismatches)
3. Default matrix from Python NEAT — used when neither is provided

**BAM MD tag requirement:**
The BAM path requires MD tags to identify reference bases at mismatch positions. Most aligners (BWA-MEM, STAR with `--outSAMattributes MD`) add them automatically. If yours does not, generate them with:
```bash
samtools calmd -b aligned.bam reference.fa > aligned_with_md.bam
samtools index aligned_with_md.bam
```

**TSV format:**
The transition matrix TSV has 4 data rows (one per reference base: A, C, G, T) and 4 whitespace-separated columns (one per read base: A, C, G, T). An optional header row (non-numeric first token) is skipped. Diagonal values are zeroed automatically. Rows are re-normalized to sum to 1.

```
A    C    G    T
0.0  0.5  0.3  0.2
0.5  0.0  0.3  0.2
0.4  0.3  0.0  0.3
0.3  0.3  0.4  0.0
```

That's it so far! Test out the features and let us know what you think!