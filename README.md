# The rusty-neat project
Welcome to `rneat`, a Rust port of NEAT (https://github.com/ncsa/neat), a genetic simulation program that creates next-gen-looking fastq files along with an accompanying golden VCF and BAM files showing all of the variants inserted by `rneat` and the ideal alignments. It also generates "noise" in the form of sequencing errors as it is writing out files. These features can help you hone in your alignment and variant calling software to your data. Training models on your data will allow `rneat` to faithfully reproduce the statistical properties of your dataset.

The most recent version of `rneat` includes adding in the utilities to allow users to generate their own models for running `rneat`. It also now outputs a golden BAM file with ideal alignments, accepts a BED file for filtering the genome down to target regions for read creation, a more accurate coverage tool, and the ability to read in custom variants from a VCF file and insert them into the code. We swapped out the temp file writing method for a file streaming method that seems to work even faster and does not take up as much space on disk, avoiding expensive disk I/O of previous versions. We've done numerous benchmarks and tests on the data, but we look forward to user feedback on how well the output reproduces their data on a statistical level,

# How to use `rneat`

## Prerequisites
You will need to install the rust toolchain to compile `rneat`, including `cargo`. Check the cargo documentation for instructions (https://doc.rust-lang.org/cargo/getting-started/installation.html). Alternatively, you can try one of the binaries on the release page. Select the one that matches your system and let us know if you run into errors. During compilation, you may run into errors, such as cmake not found. Some of the packages `rneat` uses have these dependencies. For Debian/Ubuntu this should be a simple `sudo apt install cmake` and for RHEL/Rocky type distros this should be `sudo dnf install cmake`. There may be some other requirements. Drop a comment if you need specific help.

Download the executable in the release (current version 1.5.0).

```bash
$ rneat --help
Usage: rneat [OPTIONS] [SUB-COMMAND]

SUB-COMMANDS:
  gen-reads              Generates reads for an input dataset
  filter-reads           Filters the output of gen-reads
  gen-mut-model          Generates a mutation model from real VCF data
  gen-seq-error-model    Generates a sequencing error model from real FASTQ data
  gen-frag-length-model  Generates a fragment length model from a BAM or SAM file
  gen-gc-bias-model      Generates a GC bias model from a reference FASTA and aligned BAM
  gen-bam-models         Builds multiple models (frag-length, GC bias) from one BAM in a single pass
  help                   Print this message or the help of the given subcommand(s)

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

To compile and run `rneat` yourself, besides the rust toolchain, you will need `git` installed for your operating system. You will then need to git clone and cd into the repo directory. From your home directory in Linux the process might look something like:

```bash
~/$ git clone git@github.com:ncsa/rusty-neat.git
~/$ cd rusty-neat
~/rusty-neat/$
```

For Windows and Mac users, try the binary packaged with the latest version of `rneat`.

Once in the repo, you can build the program either in debug (default) or release mode. The main difference is how much info it gives you if there is an error. Release mode also has some optimizations to run it faster.

```bash
~/rusty-neat/$ cargo build --release
```

If you prefer to run the package directly without using the binary, you can also use
```bash
~/rusty-neat/$ cargo run -- gen-reads -c my_config.yml
```

Rust will download any required packages. Compiling Rust code is the slowest part of the process. The final binary will be built and the program run immediately after in the second case. To run the program manually, from the repo main dir, run

```bash
~/rusty-neat/$ ./target/release/rneat -h
```

`rneat` uses a configuration file to read values it needs for the run. A command line execution might look like this:

```bash
~/rusty-neat/$ ./target/release/rneat -c /path/to/filled/in/config.yml
```
If you record the output in the logs of Seed string to regenerate these exact results: XXXXXXX, you should be able to use that string as input with rng_seed and reproduce your results.

Fastq Output
============
The fastq output will have a key name that identifies the block where the read was drawn from, for quick comparisons in alignments. The output BAM file will contain the original sequence and cigar string. The name will have the format `RNEAT_generated_<contig short name>_<fragment_start>_<fragment_end>/1` (or `/2` for the second read in a pair), where start and end are zero-padded to 10 digits.

```bash
@RNEAT_generated_Chromosome_0000000000_0000000353/1
CTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTTAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGC
+
>AC7<GDEGGGGEFGA<GFCGG;GGGGGF>GEGGEGGGFGFEFGCEGGGGGGCG:AEFGFFGG>FG;GDGA9$GGAGF=GFG=EFFCGGGFGGGGGC$BFGEFFAGGG9F7E@>?GFGGGG>EBFGFDG)DGC6DEDFA2EG:EGG%FFB
```
The above example comes from the chromosome in the reference (which was simply named "Chromosome") between 0 and 353, which is the size
of a fragment in the simulated DNA. If there is a pair with this read, it will have the same coordinates, though it started at index 352 instead.

BAM Output
==========
`rneat` can write a golden BAM file alongside the FASTQ output. The BAM contains the same reads as the FASTQ — same sequences, same quality scores, same variants and sequencing errors applied — with alignment information included. To enable it, set `produce_bam: true` in your `gen-reads` config. The output path is derived automatically from `output_filename` (e.g. `output_filename: my_run` → `my_run.bam`). Please note that the BAM has a higher overhead than the fastq, and may take longer to produce.

```yaml
produce_bam: true
```

The CIGAR strings in the BAM reflect the full ground-truth alignment:
- Genomic variants (SNPs, insertions, deletions) are encoded as `M`, `I`, and `D` ops.
- Sequencing error indels are also encoded: deletion errors add `D` ops for skipped reference bases; insertion errors add `I` ops for inserted bases.

The BAM is written in coordinate-sorted order. No post-processing sort is required before indexing:

```bash
samtools index my_run.bam
```

Note: heterozygous variants are applied probabilistically, identical to the FASTQ. A read drawn to the reference allele will not show the variant in either the FASTQ or the BAM.

Targeted Generation with a BED File
=====================================
For exome, targeted-panel, or any run where you only want reads over specific regions, you can supply a BED file at generation time:

```yaml
target_bed: /path/to/targets.bed
```

When `target_bed` is set, `gen-reads` skips contigs absent from the BED entirely — no blocks are read, no variants are placed, and no reads are generated for those contigs. Within covered contigs, reads and variants are generated only over the intersecting non-N regions. This is the recommended approach for targeted runs on large genomes; it is far more efficient than generating genome-wide reads and post-filtering with `filter-reads`.

The BED contig names must match the short names derived from the reference FASTA (the text after `>` up to the first whitespace character). Both `.bed` and `.bed.gz` inputs are accepted.

Custom Mutation Rate Regions with a BED File
============================================
In addition to targeting certain regions, you can use a BED file with a custom field entered in any column after the third. 

```yaml
mutation_regions: /path/to/mutation_regions.bed
```

The text pattern is "mut_rate=0.0001" followed by a delimiter (end of line, space, semicolon, comma, bar) where the number after the equal sign is any float, and will be the mutation rate for the region defined in columns 1-3. 
```text
chr1    39930   39957   Bed_info    mut_rate=0.002
chr2    0   111199  Other_bed_info  mut_rate=0.02
```
Any region outside of the regions defined in the mutation regions BED will be assigned the default mutation rate, set by the mutation model, which can be overridden a custom default defined in the config file. For example, if you only want to have mutations in exomes, but you want reads from the full genome, you could set the `mutation_rate` to 0.0 and use a bed file defining all exomes, with a column appende `mut_rate=0.0011` or whatever mutation rate you desire. This BED works in conjunction with the `target_bed`, and any region outside the target bed will be excluded from reads anyway, and thus will have no variants. 

Input Variants VCF
==================
You can supply a VCF of variants to force into the simulation:

```yaml
input_vcf: /path/to/variants.vcf.gz
```

`rneat` will place every variant from the VCF into the corresponding position in the simulated reads and the output VCF. Random variants are still generated at `mutation_rate` for all positions not covered by the input VCF; set `mutation_rate: 0.0` to disable random variants entirely and output only the provided set of variants.

**Requirements**

- The VCF must be single-sample (one sample column).
- Every record must include `GT` in the FORMAT field. `rneat` uses the genotype to determine whether to apply the variant to all reads covering the position (homozygous, e.g. `1/1`) or only a probabilistic subset (heterozygous, e.g. `0/1`). Records without `GT` are rejected.
- Contig names must match the short names derived from the reference FASTA (text after `>` up to the first whitespace character). Variants on unrecognised contigs are skipped with a warning.
- Both `.vcf` and `.vcf.gz` files are accepted.

**Supported variant types**

| Type | Condition | Handled |
|------|-----------|---------|
| SNP | REF and ALT both single base | Yes |
| Insertion | single-base REF, multi-base ALT | Yes |
| Deletion | multi-base REF, single-base ALT | Yes |
| Complex | multi-base REF **and** multi-base ALT | **No** — skipped with warning |

**Current caveats**

- *Multi-allelic records*: only the first ALT allele is used; additional alleles are silently ignored. Split multi-allelic records with `bcftools norm -m -` before passing to `rneat`.
- *REF allele verification*: `rneat` does not check that the REF field matches the reference sequence at that position. Mismatches will produce biologically incorrect output without any warning.
- *Structural variants / breakends*: records whose REF and ALT are both multi-base strings (complex variants) are skipped with a logged warning and do not appear in the output.

FASTQ Shuffling
===============
`rneat` writes reads in contig order. To shuffle the output, use `seqkit shuffle` as a post-processing step:

```bash
# single-ended
seqkit shuffle sample.fastq.gz -o sample_shuffled.fastq.gz

# paired-ended (keeps mates in sync)
seqkit shuffle -2 sample_R1.fastq.gz sample_R2.fastq.gz \
    -o sample_R1_shuffled.fastq.gz -o sample_R2_shuffled.fastq.gz
```

`seqkit` uses reservoir sampling and streams from disk, keeping memory use bounded regardless of file size (`seqkit` is an open source toolkit for FASTQ files: https://github.com/shenwei356/seqkit).

Parallel Processing
===================
`rneat gen-reads` processes contigs in parallel by default using rayon's work-stealing thread pool. Each contig is an independent unit of work — variant generation, fragment sampling, and FASTQ writing all happen concurrently across contigs — so multi-core machines see roughly linear speedup up to the number of contigs in the reference.

**Thread count:**

By default `rneat` uses all available logical cores. You can cap the thread count with the `num_threads` config key:

```yaml
# use 4 threads instead of all available cores
num_threads: 4
```
or disable parallelism entirely (useful for debugging or reproducibility testing)
```yaml
# disable parallelism entirely (useful for debugging or reproducibility testing)
num_threads: 1
```

Omit `num_threads` (or set it to `.`) to restore the default all-cores behaviour.

**BAM output is fully parallel:**

`rneat` uses a per-contig temp-file strategy: each contig worker writes its alignment records to a private temporary BAM body file, then a single concatenation pass assembles them in reference order into the final coordinate-sorted BAM.

**Reproducibility:**

Each contig's random number generator is derived deterministically from the parent seed and the contig's position in the reference, so output is identical across runs with the same seed even when the number of threads changes.

Reading Bed Data
================
`rneat` can now read bed data and filter the reads based on regions specified. This comes with some caveats. First is that rust can't handle any header lines or non-header rows in the bed file, because of potential range of possibilities. Second, the names must match what is in the fasta. Third, `rneat` can only filter, and does not use the rest of the bed information. It treats each record as a region of interest only.

Bed data will have some challenges. For example, if the contig names in the bed file don't match the assumed contig name from the fasta file, how will rust be able to know which is what? To make things work, we are going to assume, for now, that the bed file contig names match the names as derived by `rneat`. To determine this:

First, in a linux terminal, use your input fasta to find the names of the contigs like this:

```bash
$ grep "^>" input.fasta
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
>NC_001134.8 Saccharomyces cerevisiae S288C chromosome II, complete sequence
>NC_001135.5 Saccharomyces cerevisiae S288C chromosome III, complete sequence
>etc
```
This will give you the full fasta name. `rneat` determines the short name of the input fasta by skipping the initial '>' character, then taking the string up to the first whitespace delimiter. In the above example, the short names are "NC_001133.9", "NC_001134.8", etc. 

Next, we check the bed. This example awk command will look for unique values in the first column of your bed file and print out what it finds. 

```bash
$ awk '!seen[$1]++ {print $1}' file.bed
1
2
3
etc
```
In this case, the bed file uses a simple numbering scheme to number the chromosomes. We can see a mapping with the roman numeral chromosomes in the name, but it's tricky to handle consistently. So, the following command will allow you, the user, to map these chromosomes and thus instruct `rneat`. Yours will vary based on the inputs.

```bash
awk -F'\t' '{$1=($1=="1"?"NC_001133.9":$1); $1=($1=="2"?"NC_001134.8":$1); print}' OFS='\t' file.bed > file_renamed.bed
```
You will have to tailor this command to your dataset, but this is one way to map the names from the bed to the fasta in a way `rneat` will understand.

Filtering Your Data
====================

To run filter-reads, you must first copy the filter_reads_template.yml file from rusty-neat/template_config/ to a directory of your choosing. Then you can edit the file in your favorite editor. The configuration has four fields:

```bash
bed_file: /path/to/my.bed # required
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
This is the key appended to the filename, before extensions. For example, if your filename is "miscanthus_r1.fastq.gz", and you set the key as "_only_genes", the output file would be "miscanthus_r1_only_genes.fastq.gz". The default is "_filter".

```bash
overwrite_output: false # optional, default false
```
Set to `true` to allow filter-reads to overwrite existing output files. When `false`, the run will error if the output file already exists.

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

The inputs are a single-sample VCF and the reference FASTA the VCF was called against. `rneat` computes statistics for indels and SNPs and builds the trinucleotide model of SNP generation, as in the original `rneat`. The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `mutation_model` config key.

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

Caveats: Only one sample can be read at this point. Currently, high-mutation regions and common variants features from Python NEAT are not yet implemented.

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

# optional: list of Q-score bins to quantize the learned model to (e.g. NovaSeq 6000 uses
# [2, 12, 23, 37]). When set, each observed Q-score is snapped to its nearest bin before
# the model is built, and gen-reads will emit only these values. Omit for a continuous model.
# Bins must be < 94 and cannot include 31 (encodes to '@' under Phred+33).
# binned_quality_bins: [2, 12, 23, 37]

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
3. Default matrix from Python `rneat` — used when neither is provided

**BAM MD tag requirement:**
The BAM path requires MD tags to identify reference bases at mismatch positions. Most aligners (BWA-MEM, STAR with `--outSAMattributes MD`) add them automatically. If yours does not, generate them with:
```bash
samtools calmd -b aligned.bam reference.fa > aligned_with_md.bam
samtools index aligned_with_md.bam
```

**TSV format:**
The transition matrix TSV has 4 data rows (one per reference base: A, C, G, T) and 4 whitespace-separated columns (one per read base: A, C, G, T). The first row is skipped if it is non-numeric (treated as a header). Diagonal values are zeroed automatically. Rows are re-normalized to sum to 1.

```
A    C    G    T
0.0  0.5  0.3  0.2
0.5  0.0  0.3  0.2
0.4  0.3  0.0  0.3
0.3  0.3  0.4  0.0
```

**Binned quality scores:**
Modern Illumina platforms (NovaSeq 6000, NextSeq 1000/2000, the simplified MiSeq reporting modes) no longer emit the full continuous Q0–Q40 range. Instead, they quantize each per-base quality into a small set of discrete bins — for example, NovaSeq 6000 emits only `{2, 12, 23, 37}`. Set `binned_quality_bins` in the config to mimic this behaviour:

```yaml
binned_quality_bins: [2, 12, 23, 37]
```

When this field is set, `gen-seq-error-model` snaps each observed Q-score from the input FASTQ to its nearest bin (ties round down) before learning seed, transition, and global-frequency counts. The resulting model file is flagged `binned_scores: true` and lists the bin values as its `quality_score_options`. `gen-reads` then samples only those bin values when emitting reads — no extra flag needed on the read-generation side.

Validation rules:
- Bins must be non-empty integers in `[0, 94)` (entries are sorted and deduped automatically).
- Value `31` is rejected — under Phred+33 it encodes to `@`, which would corrupt FASTQ output. If you need a bin near Q31, pick 30 or 32.
- Bins with no observed counts in the training FASTQ are still kept in `quality_score_options`; their transition rows fall back to a uniform distribution, and a warning is logged listing the empty bins.

Some common platform bin sets (consult your sequencer's documentation for the authoritative list):

| Platform                       | Suggested bins     |
|--------------------------------|--------------------|
| NovaSeq 6000 (4-bin)           | `[2, 12, 23, 37]`  |
| NextSeq 1000/2000 (3-bin)      | `[2, 15, 35]`      |
| NextSeq 1000/2000 (4-bin)      | `[2, 12, 23, 37]`  |

Generating a GC Bias Model
====================
`rneat` can learn a GC bias model from a reference FASTA and an aligned BAM file using the `rneat gen-gc-bias-model` subcommand. It walks the BAM once to accumulate per-base reference coverage, tiles the reference in fixed-size windows, computes the GC% of each window, looks up the mean coverage over that window, and accumulates the data into a 101-bin weight table (one bin per integer GC percentage, 0–100%). The resulting model can be passed to `gen-reads` to make fragment start positions favour regions whose GC content matches the coverage bias observed in your data.

```bash
$ rneat gen-gc-bias-model -c gen_gc_bias_model_config.yml
```

The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `gc_bias_model` config key.

Copy `template_config/gen_gc_bias_model.yml` to a directory of your choosing and fill in the fields:

```yaml
# required: reference FASTA used to compute GC content
# Both plain and gzipped (.fa.gz) references are accepted.
reference: /path/to/reference.fa

# required: aligned BAM. Coverage is accumulated in process per reference base.
bam_file: /path/to/aligned.bam

# optional: minimum mapping quality applied while accumulating coverage.
# Default 0 matches `samtools depth` defaults.
min_mapq: 0

# optional: BED file restricting model inference to target regions
# Use "." (dot) to use the entire reference (default)
bed_file: .

# required: output path for the generated model; should end in .json.gz
output_file: /path/to/gc_bias_model.json.gz

# optional: set to true to overwrite an existing output file (default: false)
overwrite_output: false

# Window size in base pairs for GC% and coverage averaging.
# Should match the typical read length (or fragment length for long reads).
# Default: 100
window_size: 100

# Step between successive windows. Must not exceed window_size.
# Use the same value as window_size for non-overlapping windows (recommended).
# Values smaller than window_size produce overlapping windows.
# Default: window_size
window_stride: 100

# Bins with fewer windows than this receive a neutral weight of 1.0 rather than
# a learned weight. Increase this to require more evidence before trusting a bin.
# Default: 10
min_windows_per_bin: 10
```

**BAM filtering:** unmapped, secondary, and supplementary records are skipped automatically. Reads with mapping quality `<= min_mapq` are dropped before contributing to coverage. The default `min_mapq: 0` reproduces `samtools depth` defaults; set it to 20 if your downstream `gen-reads` runs assume mapq-filtered coverage.

**Large genomes:** The tool walks the BAM once and allocates a per-contig depth array (`u32` per reference base) only for contigs that receive at least one record. Peak memory for hg38 at 30× is on the order of all-contigs-summed depth arrays — see the HPC section below for sizing.

**Long reads:** Set `window_size` to approximately your typical read length (e.g. 5000–50000 for ONT or PacBio). Larger windows mean fewer windows per contig and faster model building. Non-overlapping windows (`window_stride: window_size`) are recommended so each observation is independent.

**If the model has no effect in gen-reads:** If every GC% bin in your reference has fewer observations than `min_windows_per_bin`, all weights will be neutral (1.0) and no bias will be applied. This is logged as a warning. To diagnose, lower `min_windows_per_bin` or increase the region covered (use a larger reference or remove the BED restriction).

**Using the model in gen-reads:**

```yaml
# in gen_reads_config.yml
gc_bias_model: /path/to/gc_bias_model.json.gz

# optional: inflate fragment count to compensate for low-weight regions
# true (default): total coverage stays close to the requested depth
# false: coverage in low-GC regions will be lower than requested
gc_bias_normalize_coverage: true
```

Generating a Fragment Length Model
====================
`rneat` can learn a fragment length model from real paired-end alignment data using the `rneat gen-frag-length-model` subcommand. It reads a BAM or SAM file, collects the template lengths (TLEN) of confidently-mapped concordant read pairs, filters out rare and extreme-outlier lengths, then fits a normal distribution to the result and writes a `FragmentLengthModel` that `gen-reads` can use.

```bash
$ rneat gen-frag-length-model -c gen_frag_length_model_config.yml
```

The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `fragment_model` config key.

Copy `template_config/gen_frag_length_model_template.yml` to a directory of your choosing and fill in the fields:

```yaml
# required: path to input BAM or SAM file
# The file must contain paired-end aligned reads. BAM is strongly recommended.
input_file: /path/to/aligned.bam

# required: output path for the generated model; should end in .json.gz
output_file: /path/to/fragment_length_model.json.gz

# optional: minimum number of reads for a fragment length to be included in the model
# Default is 2 to handle smaller datasets. Set to 0 to disable filtering entirely.
# For larger datasets (e.g. whole-genome), try 100 and adjust from there.
min_reads: 2

# optional: set to true to overwrite an existing output file (default: false)
overwrite_output: false
```

**Read filtering rules:**  
Only reads that satisfy all of the following are used:
- Paired-end and first in pair
- Mapped (not unmapped, secondary, or supplementary)
- Mate mapped to the same reference sequence
- Mapping quality > 10

**Outlier filtering:**  
Fragment lengths that exceed `median + 10 × MAD` (median absolute deviation) are removed before fitting. This mirrors the filtering in Python NEAT's fragment length modeler. If no lengths survive the filter, lower `min_reads` or set it to `0`.

Building multiple BAM-derived models in one pass
====================
If you need both a fragment length model and a GC bias model from the same BAM, `rneat gen-bam-models` walks the BAM once and feeds both observers from the single pass — half the I/O of running the two per-tool commands back-to-back.

```bash
$ rneat gen-bam-models -c gen_bam_models_config.yml
```

Copy `template_config/gen_bam_models.yml` and fill in the sections for whichever models you want:

```yaml
bam_file: /path/to/aligned.bam
min_mapq: 0   # walker-level filter; FragLengthObserver enforces MAPQ > 10 internally

frag_length:
  output_file: /path/to/frag_length_model.json.gz
  overwrite_output: false
  min_reads: 2

gc_bias:
  reference: /path/to/reference.fa
  output_file: /path/to/gc_bias_model.json.gz
  overwrite_output: false
  bed_file: .
  window_size: 100
  window_stride: 100
  min_windows_per_bin: 10
```

Both sections are optional but at least one must be present. The output models are interchangeable with what `gen-frag-length-model` and `gen-gc-bias-model` would produce on their own, so the resulting files plug into `gen-reads` exactly the same way.

Running on HPC
==============

`rneat` runs as a single process and fits naturally onto a single HPC compute node. No MPI, distributed computing, or special environment setup is required. The notes below cover resource budgets for whole-genome human-scale runs.

### gen-reads

`gen-reads` is already multi-threaded via rayon (see [Parallel Processing](#parallel-processing) above). Set `num_threads` in your config to match your CPU allocation:

```yaml
num_threads: 16
```

Memory scales with `num_threads × largest_contig_size`, not total genome size — each thread processes one contig at a time. For hg38, chr1 is ~249 Mbp, so 16 simultaneous threads need roughly 16 GB for sequence data plus overhead. Budget additional scratch space for per-contig temp files before the final assembly: roughly 3× the expected output size.

| Genome | Threads | Recommended RAM | Scratch space | Typical wall time |
|--------|---------|-----------------|---------------|-------------------|
| Bacterial (~5 Mbp) | 4 | 4 GB | 5 GB | < 1 min |
| Human hg38, 10× SE | 16 | 32 GB | 100 GB | ~20 min |
| Human hg38, 30× PE | 16 | 48 GB | 300 GB | ~60 min |

Example SLURM header:

```bash
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --tmp=300G
```

### gen-gc-bias-model

Single-threaded; walks the BAM once and accumulates per-base reference coverage as a `Vec<u32>` per observed contig. Peak memory is roughly `4 bytes × total reference bases that received at least one record`. For a full human genome (~3.2 Gbp) this is ~13 GB. For a single chromosome or a targeted/exome BAM it is much smaller.

If memory is tight on whole-genome runs, split the BAM by chromosome and run `gen-gc-bias-model` once per chromosome on a representative subset that spans your GC range of interest, then average or combine the resulting models.

Recommended SLURM header (whole-genome BAM):

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
```

### gen-mut-model

Single-threaded; processes the VCF one chromosome at a time. A full human genome WGS VCF (~10 GB) typically completes in 30–60 minutes with peak RSS under 4 GB.

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:30:00
```

### gen-seq-error-model

Single-threaded; streams through the FASTQ. For a full genome FASTQ (~600 M reads), expect 30–60 minutes. For most purposes, training on a representative subset produces a statistically equivalent model in a fraction of the time:

```yaml
max_reads: 5000000   # 5 M reads is sufficient; set to 0 for unlimited
```

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00   # with max_reads=5M; scale up for unlimited
```

### gen-frag-length-model

Single-threaded; streams TLEN fields from the BAM. A full human genome BAM at 30× typically completes in under 15 minutes with peak RSS under 2 GB.

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00
```

### Environment notes

The `rneat` binary has no runtime dependencies beyond a standard C library (glibc), which is present on all Linux HPC systems. No module loads or conda environments are required — copy the release binary to your scratch or project directory and run it directly. If you compile from source on the cluster, ensure `cmake` is available (`module load cmake` on most systems) before running `cargo build --release`.

**Current Benchmarks**
We ran some benchmarks on `rneat` on my home desktop, a very basic Linux desktop, using data pulled from public resources on the internet (mostly ncbi). The "yeast" is brewer's yeast. These are the results:

[16:26:39] Build complete.
rneat Benchmark Report
Run:        20260510_162639
Machine:    pop-os  |  8 logical CPUs
Coverage:   10x    |  Read length: 151 bp
Binary:     /home/joshfactorial/code/rusty-neat/target/release/rneat
──────────────────────────────────────────────────────────────────────
SECTION 1: Single-ended read generation
Genome          Size(MB)   Wall time  Peak RSS(MB)      CPU%
──────────────────────────────────────────────────────────────────────
[16:26:40] Single-ended: ecoli  (4.6 MB)
ecoli                4.6     0:14.83         416.9       99%
[16:26:54]   Done: wall=0:14.83  peak_rss=416.9 MB  cpu=99%
[16:26:54] Single-ended: pneumonia  (2.1 MB)
pneumonia            2.1     0:05.39         197.0       99%
[16:27:00]   Done: wall=0:05.39  peak_rss=197.0 MB  cpu=99%
[16:27:00] Single-ended: yeast  (12 MB)
yeast                 12     0:10.94         981.7      352%
[16:27:11]   Done: wall=0:10.94  peak_rss=981.7 MB  cpu=352%
[16:27:11] Single-ended: c_elegans  (97 MB)
c_elegans             97     7:35.23        8591.2      387%
[16:34:46]   Done: wall=7:35.23  peak_rss=8591.2 MB  cpu=387%
──────────────────────────────────────────────────────────────────────
Section 1 — Passed: 4 / 4  |  Failed: 0
──────────────────────────────────────────────────────────────────────
SECTION 2: Paired vs single comparison  (genome: yeast, 10x)
  Paired-end uses fragment_mean=400, fragment_st_dev=50
  Same reference, same coverage, same seed — only the read-generation mode differs.
Mode               Wall time  Peak RSS(MB)      CPU%
──────────────────────────────────────────────────────────────────────
[16:34:46] Paired comparison — single-ended pass
single-ended         0:10.48         970.4      365%
[16:34:57]   Single done: wall=0:10.48  peak_rss=970.4 MB  cpu=365%
[16:34:57] Paired comparison — paired-ended pass
paired-ended         0:13.30        1126.8      354%
[16:35:10]   Paired done: wall=0:13.30  peak_rss=1126.8 MB  cpu=354%
──────────────────────────────────────────────────────────────────────
  Paired/single wall-time ratio: 1.27x  (10.48s → 13.30s)
  Paired/single peak-RSS ratio:  1.16x  (970.4 MB → 1126.8 MB)
──────────────────────────────────────────────────────────────────────

**Zenodo release**

[![DOI](https://zenodo.org/badge/765847780.svg)](https://doi.org/10.5281/zenodo.20100558)