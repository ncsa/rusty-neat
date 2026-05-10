# The rusty-neat project
Welcome to `rneat`, a Rust port of NEAT (https://github.com/ncsa/neat), a genetic simulation program that creates next-gen-looking fastq files along with an accompanying golden VCF and BAM files showing all of the variants inserted by `rneat` and the ideal alignments. It also generates "noise" in the form of sequencing errors as it is writing out files. These features can help you hone in your alignment and variant calling software to your data. Training models on your data will allow `rneat` to faithfully reproduce the statistical properties of your dataset.

The most recent version of `rneat` includes adding in the utilities to allow users to generate their own models for running `rneat`. It also now outputs a golden BAM file with ideal alignments, accepts a BED file for filtering the genome down to target regions for read creation, a more accurate coverage tool, and the ability to read in custom variants from a VCF file and insert them into the code. We swapped out the temp file writing method for a file streaming method that seems to work even faster and does not take up as much space on disk, avoiding expensive disk I/O of previous versions. We've done numerous benchmarks and tests on the data, but we look forward to user feedback on how well the output reproduces their data on a statistical level,

# How to use `rneat`

## Prerequisites
You will need to install the rust toolchain to compile `rneat`, including `cargo`. Check the cargo documentation for instructions (https://doc.rust-lang.org/cargo/getting-started/installation.html). Alternatively, you can try one of the binaries on the release page. Select the one that matches your system and let us know if you run into errors. During compilation, you may run into errors, such as cmake not found. Some of the packages `rneat` uses have these dependencies. For Debian/Ubuntu this should be a simple `sudo apt install cmake` and for RHEL/Rocky type distros this should be `sudo dnf install cmake`. There may be some other requirements. Drop a comment if you need specific help.

Download the executable in the release (current version 1.4.1).

```bash
$ rneat --help
Usage: rneat [OPTIONS] [SUB-COMMAND]

SUB-COMMANDS:
  gen-reads              Generates reads for an input dataset
  filter-reads           Filters the output of gen-reads
  gen-mut-model          Generates a mutation model from real VCF data
  gen-seq-error-model    Generates a sequencing error model from real FASTQ data
  gen-frag-length-model  Generates a fragment length model from a BAM or SAM file
  gen-gc-bias-model      Generates a GC bias model from a reference FASTA and coverage file
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

Input Variants VCF
==================
You can supply a VCF of variants to force into the simulation:

```yaml
input_vcf: /path/to/variants.vcf.gz
```

`rneat` will place every variant from the VCF into the corresponding position in the simulated reads and the output VCF. Random variants are still generated at `mutation_rate` for all positions not covered by the input VCF; set `mutation_rate: 0` to disable random variants entirely and simulate only the provided set.

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
`rneat` globally shuffles all generated reads before writing the final FASTQ by default, so no chromosome ordering is visible in the output — matching real sequencer output. Set `shuffle_fastq: false` in `gen-reads` config to preserve the generated ordering and avoid the in-memory global shuffle.

**Large genome note:** The global shuffle loads every read record into RAM at once. This is fine for small to moderate genomes (viral, bacterial, small eukaryotes), but becomes impractical for mammalian-scale genomes at typical coverage depths (e.g. human 30× ≈ 900 M reads ≈ hundreds of GB). For those cases, run the post-processing shuffle instead:

```bash
# single-ended
seqkit shuffle sample.fastq.gz -o sample_shuffled.fastq.gz

# paired-ended (keeps mates in sync)
seqkit shuffle -2 sample_R1.fastq.gz sample_R2.fastq.gz \
    -o sample_R1_shuffled.fastq.gz -o sample_R2_shuffled.fastq.gz
```

`seqkit` uses reservoir sampling and streams from disk, so its memory use is bounded regardless of file size. `rneat` will emit a warning at runtime when the reference genome exceeds 500 Mbp as a reminder that post-processing may be preferable.

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

Generating a GC Bias Model
====================
`rneat` can learn a GC bias model from a reference FASTA and a per-base coverage file using the `rneat gen-gc-bias-model` subcommand. It tiles the reference in fixed-size windows, computes the GC% of each window, looks up the mean coverage over that window, and accumulates the data into a 101-bin weight table (one bin per integer GC percentage, 0–100%). The resulting model can be passed to `gen-reads` to make fragment start positions favour regions whose GC content matches the coverage bias observed in your data.

```bash
$ rneat gen-gc-bias-model -c gen_gc_bias_model_config.yml
```

The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `gc_bias_model` config key.

Copy `template_config/gen_gc_bias_model.yml` to a directory of your choosing and fill in the fields:

```yaml
# required: reference FASTA used to compute GC content
# Both plain and gzipped (.fa.gz) references are accepted.
reference: /path/to/reference.fa

# required: per-base coverage file produced by samtools or bedtools
coverage_file: /path/to/coverage.txt

# required: format of the coverage file
# samtools-depth        → samtools depth output  (CHROM  1-based-pos  depth)
# bedtools-genomecov-d  → bedtools genomecov -d  (CHROM  1-based-pos  depth)
# bedtools-genomecov-dz → bedtools genomecov -dz (CHROM  0-based-pos  depth, nonzero only)
coverage_format: samtools-depth

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

**Generating the coverage file:**

The recommended command is:

```bash
samtools depth -a -q 20 -F 1796 aligned.bam > coverage.txt
```

- `-a` outputs all positions including zero-depth sites (required for `samtools-depth` and `bedtools-genomecov-d` formats)
- `-q 20` excludes reads with mapping quality below 20
- `-F 1796` excludes unmapped, secondary, duplicate, and supplementary reads

Getting these filters wrong silently produces a biased model, so use this command as-is unless you have a specific reason to deviate.

If you prefer bedtools or need smaller files:

```bash
# bedtools genomecov -d (1-based, all positions)
bedtools genomecov -d -ibam aligned.bam > coverage.txt

# bedtools genomecov -dz (0-based, nonzero positions only; smaller files)
bedtools genomecov -dz -ibam aligned.bam > coverage.txt
```

Note that bedtools does not expose MAPQ or flag filtering as simply as samtools; if duplicate removal or MAPQ filtering matters for your dataset, pre-filter the BAM first:

```bash
samtools view -b -q 20 -F 1796 aligned.bam > filtered.bam
samtools index filtered.bam
bedtools genomecov -dz -ibam filtered.bam > coverage.txt
```

Set `coverage_format` in your config to match whichever command you used.

**Coverage file must be plain text.** Gzipped coverage files are not supported. If yours is compressed, decompress it first:

```bash
gunzip coverage.txt.gz
```

**Large genomes:** The tool is designed to handle human-scale genomes without excessive memory use. It indexes the coverage file once, then loads each contig's data separately rather than holding the entire genome in RAM. Peak memory use is proportional to the longest contig, not total genome size.

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

**Zenodo release**

[![DOI](https://zenodo.org/badge/765847780.svg)](https://doi.org/10.5281/zenodo.20100558)