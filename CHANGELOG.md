# RustyNeat v1.0.0

9/5/2023
=========
- Initialized Cargo package with hello world default script.
- Added function to select a random integer from a list (for selecting positions)

9/9/2012
=========
- Created a first pass program that mutates a hardcoded DNA string at one random non-N position.
- Added a logging function to print out time-stamped log items both to stdout and to file.
- Need to improve error handling all around, taking advantage of Rusts error handling system.
- Design decision: After toying with some options to create a standardized RNG to use throughout the program, it seems the best option is to use an on-the-fly rng per function for now. This eliminates reproducibility between runs, but should enable parallelism of processing.
- Redesign from the ground up. Trying to use KISS as much as possible, and only introduce complexity as it becomes necessary, or as I learn Rust tricks that may help this endeavor.

## rusty-neat v0.1

Goal: To release a version of NEAT in rust capable of reading an arbitrary fasta file, 
and outputting a version of that file with 1% of the bases mutated as SNPs.

9/11/2023
=========
- Found a sqlite crate for rust. This might be a good way to store the reference information for our genome. We can have one row per chromosome. Name, info, sequence. Read the data in once, store it for long term retrieval. Website: https://docs.rs/sqlite/latest/sqlite/

## rusty-neat v0.1.0

Moved the code from a brach under NEAT (https://github.com/ncsa/neat) to its own repo https://github.com/ncsa/rusty-neat. This project will be considered a separate, related project. We will use the same underlying design ideas as NEAT, but it will essentially be different and won't be guaranteed in any way to be able to reproduce what the original NEAT did (though we aim to get it as close as possible, the differences in Python and Rust would make byte-by-byte compatability very difficult). I have implemented an idea I prototyped in a private repo (https://github.com/joshfactorial/basic_neat), of reading in the file in what I hope is a straightforward way. I was looking for stuff like running out of memory and speed. It can read in the `E. coli` fasta file in about 30 seconds on my laptop, which is pretty fast, compared to Python. What's more is it can store the whole thing in memory, due to the efficient way Rust uses memory. Instead of copying the entire data set every time in memory, quickly eating it up, Rust forces you to decide when to pass a read-only version, copy the dataset, or use up the dataset once and lose it forever (all valid choices depending on context). So I'm feeling hopeful that this will work. I am publishing an Alpha called v0.1.0 for people to test, if they feel so inclined. See the README for more details.

9/9/2025
=========

## rusty-neat v1.0.0
- It's been a couple years since I had an update to rusty-neat. This new vrison is version 1.0.0! It is not yet fully functional, as it is missing the BAM creation, though it does work on some test data: H1N1.fa and ecoli.fa. It's very fast for H1N1, but that's not too surprising. It's slow right now on fastq creation. I already know what I need to do to fix that part, I put an [issue](https://github.com/ncsa/rusty-neat/issues/65) up here that I am working on in develop. But for now, rusty-neat functions and can be executed on some arbitrary data. Will update the readme as well. Once I get the fastq generation idea to work (will require changes to fastq_tools and runner scripts), next step will be parallelizing the main contig loop in runner. That should give us some speedup. Also, we went from using way too much memory to using almost no memory. I'm not sure which is better. It would be nice to find a balance. I still think that a database could work. I'm wondering, maybe we index on position. Just use like a json-based database: {"contig": "chr1", "sequence": {"1":A, "2": C, etc}} We'd want something that had fast retrieval for writing. Something to investigate, maybe.
- The code is significantly different now, I read in all the old models from the original NEAT2.0, and stored them as default models. I used a lot of the same algorithms that they used in NEAT 2.0, figuring it was closer to rust than NEAT 3.X. 
- rusty-neat writes out a ton of files. It should clean them up afterward, but I'm not sure it does yet, may need to investigate if there is an added command required. But it stores them in temp and linux cleans that up. Basically, we are splitting up the input fasta into smaller chunks (tbd on if there is an optimal size) then writing those to file for later retrival. Need to test various chunk sizes. I am basically doing blocked gzip-type things, but with some twists. If someone really understood bgzip format, it could probably be used directly. I'm not there yet.
- My runtimes running on a WSL2 Ubuntu setup with the latest rust installed: 311 milliseconds for H1N1.fa, 13 minutes for ecoli. It's all the IO during that step. I knew it was ugly as I was writing it but it's a rough draft. Obviously all the repeat code in fastq_tools is not ideal.
- The output files, so far have looked good. I realized I haven't tested paired-ended yet.
- I decided to map out the N-regions then overwrite them with random valid bases. The main reason being that it kept dipping into N regions and trying to use them for snp creation and that isn't going to work. NEAT 2.0 had a data structure that kept track of how long N-regions were, then would skip them if they were longer than a read length or overwrite them if they were short. That might work. Potential room for improvement, but it's fast so far.
- Eliminated bam creation for now. The option is in there, but it warns you that it doesn't work yet. Once fastqs are sorted and maybe after parallelism, it'll be next up.

9/11/2025
=========

## rusty-neat v1.1.0
- We have refactored the fastq code. Now a temporary fastq is written to file per block read in by the fasta reader. This should allow us to parallelize the temporary fastqs. Then it's simply a matter of concatenating them together. Works great! I got through single-ended E coli depth 10 in 36 seconds. I got through paired-ended E coli in 78 seconds. 
- Caveats: the fastqs are sorted by read start. Included in each name is the exact location where the sequence was drawn from. Now, when you run an alignment, the reads will be named in such a way that will have <chromosome name>_start-pos_length in the name, so you can shuffle the reads and run an alignment, and easily check to see if your alignment is in the correct contig and location.

9/13/2025
=========

## rust-neat v1.1.2
- I missed an intermediate version. Basically, it was converting to the command `neat gen-reads` instead of `generate-reads`, which will allow us to add subcommands, similar to how python NEAT worked. I also fixed some bugs in the fastq generation and paired-ended creation.
- In this release, we are adding a bed reader that can limit the fasta regions of interest. Basically, it will take an input bed and only generate reads over the sections in the bed.

9/19/2025
=========

## rneat v1.2.0
- Executable name change to `rneat`. In order to differentiate rusty-neat from the main neat, we are changing the binary name to `rneat` and we have expanded subcommands and made them more robust.
- Added bed filtering functions. Now if you supply a "filter_output" bed file to your gen-reads run, it will produce the original files plus a filtered version showing only reads and variants that overlap with the regions in the bed.
- It's also possible to run this utility separately with `rneat filter-reads`

5/9/2026
=========

## rneat v1.4.0
- Streaming FASTA reader (`FastaStream`) moved from `gen-gc-bias-model` to the `common` crate, making it available to all subcommands. `gen-reads` now streams the reference one contig at a time instead of writing temp JSON block files, eliminating up to ~570 MB of disk I/O for human-scale genomes and capping peak memory to a single contig rather than the full genome.
- Contig-level parallelism added to `gen-reads` via rayon `par_bridge`. When `produce_bam: false`, all contigs are processed concurrently using rayon's work-stealing thread pool; the output ordering is preserved by sorting results by contig index before writing. BAM output continues to use a sequential loop because coordinate-sorted writes must occur in reference order.
- New `num_threads` config key in `gen-reads`. Set to an integer to cap the rayon thread pool used for parallel contig processing. Omit (or set to `.`) to use all available cores (default). Set to `1` to disable parallelism entirely. Has no effect when `produce_bam: true`.
- `NeatRng::derive_child(idx)` added to the common RNG: each parallel contig task derives a deterministic per-contig seed from the parent RNG state and the contig index, so results are fully reproducible regardless of thread scheduling.

5/3/2026
=========

## rneat v1.3.0
- Added gen_mut_model mimicking the functionality in python Neat. The purpose is to allow users to create custom models based on real data.
- Added `gen-seq-error-model` subcommand. Reads FASTQ (.fastq or .fastq.gz) input and outputs a serialized `SequencingErrorModel` + `QualityScoreModel` (JSON, gzip-compressed). Supports optional BAM file input to infer a SNP transition matrix from real aligned reads using MD tags, and an optional TSV file to supply a fully custom transition matrix. Priority order: TSV override > BAM inference > default matrix.
- BAM-based transition matrix inference reads MD tags from aligned reads (via noodles), walking CIGAR + MD tokens to identify reference/alt base pairs at each SNP position and accumulating per-context counts into the 4×4 substitution matrix.
- Added FASTQ read shuffle in gen-reads: output reads are now randomly shuffled so alignment-tool ordering assumptions are not baked in. Simplified read naming to reduce verbosity.
- Added `template_config/gen_seq_error_model_template.yml` documenting all config fields with inline comments.
- Updated README with a full "Generating a Sequencing Error Model" section covering config schema, BAM MD-tag requirement, TSV format, and priority order.
- Bug fixes: degenerate zero-row handling in transition matrix normalization; edge-position SNP handling in MD-tag walker.
- Added `gen-frag-length-model` subcommand. Reads a paired-end BAM or SAM file and outputs a `FragmentLengthModel::Normal` (gzipped JSON) fitted to the observed template lengths. Filters using a median + 10×MAD outlier ceiling and a configurable per-length `min_reads` count floor (default 2; set 0 to disable). Added `template_config/gen_frag_length_model_template.yml` and README section.
- Implemented `produce_bam` / `output_bam` in `gen-reads` config. Setting `produce_bam: true` now writes a golden BAM file whose reads exactly match the generated FASTQs (same sequence with variants and sequencing errors applied, same quality scores). CIGAR strings reflect both variant positions and sequencing error indels: deletion errors contribute `D` ops and insertion errors contribute `I` ops, so the BAM alignment accurately represents the simulated read relative to the reference. The BAM is written in coordinate-sorted order using a windowed flush strategy (inspired by NEAT 2.1): records are buffered per block, flushed in sorted order at each block boundary, and any cross-block paired-end R2 reads are carried forward to the correct position. No `samtools sort` is needed; run `samtools index` directly on the output. Note: heterozygous variants use the same probabilistic application as the FASTQ, so reads drawn to the reference allele will not show the variant in either file.
- FASTQ shuffle is now global across all contigs rather than per-contig, so no chromosome ordering is visible in the output. All per-block temp files from every contig are collected into memory together and shuffled in a single pass before writing. For large genomes (> 500 Mbp) a runtime warning is emitted recommending `seqkit shuffle` as a memory-bounded post-processing alternative.
- Added `target_bed` config key to `gen-reads`. When supplied, reads and variants are generated only over regions in the BED; contigs absent from the BED are skipped entirely (no blocks read, no temp files written). This is the recommended approach for exome or targeted-panel runs on large genomes, where post-hoc `filter-reads` would otherwise process the full genome first. Both `.bed` and `.bed.gz` inputs are accepted.
- Added `input_vcf` config key to `gen-reads`. A single-sample VCF (`.vcf` or `.vcf.gz`) can now be supplied to force specific variants into the simulation. Each record must carry a `GT` field; SNPs and indels are applied at the given position before random variant generation runs on the remaining mutable bases. Complex variants (multi-base REF and ALT) are skipped with a warning. Random variants still fire at `mutation_rate` for positions not covered by the input VCF; set `mutation_rate: 0` to suppress them entirely. Caveats: only the first ALT allele of multi-allelic records is used; REF alleles are not verified against the reference sequence.
- Added `gen-gc-bias-model` subcommand. Reads a reference FASTA and a per-base coverage file (samtools depth, bedtools genomecov -d, or bedtools genomecov -dz), tiles the reference in sliding windows, and accumulates mean coverage per integer GC% bin to produce a `GcBiasModel` (gzipped JSON). Bins with fewer than `min_windows_per_bin` observations receive a neutral weight of 1.0 so sparse GC% values do not distort the model. Designed to scale to large genomes: coverage is indexed once and loaded one contig at a time; the reference is streamed without temp files; both GC counting and coverage summing use O(1)-amortized sliding windows so peak memory is bounded by the longest contig rather than total genome size. Added `template_config/gen_gc_bias_model.yml` documenting all config fields.
- Added GC bias simulation to `gen-reads`. When a `gc_bias_model` is supplied, fragment start positions are sampled directly from a per-position probability distribution derived from GC content, replacing the previous rejection-sampling approach. The model's `window_size` is embedded in the model file so the application-time window always matches the training-time window. Set `gc_bias_normalize_coverage: true` (default) to inflate fragment count to compensate for low-weight regions; set `false` to preserve the raw requested coverage.