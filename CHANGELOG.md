5/20/2026
=========
## rneat v1.5.1
- Added a binned scoring feature. Activated by user input when forming the quality score model, this will force rneat to bin the quality scores of the input data (if they are already binned, this is trivial) to user-defined bins. When this model is used in gen_reads, it creates a binned output.
- Added a number of unit tests and integration tests.

## rneat v1.5.0

### Memory improvements in `gen-reads`
- Replaced the per-position `Vec<f64>` bias map (one `f64` per chromosome base) with a compact `Vec<(start, end, rate)>` segment list. Measured peak RSS for rice (373 Mbp, 12 chromosomes) dropped from 3.06 GB → 1.05 GB (~3× reduction); major page faults dropped from 109 K → 0. Genomes previously too large to run safely (miscanthus ~11.9 GB estimated) now fit comfortably (~5.2 GB estimated).
- `cover_dataset` (single-ended path): fragment pool reduced from `num_frags` entries (`O(coverage × contig_length)`) to a single-entry pool with a cyclic index. Eliminates up to ~23 MB of VecDeque allocation per large chromosome and millions of redundant RNG shuffle calls.

### `gen-mut-model` fixes for real-world VCF data
- VCF parser no longer panics on missing allele calls (`./1`, `./.`) in the GT field. Missing alleles are skipped during genotype determination; variants with no allele data are treated as homozygous alt.
- Multi-allelic ALT fields (comma-separated, e.g. `T,TAATGGAATGG`) are now skipped with a debug log rather than being stored as garbled Complex variants.
- Variants whose FORMAT or SAMPLE fields are unparseable now warn-and-skip rather than aborting the model build, making the parser robust to mixed-format real-world VCFs.
- Soft-masked FASTA bases (`Maskeda/c/g/t`) are now canonicalized to uppercase (`A/C/G/T`) before comparing against the VCF REF allele in the trinucleotide context pass. This previously caused a hard `BaseMismatch` error on any variant in a repeat-annotated region of a soft-masked reference (e.g. hg38). Genuine mismatches (different base, not just masking) now warn-and-skip instead of aborting.
- Optional `bed_file` key in `gen-mut-model` config now handled safely; previously a missing key triggered a HashMap index panic.
- Per-variant "Found genotype in vcf file" log message demoted from `INFO` to `DEBUG`, eliminating millions of log lines when building models from whole-genome VCFs (NA12878 WGS produced ~4 M `INFO` lines before this fix).
- Human NA12878 whole-genome VCF (hg38) now builds successfully in ~2:15 (wall clock), peak RSS 5.6 GB, yielding `mutation_rate: 0.001457`, `homozygous_frequency: 0.390`, variant distribution 86% SNP / 6.6% INS / 7.4% DEL.

5/9/2026
=========

## rneat v1.4.1
- `gen-reads` new feature: `long_reads` configuration flag (boolean, default `false`). When `true`, fragments sampled from the fragment length model that are shorter than `read_len` are accepted rather than discarded, and the resulting read is truncated to the actual fragment length. This reproduces the realistic read-length distribution seen in real long-read datasets (PacBio HiFi, Oxford Nanopore), where size selection is imperfect and a tail of shorter reads is always present. When `false` (the default, appropriate for Illumina short reads), fragments shorter than `read_len` continue to be rejected — for Illumina, short fragments produce adapter read-through rather than truncated reads, which is outside the scope of this simulator. The flag affects both the uniform-coverage fragment path and the GC-bias-weighted path, as well as the FASTQ and BAM writers, which compute an effective read length of `min(fragment_len, read_len)` per fragment when `long_reads: true`. Added to `template_config/gen_reads_template.yml` with inline documentation.

5/9/2026
=========

## rneat v1.4.0
- `gen-reads` bug fix: `generate_read` in `fastq_tools.rs` was rejecting fragments whose length equals `read_length` (`<=` guard); all single-ended reads were silently discarded. Changed to `<` so fragments of exactly `read_length` are accepted.
- `gen-reads` bug fix: `map_buffer` in `fasta_stream.rs` was treating soft-masked bases (`Maskeda/c/g/t`) as N-regions, excluding repeat-annotated but valid sequence from read generation. Only `N` and `X` now delimit true gap regions.
- `gen-reads` bug fix: `generate_fragments` was drawing fragment lengths from the fragment model for single-ended runs, causing severe under-coverage when the model mean exceeds `read_length`. Single-ended fragments are now fixed at exactly `read_length`; the fragment model applies only to paired-ended runs.
- `gen-reads` bug fix: `cover_dataset` accumulated `start` across loop iterations instead of resetting to `temp_end`, causing reads to pile up at the beginning of each region rather than spanning it evenly.
- `gen-reads` / `mutation_model` bug fix: soft-masked nucleotides (`Maskeda/c/g/t`) were leaking into VCF REF/ALT fields and FASTQ read sequences as lowercase letters (`a/c/g/t`). `generate_mutation` now applies `check_base()` to the extracted `ref_base` and to each base of the deletion reference slice, unmasking all bases before they are written to output.
- New integration tests in `gen-reads`: `test_paired_ended_bam_flags_correct` (SAM flags on every R1/R2 record), `test_paired_ended_insert_size_matches_model` (TLEN mean within ±20% of configured fragment mean), `test_input_vcf_snp_appears_in_bam_reads` (seeded input-VCF SNP visible in BAM reads at the correct position).
- `gen-seq-error-model`: FASTQ input is now streamed one record at a time instead of being fully loaded into memory. Peak memory is bounded by a single record rather than the entire file, which matters for large (multi-GB) FASTQ inputs. The `max_reads` early-exit now also fires without reading the rest of the file.
- `gen-mut-model`: Trinucleotide probability computation is now O(n) instead of O(n²). Transition counts are pre-grouped by reference frame once before the probability loop, eliminating a redundant linear scan per trinucleotide context.
- `filter-reads` (`filter_fastq` and `filter_vcf`): Replaced `format!("{}\n", line)` heap allocations on every written line with two `write_all` calls (bytes + newline). Eliminated double HashMap lookups (`contains_key` + index) in favour of a single `get` call. Contig name reconstruction changed from a manual push loop to `join("_")`. The read-name prefix check simplified from `find` + index comparison to `starts_with`.
- `filter-reads` config: Removed a redundant identity `match` on a bool.
- `gen-frag-length-model`: `filter_lengths` no longer clones its input `Vec` before sorting — the owned value is sorted in place, eliminating one allocation proportional to the fragment count.
- `gen-gc-bias-model`: `overall_mean` computation replaced two separate filtered iterations with a single `fold`, scanning the 101-bin array once instead of twice.
- Streaming FASTA reader (`FastaStream`) moved from `gen-gc-bias-model` to the `common` crate, making it available to all subcommands. `gen-reads` now streams the reference one contig at a time instead of writing temp JSON block files, eliminating up to ~570 MB of disk I/O for human-scale genomes and capping peak memory to a single contig rather than the full genome.
- Contig-level parallelism added to `gen-reads` via rayon `par_bridge`. All contigs are now processed concurrently using rayon's work-stealing thread pool regardless of whether BAM output is requested; the output ordering is preserved by sorting results by contig index before writing.
- BAM creation is now fully parallel. Each contig worker writes its alignment records to a private temp BAM body file (no header, plain BGZF blocks). After all contigs finish, a single concatenation pass strips the BGZF EOF block from each intermediate file and assembles them in reference order into the final coordinate-sorted BAM. No `samtools sort` step is required.
- New `num_threads` config key in `gen-reads`. Set to an integer to cap the rayon thread pool used for parallel contig processing. Omit (or set to `.`) to use all available cores (default). Set to `1` to disable parallelism entirely.
- `NeatRng::derive_child(idx)` added to the common RNG: each parallel contig task derives a deterministic per-contig seed from the parent RNG state and the contig index, so results are fully reproducible regardless of thread scheduling.
- Removed `fasta_reader` and `fasta_map` modules. `SequenceBlock`, `SequenceMap`, and `RegionType` now live in `common::structs::sequence_block`; `map_buffer` and `apply_n_substitution` moved to `common::file_tools::fasta_stream`. Eliminates the temp JSON block files that were written to disk for every contig.
- `gen-mut-model` migrated to a single-pass `FastaStream` loop: trinucleotide counting and variant processing now share one iteration over each contig's sequence instead of two separate passes.

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

9/19/2025
=========

## rneat v1.2.0
- Executable name change to `rneat`. In order to differentiate rusty-neat from the main neat, we are changing the binary name to `rneat` and we have expanded subcommands and made them more robust.
- Added bed filtering functions. Now if you supply a "filter_output" bed file to your gen-reads run, it will produce the original files plus a filtered version showing only reads and variants that overlap with the regions in the bed.
- It's also possible to run this utility separately with `rneat filter-reads`

9/13/2025
=========

## rust-neat v1.1.2
- I missed an intermediate version. Basically, it was converting to the command `neat gen-reads` instead of `generate-reads`, which will allow us to add subcommands, similar to how python NEAT worked. I also fixed some bugs in the fastq generation and paired-ended creation.
- In this release, we are adding a bed reader that can limit the fasta regions of interest. Basically, it will take an input bed and only generate reads over the sections in the bed.

9/11/2025
=========

## rusty-neat v1.1.0
- We have refactored the fastq code. Now a temporary fastq is written to file per block read in by the fasta reader. This should allow us to parallelize the temporary fastqs. Then it's simply a matter of concatenating them together. Works great! I got through single-ended E coli depth 10 in 36 seconds. I got through paired-ended E coli in 78 seconds.
- Caveats: the fastqs are sorted by read start. Included in each name is the exact location where the sequence was drawn from. Now, when you run an alignment, the reads will be named in such a way that will have <chromosome name>_start-pos_length in the name, so you can shuffle the reads and run an alignment, and easily check to see if your alignment is in the correct contig and location.

9/9/2025
=========

## rusty-neat v1.0.0
- It's been a couple years since I had an update to rusty-neat. This new version is version 1.0.0! It is not yet fully functional, as it is missing the BAM creation, though it does work on some test data: H1N1.fa and ecoli.fa. It's very fast for H1N1, but that's not too surprising. It's slow right now on fastq creation. I already know what I need to do to fix that part, I put an [issue](https://github.com/ncsa/rusty-neat/issues/65) up here that I am working on in develop. But for now, rusty-neat functions and can be executed on some arbitrary data. Will update the readme as well. Once I get the fastq generation idea to work (will require changes to fastq_tools and runner scripts), next step will be parallelizing the main contig loop in runner. That should give us some speedup. Also, we went from using way too much memory to using almost no memory. I'm not sure which is better. It would be nice to find a balance. I still think that a database could work. I'm wondering, maybe we index on position. Just use like a json-based database: {"contig": "chr1", "sequence": {"1":A, "2": C, etc}} We'd want something that had fast retrieval for writing. Something to investigate, maybe.
- The code is significantly different now, I read in all the old models from the original NEAT2.0, and stored them as default models. I used a lot of the same algorithms that they used in NEAT 2.0, figuring it was closer to rust than NEAT 3.X.
- rusty-neat writes out a ton of files. It should clean them up afterward, but I'm not sure it does yet, may need to investigate if there is an added command required. But it stores them in temp and linux cleans that up. Basically, we are splitting up the input fasta into smaller chunks (tbd on if there is an optimal size) then writing those to file for later retrival. Need to test various chunk sizes. I am basically doing blocked gzip-type things, but with some twists. If someone really understood bgzip format, it could probably be used directly. I'm not there yet.
- My runtimes running on a WSL2 Ubuntu setup with the latest rust installed: 311 milliseconds for H1N1.fa, 13 minutes for ecoli. It's all the IO during that step. I knew it was ugly as I was writing it but it's a rough draft. Obviously all the repeat code in fastq_tools is not ideal.
- The output files, so far have looked good. I realized I haven't tested paired-ended yet.
- I decided to map out the N-regions then overwrite them with random valid bases. The main reason being that it kept dipping into N regions and trying to use them for snp creation and that isn't going to work. NEAT 2.0 had a data structure that kept track of how long N-regions were, then would skip them if they were longer than a read length or overwrite them if they were short. That might work. Potential room for improvement, but it's fast so far.
- Eliminated bam creation for now. The option is in there, but it warns you that it doesn't work yet. Once fastqs are sorted and maybe after parallelism, it'll be next up.

9/11/2023
=========

## rusty-neat v0.1.0

Moved the code from a branch under NEAT (https://github.com/ncsa/neat) to its own repo https://github.com/ncsa/rusty-neat. This project will be considered a separate, related project. We will use the same underlying design ideas as NEAT, but it will essentially be different and won't be guaranteed in any way to be able to reproduce what the original NEAT did (though we aim to get it as close as possible, the differences in Python and Rust would make byte-by-byte compatibility very difficult). I have implemented an idea I prototyped in a private repo (https://github.com/joshfactorial/basic_neat), of reading in the file in what I hope is a straightforward way. I was looking for stuff like running out of memory and speed. It can read in the `E. coli` fasta file in about 30 seconds on my laptop, which is pretty fast, compared to Python. What's more is it can store the whole thing in memory, due to the efficient way Rust uses memory. Instead of copying the entire data set every time in memory, quickly eating it up, Rust forces you to decide when to pass a read-only version, copy the dataset, or use up the dataset once and lose it forever (all valid choices depending on context). So I'm feeling hopeful that this will work. I am publishing an Alpha called v0.1.0 for people to test, if they feel so inclined. See the README for more details.

## rusty-neat v0.1

Goal: To release a version of NEAT in rust capable of reading an arbitrary fasta file,
and outputting a version of that file with 1% of the bases mutated as SNPs.

9/5/2023
=========
- Initialized Cargo package with hello world default script.
- Added function to select a random integer from a list (for selecting positions)

9/9/2023
=========
- Created a first pass program that mutates a hardcoded DNA string at one random non-N position.
- Added a logging function to print out time-stamped log items both to stdout and to file.
- Need to improve error handling all around, taking advantage of Rusts error handling system.
- Design decision: After toying with some options to create a standardized RNG to use throughout the program, it seems the best option is to use an on-the-fly rng per function for now. This eliminates reproducibility between runs, but should enable parallelism of processing.
- Redesign from the ground up. Trying to use KISS as much as possible, and only introduce complexity as it becomes necessary, or as I learn Rust tricks that may help this endeavor.
