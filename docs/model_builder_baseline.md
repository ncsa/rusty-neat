# Model-builder resource baseline & output fidelity

The model builders (`gen-seq-error-model`, `gen-frag-length-model`, `gen-gc-bias-model`,
`gen-mut-model`, `gen-bam-models`) are what a first-contact user runs on *their own*
data. `rneat/tests/model_parity.rs` pins the fit algorithms on tiny fixtures; this
document records how they behave on real, full-size input — the resource envelope
and the proof that a fitted model actually shapes gen-reads output.

Harness: `scripts/delta/model_builders.sbatch`, fed by `scripts/delta/stage_soy.sh`.

## Resource envelope (Delta, soybean Wm82, ~978 Mb, ~200 contigs)

Two staging modes, same builders. Chr mode scopes to the largest contig (fast
correctness check); `FULL_GENOME=1` stages the whole messy assembly (the real
stress). Numbers are wall-clock and peak RSS from `/usr/bin/time -v`.

| builder      | chr wall / RSS | full-genome wall / RSS | notes |
|--------------|----------------|------------------------|-------|
| seq_error    | 0:06.7 / 19 MB | 0:53.8 / **18 MB**     | streams the FASTQ — memory-flat at any scale |
| frag_length  | 0:03.3 / 50 MB | 0:42.0 / 748 MB        | BAM TLEN scan |
| gc_bias      | 0:04.5 / 323 MB| 1:03.5 / **3.6 GB**    | loads reference + per-GC windows across all contigs |
| mut_model    | 0:03.4 / 228 MB| 0:39.9 / 1.3 GB        | reference + genome-wide VCF |
| bam_models   | 0:04.5 / 352 MB| 1:04.3 / **4.2 GB**    | frag + GC in one BAM pass — heaviest |
| roundtrip    | 177,275 reads  | 183,583 reads          | gen-reads consumes the built models, emits valid pairs |

Source runs: chr = job 19887967, full-genome = job 19888625. Both `overall: PASS`.

### Memory scales with reference size — plan accordingly

The heavy builders (`gc_bias`, `bam_models`) hold the reference plus per-window
state, so peak RSS tracks **genome size**, not read count:

- ~4.2 GB peak at ~1 Gb (soybean).
- Extrapolated **~12–14 GB for a 3.1 Gb human genome.**

Trivial on an HPC node, but a user building `bam_models`/`gc_bias` against a human
reference on a 16 GB laptop is near the edge. `seq_error` is the exception — it
streams and stays flat (~18 MB) regardless of input size.

Nothing approached the harness's 128 GB / 8 h request at soybean scale; no OOM, no
runtime cliff, clean on a ~200-contig assembly.

## Output fidelity — do the fitted numbers reach the reads?

A model that parses is not the same as a model that shapes the output. Verified
links (all are real, runnable tests, not assertions on paper):

| model → output | evidence |
|----------------|----------|
| explicit `fragment_mean` → BAM insert size (TLEN) | `gen_reads … test_paired_ended_insert_size_matches_model` |
| **built `fragment_model` file → BAM insert size** | `tests/model_fragment_fidelity.rs` — fit a model on a 200 bp-mean BAM, simulate with the file, output insert mean = **199.8 bp** (fitted 200.0) |
| input VCF variant → output reads | `gen_reads … test_input_vcf_snp_appears_in_bam_reads` |

The `model_builders.sbatch` round-trip is only a *usability* smoke test (reads exist
and are correctly-lengthed); the fragment fidelity test above is what proves the
build→simulate path is numerically faithful through the model-file load path.

### Not yet numerically verified (follow-ups)

The same "built model file drives output" proof is still open for the other three
builders — worth closing the same way:

- **seq_error**: built error/quality profile → output FASTQ quality distribution.
- **gc_bias**: built bias curve → output coverage-vs-GC.
- **mut_model**: built mutation rate/spectrum → output VCF variant rate.
