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
| **built `fragment_model` file → BAM insert size** | `tests/model_fragment_fidelity.rs` — fit on a 200 bp-mean BAM, simulate with the file, output insert mean = **199.8 bp** (fitted 200.0) |
| **built `sequence_error_model` file → read qualities** | `tests/model_output_fidelity.rs` — train on a uniform-Phred-35 FASTQ, simulate, output mean quality = **35.0** |
| **built `mutation_model` file → output variant count** | `tests/model_output_fidelity.rs` — a 139× rate difference between two fitted models yields 0 vs **405** output variants (405 ≈ rate 3.17e-2 × 13,133 bp) |
| **built `gc_bias_model` file → which regions get sequenced** | `tests/model_output_fidelity.rs` — train on a BAM covering a 20%-GC contig but not an 80%-GC one (fitted `w[20]=1.98`, `w[80]=0.02`), simulate: output is **506 low-GC reads / 0 high-GC** |
| input VCF variant → output reads | `gen_reads … test_input_vcf_snp_appears_in_bam_reads` |

The `model_builders.sbatch` round-trip is only a *usability* smoke test (reads exist
and are correctly-lengthed); the fidelity tests above prove the build→simulate path
is numerically faithful through the model-file load path — for **all four** builders.

Because H1N1's ~40% GC can't populate distinct GC bins, the gc_bias test builds a
synthetic 2-contig reference (a 20%-GC and an 80%-GC contig) so the weight table is
deterministic; the others reuse the bundled H1N1 fixture.
