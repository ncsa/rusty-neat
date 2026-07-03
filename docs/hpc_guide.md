# Running rneat on an HPC cluster

This guide is for users who want to run `rneat` at scale — from a bacterial genome
on one node up to a whole human (or larger) genome across many nodes — on a SLURM
cluster. It distills what we learned validating rneat on NCSA Delta (see
`docs/access_report_draft.md` §3.6 for the underlying measurements) into
cluster-agnostic guidance.

The SLURM scripts referenced here live in `scripts/delta/`. They carry Delta
specifics (module names, an `--account`, conda env names); adapt those three
things to your site and the rest is portable.

---

## 1. The one-paragraph mental model

rneat is **single-thread efficient and memory-bandwidth-bound**, with **low, flat
memory**. In-process (rayon) threading is *break-even at best* — a single rneat
process already saturates a surprising amount of a node's memory bandwidth, so
adding threads to one process buys little and can even regress. The way to use a
big machine is therefore **not** "one process, many threads" but **many
single-threaded processes**, each simulating a different region of the genome,
spread across nodes. This is *region sharding*, and it is the entire HPC story.

Three practical consequences follow, and the rest of this guide is just their
detail:

1. **Run one thread per process** (`num_threads: 1`). Pack multiple processes onto
   a node instead of threading one.
2. **Run on exclusive nodes.** rneat competes for memory bandwidth; a co-tenant's
   bandwidth-hungry job on a shared node can inflate an 8-minute window to 170
   minutes with no warning (measured — see §5).
3. **Match packing density to node availability**, not to core count. More
   processes per node ≠ proportionally more throughput past ~2 (the bandwidth
   ceiling); it just makes each one slower.

---

## 2. Before anything at scale: the log-level footgun

**Set `--log-level info` (the default as of v1.19.1) and never run a large job at
`debug` or `trace`.**

rneat's per-base `debug!`/`trace!` events fire on the order of
`coverage × read_length × reference_bp` times. At `trace`, a single chromosome-22
run wrote a multi-gigabyte `.neat.log` and spent *most of its wall-clock on log
I/O*. Before v1.19.1 the default was `trace`; if you are on an older build, pass
`--log-level info` explicitly on every run. On any version, only turn up verbosity
for a *small* debugging run, never a production one.

```bash
rneat --log-level info gen-reads -c sim.yml    # good for production
# rneat --log-level debug ...                  # ONLY on a tiny debugging input
```

---

## 3. Quickstart: a small genome on one node

For anything that fits comfortably in one node's memory and finishes in a
reasonable single-process wall-clock (bacteria, fungi, small invertebrates), you
do not need sharding at all. One single-threaded process:

```yaml
# sim.yml
reference: /path/to/genome.fa
read_len: 151
coverage: 30
ploidy: 2
paired_ended: true
fragment_mean: 350
fragment_st_dev: 50
produce_fastq: true
produce_vcf: true
produce_bam: false
overwrite_output: true
output_dir: /scratch/$USER/run1
output_filename: sample
rng_seed: my run 1
num_threads: 1
```

```bash
#SBATCH --partition=<your-cpu-partition>
#SBATCH --nodes=1
#SBATCH --exclusive          # or a generous --cpus-per-task + --mem if you can't get exclusive
#SBATCH --time=01:00:00
rneat --log-level info gen-reads -c sim.yml
```

If you want more throughput on that node, do not raise `num_threads` — instead run
several of these processes concurrently on different regions (see the sharding
approach below), or several independent samples at once.

**Memory:** rneat's resident memory is dominated by the reference contigs it loads,
not by coverage. Reference loading is `target_bed`-aware and loads only the contigs
a run touches, so peak RSS per process ≈ *the size of the largest contig that
process simulates* plus a small overhead. Budget accordingly: a 25 Mb window on
human chr1 still loads all of chr1 (~250 MB); a whole small genome loads all of it.

---

## 4. Whole-genome scale: region-sharded array jobs

The genome is partitioned into contiguous **anchor-windows**. rneat's `target_bed`
is a *generation-time* filter in absolute reference coordinates that anchors each
fragment to exactly one window (reads may extend past a window edge; ownership is
by anchor point), so the union of all shards reconstructs a whole-genome run **with
no double-counting and no gaps**, and the per-shard outputs **merge by simple
concatenation**. Each window is one single-threaded rneat process; a SLURM array
spreads them across nodes.

The tooling (in `scripts/delta/`):

| Script | Role |
|---|---|
| `make_shard_beds.sh` | partition `REF.fa.fai` into `shard_NNNN.bed` windows; prints the shard count |
| `genome_array.sbatch` | SLURM array; each task owns one **exclusive node** and runs `SHARDS_PER_NODE` windows on it concurrently. Idempotent + `--requeue`. |
| `run_genome_reps.sh` | driver: sweep packing density `K` and/or run replicates with distinct seeds |
| `aggregate_reps.sh` | fold rep manifests into mean ± sd wall-clock tables |
| `merge_shards.sh` | concatenate per-shard FASTQ + `bcftools concat|sort` the truth VCF into one dataset |

**Do a fail-fast smoke test first** (see §7) so you never commit a big node
request to a config that dies in the first shard.

Then the full run:

```bash
# 1. index + partition (25 Mb windows → ~124 shards for GRCh38 primary assembly)
samtools faidx REF.fa
N=$(bash scripts/delta/make_shard_beds.sh REF.fa.fai 25000000 $SCRATCH/shards)

# 2. choose packing density K and derive node-tasks
K=8; T=$(( (N + K - 1) / K ))

# 3. submit the exclusive-node array (%20 caps concurrent nodes — tune to availability)
SHARD_DIR=$SCRATCH/shards OUTROOT=$SCRATCH/wg_run REFERENCE=REF.fa SHARDS_PER_NODE=$K \
  sbatch --array=0-$((T-1))%20 scripts/delta/genome_array.sbatch

# 4. after the array finishes, merge (on a compute node — it's I/O-heavy)
SHARD_OUTROOT=$SCRATCH/wg_run sbatch scripts/delta/merge_shards.sh
```

For a **simulation-only timing** run (no downstream alignment/calling — the number
users actually want to plan around) this is the whole job: steps 1–3 produce the
FASTQ + truth VCF and the array's `time.txt` files hold the per-window wall-clock.

---

## 5. Choosing packing density (K) and whether to request a reservation

This is the one genuinely non-obvious decision, and it is driven by three
trade-offs we measured on GRCh38 (306 × 25 Mb shards):

- **Bandwidth saturates by K≈2.** A 25 Mb human window runs in ~8 min alone (K=1),
  **29 min at K=2, and 101 min at K=4** on the same node — superlinear slowdown.
  Packing more shards per node does *not* raise node throughput past ~2; it only
  makes each shard slower.
- **Low-K-many-nodes is the fastest compute but the hardest to schedule.**
  `--exclusive` claims a whole node even for a 1-core task, so minimizing
  wall-clock means grabbing dozens of whole nodes at once — which queues behind
  everyone on a busy cluster.
- **A tight per-task walltime is harmful once you're exclusive.** It silently kills
  legitimately-slow packed windows. Use a *generous* walltime + `--requeue`
  (`genome_array.sbatch` already does; `run_genome_reps.sh` sizes it `30 + 20K` min).

Measured whole-genome results (each run complete, 306/306 shards):

| Configuration | Whole-genome wall | Heaviest window | Notes |
|---|---|---|---|
| Shared partition (16-core slice) | 3 h 41 m | **170 min** | contention-bound, irreproducible — avoid |
| Exclusive, **K=2** | **2 h 30 m** | 29 min | ~30 nodes in waves; beat shared on every axis |
| Exclusive, **K=4** | 12 h 32 m | 101 min | slow = *scheduling* (only ~10 nodes free), not compute |

**Rule of thumb:**

- **Cluster busy, no reservation, and you can wait:** exclusive, **K=8–16**. Fewer
  whole-node grabs schedule sooner; each shard is slower but you're not paying for
  idle cores and the run completes unattended. This is the pragmatic default.
- **Cluster idle or you hold a reservation:** exclusive, **K=1–2** across as many
  nodes as you can get. This is the ~30-min compute floor — a whole human genome at
  30× in a single ~29-min wave when ~150 nodes are actually available.

**Do you need a reservation?** For most work, **no.** A reservation only buys the
low-K compute floor (the "30-minute genome" headline). If you have no hard
deadline and are fine letting a moderate-K run schedule and complete on its own
time, skip the reservation and run K=8–16 exclusive. Request a reservation only
when you specifically need many whole nodes *simultaneously* — i.e. you are chasing
the minimum wall-clock, not just a finished dataset.

**Cost intuition (SLURM node-hours):** a whole-genome 30× simulation is *cheap* —
on the order of tens of node-hours of actual compute (each 25 Mb window is
minutes). The binding constraint at scale is **per-run wall-clock and node
availability, not accounting credits.** Plan around scheduling, not budget.

---

## 6. Reproducibility scope

rneat is deterministic and parallelism-invariant *where it matters*, with one
documented boundary:

- **Determinism:** same seed + same reference + same params → byte-identical output.
- **Thread-invariant:** 1 thread vs N threads → byte-identical. Parallelism never
  perturbs results.
- **Shard-order independent + disjoint:** the sharded whole-genome path gives the
  same merged result regardless of shard completion order, with zero duplicate
  variants across windows. This is why sharding is safe.
- **Contig-order sensitive (by design, monolithic runs only):** a single
  *monolithic* (non-sharded) run consumes one global RNG stream in contig order, so
  reversing contig order in the FASTA changes the realization (even the variant
  count). **Consequence:** to reproduce a monolithic run exactly, use the
  *byte-identical* reference (same contig order), not merely the same genome. The
  sharded path is immune to this — each shard seeds independently from an
  absolute-coordinate window. (Tracked as #322: optional per-contig seeding to make
  even monolithic output contig-order-independent.)

Practical guidance: **for reproducible large runs, prefer the sharded path and
record your seed(s) + the exact reference file** (or its checksum).

---

## 7. Fail-fast before you commit big node requests

A large array that dies on a bad config or a malformed-output bug wastes both your
time and your queue priority. Two cheap guards:

1. **Smoke-test one window first.** `scripts/delta/wg_smoke.sh` runs a single small
   window as a normal (non-exclusive, few-minute) job and then *validates the output
   the way a downstream aligner would* — it checks that every FASTQ record has
   `len(seq) == len(qual)` and that reads and truth variants were actually produced.
   This is the exact class of defect (a malformed FASTQ that `wc -l` can't see but
   bwa-mem2 silently truncates on) that cost us a full at-scale run during adapter
   validation. Run it whenever you change the build, the config, or the reference:

   ```bash
   REFERENCE=$SCRATCH/neat_data/genome.fa bash scripts/delta/wg_smoke.sh
   ```

2. **Throttle the first real array** (`--array=...%N` with a small `N`) so a
   systematic problem surfaces on the first few nodes, not all of them at once.

`genome_array.sbatch` is idempotent (a completed shard is skipped on retry) and
uses `--requeue`, so a task killed by a sick node reruns elsewhere and only the
unfinished windows in its batch redo — you never lose completed work.

---

## 8. Shared-filesystem hygiene (Lustre / GPFS)

Large simulation runs move a lot of FASTQ/BAM. On a shared parallel filesystem:

- **Stage heavy I/O on node-local disk** (`$SLURM_TMPDIR` / local NVMe) when your
  site provides it, and copy only the small final artifacts back. A dropped Lustre
  OST otherwise wedges jobs in uninterruptible I/O wait. (`germline_e2e.sbatch`
  demonstrates this pattern.)
- **Keep manifests and bookkeeping on durable storage, never scratch** — scratch can
  throw transient EIO or hit quota mid-run, and you don't want a driver half-submit
  because a tiny write failed. (`run_genome_reps.sh` writes its manifest to `$HOME`.)
- **Prune intermediates.** Whole-genome 30× FASTQ is ~35–40 GB per run; delete BAMs
  and FASTQ after scoring if you only need the metrics (`PRUNE=1` in the harnesses).
- **Watch your scratch quota** — on Delta it is a *project* quota, not per-user, so a
  teammate's data counts against yours.

---

## 9. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Run is far slower than expected, huge `.neat.log` | `--log-level trace/debug` | use `--log-level info` (default ≥ v1.19.1) |
| Whole-genome wall-clock wildly variable, some windows 10–20× slow | shared-partition co-tenant contention | request `--exclusive` nodes |
| Array tasks stuck `PD (Resources/Priority)` forever | low-K wants too many whole nodes at once | raise K (8–16), or request a reservation |
| Packed windows all hit the walltime and finish incomplete | walltime ceiling below packed-window time | raise `--time` (never a tight ceiling once exclusive) |
| Downstream aligner recovers almost no reads despite full read count | malformed FASTQ record | run `wg_smoke.sh`'s integrity check; upgrade to ≥ v1.19.0 |
| Can't reproduce a previous run's variants | different reference contig order | use the byte-identical reference; or use the sharded path |
| Peak RSS higher than expected | a shard's largest contig is large | lower K so fewer big contigs load per node simultaneously |

---

*Measurements in this guide come from the rneat validation campaign on NCSA Delta;
see `docs/access_report_draft.md` (§3.6 for scaling, §3.10 for reproducibility,
§3.11 for the malformed-FASTQ story) and `docs/regression_protocol.md` for how the
performance and fidelity baselines are frozen and gated.*
