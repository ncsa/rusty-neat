# rusty-neat (rneat) — working practices for agents

Hand-written guidance below; the GitNexus block that follows is auto-generated
(regenerated between its `gitnexus` markers — keep edits **above** the markers).

## Before you touch anything
- **Sync first.** `git fetch --all` and confirm the branch is current before working.
  `develop` is the integration branch (PRs target `develop`, not `main`) and has
  silently diverged from `main` before — don't build on a stale checkout.
- **After a PR merges, verify your commits actually landed.** Late pushes to an open
  PR have missed the merge more than once. Check with
  `git merge-base --is-ancestor <sha> origin/develop`; recover a missed commit with
  `git cherry-pick`.

## Don't trust a green result — read the artifact
- Harness "overall: PASS" / summary lines can be **false passes** — e.g. a run that
  found no inputs still printed PASS over an empty `summary.tsv` (job 19887446). Open
  the actual output (per-item rows, counts, sizes) before believing a success.

## Delta / HPC (`scripts/delta/`)
- Real-data validation runs on **NCSA Delta** (SLURM, account `bhrd-delta-cpu`). The
  cluster filesystem is **not reachable from this workstation** — the user runs jobs
  and pastes results; artifacts archive to `/projects/bhrd/jallen17/rneat-access-results/`.
- Staging: `fetch_validation_data.sh` (references), `stage_soy.sh` (align + call a
  self-consistent ref/BAM/VCF; `FULL_GENOME=1` for the whole-genome stress vs the
  fast single-chromosome default). `model_builders.sbatch` exercises the builders.
- **Shell footgun:** `set -euo pipefail` + `zcat … | head` makes `zcat` take SIGPIPE
  (exit 141) when `head` closes the pipe early, and `set -e` then aborts a step that
  actually succeeded. Wrap head-truncated pipelines in `set +o pipefail` … `set -o pipefail`.

## Testing (Rust)
- `cargo test` builds a **fresh** binary (`assert_cmd::cargo_bin`) — no staleness worry.
  Integration tests live in `rneat/tests/` (`mod common;` for shared helpers).
- `model_parity.rs` pins builder output byte-for-byte on the H1N1 fixture
  (`BLESS_BASELINES=1` to regenerate after an intentional change).
- **Model fidelity** (built model file actually shapes gen-reads output) is covered by
  `model_fragment_fidelity.rs` and `model_output_fidelity.rs`; `docs/model_builder_baseline.md`
  records the Delta resource envelope and the fidelity status.

## Git / GitHub mechanics
- **`gh pr edit` is broken on this repo** (deprecated Projects-classic GraphQL). Patch a
  PR with `gh api -X PATCH repos/ncsa/rusty-neat/pulls/<N> -f body=…`. `gh issue edit` works.
- End commit messages with `Co-Authored-By: Claude <noreply@anthropic.com>`.

## GitNexus block (below)
- Auto-generated code-intelligence MCP guidance. Genuinely useful for impact analysis
  when editing Rust **symbols**; treat its "MUST" language as scoped to symbol edits —
  it does not apply to shell scripts, docs, or test-only changes.

<!-- gitnexus:start -->
# GitNexus — Code Intelligence

This project is indexed by GitNexus as **rusty-neat** (2845 symbols, 7093 relationships, 241 execution flows). Use the GitNexus MCP tools to understand code, assess impact, and navigate safely.

> If any GitNexus tool warns the index is stale, run `npx gitnexus analyze` in terminal first.

## Always Do

- **MUST run impact analysis before editing any symbol.** Before modifying a function, class, or method, run `gitnexus_impact({target: "symbolName", direction: "upstream"})` and report the blast radius (direct callers, affected processes, risk level) to the user.
- **MUST run `gitnexus_detect_changes()` before committing** to verify your changes only affect expected symbols and execution flows.
- **MUST warn the user** if impact analysis returns HIGH or CRITICAL risk before proceeding with edits.
- When exploring unfamiliar code, use `gitnexus_query({query: "concept"})` to find execution flows instead of grepping. It returns process-grouped results ranked by relevance.
- When you need full context on a specific symbol — callers, callees, which execution flows it participates in — use `gitnexus_context({name: "symbolName"})`.

## Never Do

- NEVER edit a function, class, or method without first running `gitnexus_impact` on it.
- NEVER ignore HIGH or CRITICAL risk warnings from impact analysis.
- NEVER rename symbols with find-and-replace — use `gitnexus_rename` which understands the call graph.
- NEVER commit changes without running `gitnexus_detect_changes()` to check affected scope.

## Resources

| Resource | Use for |
|----------|---------|
| `gitnexus://repo/rusty-neat/context` | Codebase overview, check index freshness |
| `gitnexus://repo/rusty-neat/clusters` | All functional areas |
| `gitnexus://repo/rusty-neat/processes` | All execution flows |
| `gitnexus://repo/rusty-neat/process/{name}` | Step-by-step execution trace |

## CLI

| Task | Read this skill file |
|------|---------------------|
| Understand architecture / "How does X work?" | `.claude/skills/gitnexus/gitnexus-exploring/SKILL.md` |
| Blast radius / "What breaks if I change X?" | `.claude/skills/gitnexus/gitnexus-impact-analysis/SKILL.md` |
| Trace bugs / "Why is X failing?" | `.claude/skills/gitnexus/gitnexus-debugging/SKILL.md` |
| Rename / extract / split / refactor | `.claude/skills/gitnexus/gitnexus-refactoring/SKILL.md` |
| Tools, resources, schema reference | `.claude/skills/gitnexus/gitnexus-guide/SKILL.md` |
| Index, status, clean, wiki CLI commands | `.claude/skills/gitnexus/gitnexus-cli/SKILL.md` |

<!-- gitnexus:end -->