# Delegation and Cluster Playbook

Operational reference for how work is routed in this repository. Load this when
delegation or cluster work starts. The durable rule lives in `CLAUDE.md` /
`AGENTS.md` ("Delegation And Cost Discipline"); this file is the how-to.

Principle: route each piece of work to the cheapest agent that can do it well,
even when the lead model is Opus 4.8 or Fable. The lead's scarce budget goes to
economics, identification, calibration judgment, specification, and review.

## Routing table

| Work | Owner | How |
|---|---|---|
| Economics, identification, calibration judgment, spec, final review, decisions | **Lead** (Opus 4.8 / Fable) | — |
| Codebase search & grounding | **Explore** agents | dispatch read-only agents in parallel; keep the conclusion, not the file dumps |
| Heavy or independent coding, second-opinion diagnosis, stuck-point rescue | **Codex (GPT-5.4)** | `codex:rescue` skill / `codex:codex-rescue` agent |
| Bulk mechanical edits, log scraping, table building, many independent tasks | **Haiku / Sonnet / Fable** subagents | `Agent` with a `model` override; structured output where possible |
| Any run > ~30 min, calibration sweeps, grid / DE searches | **Torch cluster** | `code/cluster/torch.sh`, smoke-tested first |
| Model-critical numerics (SMM objective, solver, `(1 - phi)` DP threshold, KFE, target measurement) | **Lead writes, or lead verifies against the math** | never blindly accept a delegated diff here |

Balanced routing, not blind delegation: the lead never hands off identification
or calibration judgment, and never reports a delegated calibration conclusion
without checking the target-fit table and identification against `CLAUDE.md`.

## When to reach for each agent

- **Explore** — "where is X", "which files touch Y", "does convention Z exist".
  Broad fan-out reads. Cheap, parallel, read-only. Prefer over doing the sweep
  yourself when it means reading across many files.
- **Codex (`codex:rescue`)** — a substantial or well-boxed coding task, a second
  independent implementation, a deeper root-cause pass when you are stuck, or a
  diagnosis you want cross-checked. Scope the handoff tightly: files, goal,
  constraints, and what "done" means.
- **Cheap-model subagents** — mechanical, high-volume, or embarrassingly
  parallel work: renaming, reformatting, scraping logs into a table, running the
  same check across many inputs. Use `dispatching-parallel-agents` for 2+
  independent tasks.
- **Cluster** — see below. Anything long or search-shaped.
- **Lead keeps** — the model spec, the objective, identification, and the final
  read of any calibration result.

## Torch cluster SOP

Prereqs: the `torch` SSH alias maps to `login.torch.hpc.nyu.edu` as `td2248`.
Account rule: always `torch_pr_570_general`, never `torch_pr_571_general` (it is
baked into the submit scripts). Working dir on Torch:
`$SCRATCH/projects/Fertility_Spring26/code/cluster`, unless `CALIBRATION_STATUS.md`
names a dated scratch snapshot — then use that exact path.

Use `code/cluster/torch.sh` as the default interface. It runs locally, shells
over `ssh torch`, encodes the account rule and working dir, and fails loudly if
the SSH probe fails (it never works around credentials).

```bash
# 1. auth probe + queue (do this first every session)
code/cluster/torch.sh status

# 2. ALWAYS smoke-test the exact loop before a real run
code/cluster/torch.sh smoke intergen-de smoke_$(date +%Y%m%d)      # tiny array, short budget

# 3. real submit — explicit run tag, walltime, array size, seed base
code/cluster/torch.sh submit intergen-de intergen_<run_tag> --array=1-40%40 --time=3:00:00

# 4. health: confirm best.json + cases.jsonl are being written
code/cluster/torch.sh health intergen-de intergen_<run_tag>

# 5. monitor
code/cluster/torch.sh status
code/cluster/torch.sh logs intergen-de              # last N lines of recent logs
code/cluster/torch.sh logs intergen-de --follow     # tail -f newest

# 6. collect when finished
code/cluster/torch.sh collect intergen-de intergen_<run_tag>
```

Families: `intergen-twohour`, `intergen-de`, `intergen-polish` (needs
`INTERGEN_SEED_THETA_JSON`), `direct`. Point at a snapshot with
`TORCH_REMOTE_DIR=/scratch/td2248/projects/<SNAPSHOT>/code/cluster`. Preview any
command without executing via `TORCH_DRYRUN=1`. Extra env goes through
`TORCH_ENV="INTERGEN_J=20 INTERGEN_NB=80"`.

Long-run safety (`CLAUDE.md` "Long-Run Search Safety") is mandatory: estimate run
size, smoke the loop (not one solve), define time/round/stop budgets, write
progress every case or 5 minutes, and treat 30 minutes with no heartbeat as
unhealthy. Do not launch a long search "to see what happens."

**Volatile — always confirm in `CALIBRATION_STATUS.md`, never assume:** the
active target set, the run-tag convention, the lifecycle grid (`J`), the current
seed/best `theta`, and any dated scratch snapshot path. The submit-script env-var
defaults (target set, `J`, seeds) drift; the READMEs may lag. `torch.sh` passes
your explicit tag and forwards overrides but does not pin these — that is on you
per run.

### Raw commands (reference, when torch.sh does not cover a case)

Submit scripts read env vars and default the SBATCH account/array/walltime.
Collectors run from `code/model`:

```bash
# intergen (twohour / de / polish share results dir + collector)
python tools/collect_intergen_panel_results.py \
  --results-dir ../cluster/results_intergen_housing_fertility_<RUN_TAG>

# direct geometry
python tools/collect_direct_geometry_results.py \
  --results-dir ../cluster/results_python_direct_geometry_<RUN_TAG>
```

Shutdown snapshot / local pull helpers already exist:
`collect_intergen_shutdown_snapshot.sh <SNAPSHOT_TAG>` and
`pull_intergen_shutdown_snapshots_local.sh <TARGET_TIME> [ATTEMPTS] [SLEEP]`.

## Codex handoff

`codex:rescue` delegates to Codex (GPT-5.4). Reach for it when a coding or
diagnosis task is substantial, when you want an independent second
implementation, or when a bug resists the first pass. Give it: the exact files,
the goal, the constraints (do not touch the objective/solver unless asked), and
the acceptance check. Then **verify its output** — run the check, and for
anything touching model-critical numerics, read the diff against the math. A
plausible-but-wrong Codex diff to the objective is exactly the failure mode the
lead exists to catch.

## Verify delegated output

Nothing delegated is "done" until the lead has verified it
(`verification-before-completion`): run the smallest check that proves the change,
read the diff on anything model-critical, and confirm calibration claims against
the full target-fit table and identification. Report what was and was not run.
