# Delegation and Cluster Playbook

Operational reference for how work is routed in this repository. Load this when
delegation or cluster work starts. The durable rule lives in `CLAUDE.md` /
`AGENTS.md` ("Delegation And Cost Discipline"); this file is the how-to.

Principle: route each piece of work to the least expensive adequate worker. The
lead's scarce budget goes to economics, identification, calibration judgment,
specification, and final review. Model-selected CLI workers handle bounded
supporting work; Torch handles long searches.

## Routing table

| Work | Owner | How |
|---|---|---|
| Economics, identification, calibration judgment, spec, final review, decisions | **Lead** | Keep central judgment in the main task |
| Codebase search, grounding, log parsing, table building | **`explorer_fast`** | Model-selected, low-reasoning, read-only CLI worker |
| Tiny one-shot mechanical work | **`mechanic_fast`** | Fastest profile; use only when deep validation is unnecessary |
| Scoped implementation or independent diagnosis | **`worker_fast`** | Model-selected CLI worker with exclusive file ownership; `--write` only when needed |
| Strong read-only second opinion | **`reviewer_strong`** | Higher-reasoning model-selected CLI worker |
| In-thread parallelism with shared context | **Built-in subagent** | Useful for coordination, but the spawn tool does not guarantee a cheaper model |
| Any run > ~30 min, calibration sweeps, grid / DE searches | **Torch cluster** | `code/cluster/torch.sh`, smoke-tested first |
| Model-critical numerics (SMM objective, solver, `(1 - phi)` DP threshold, KFE, target measurement) | **Lead writes, or lead verifies against the math** | never blindly accept a delegated diff here |

Balanced routing, not blind delegation: the lead never hands off identification
or calibration judgment, and never reports a delegated calibration conclusion
without checking the target-fit table and identification against `CLAUDE.md`.

## Autonomous dispatch protocol

The lead, not the user, selects the route. A user request needs an outcome and
constraints; deadline or urgency is helpful, but asking the user to name a
model/profile defeats the purpose of the routing system.

Classify by the judgment and risk the answer needs, not by prompt length. Use
the following default ladder; a task can move up when its evidence warrants it.

| Request | First route | Lead's acceptance check |
|---|---|---|
| Answer already in current context or a tiny one-step edit | Lead, briefly | Confirm the cited fact or diff; avoid worker setup overhead |
| Find files, extract facts, reduce logs, build a table | `explorer_fast`, bounded read-only | Open the decisive cited lines or source rows; reject unsupported conclusions |
| Small mechanical change | `mechanic_fast`, scoped write | Inspect the diff and run the smallest relevant check |
| Well-specified implementation or diagnosis | `worker_fast`, exclusive file ownership | Inspect the diff, run a targeted check, and examine edge cases |
| Conflicting evidence or a disputed technical claim | Lead, with `reviewer_strong` for a bounded second opinion when useful | Resolve the conflict from source evidence; state what remains unknown |
| Economic mechanism, identification, calibration decision, model contract, or model-critical numerics | Lead | Delegate extraction or implementation if useful; independently verify the scientific conclusion |

Give workers a stopping rule in the task prompt. They should return a short
answer with file/line evidence or a diff, checks performed, and uncertainty.
If the answer needs an unstated economic assumption, source conflicts, a wider
search, or a contract change, they stop and report it. The lead decides whether
to answer from the available evidence, send a newly scoped pass, or
ask the author about a genuine preference or scientific decision. Escalation
responds to observed difficulty; it is not an automatic second model call.

Bound search before dispatch: name the likely folder or source, the question to
settle, and what finding is sufficient to stop. Expand in bounded steps when
the initial evidence is insufficient, with each new scope stated explicitly. For
web questions, use current authoritative sources when the claim can change;
do not launch a broad crawl for a simple answer. Keep worker output compressed
to conclusions and decisive evidence rather than raw dumps.

1. Keep the work with the lead if the likely handoff cost exceeds the work.
2. Otherwise dispatch as many least-cost adequate workers as have useful
   independent tasks or distinct verification roles, up to the available agent
   limit. Tommaso authorizes this proactively; there is no one-worker default.
   Give each a narrow scope, required output, verification, time limit, and
   stop condition.
3. Give each worker a separate deliverable and each write-capable worker
   exclusive file ownership. Record how the lead will verify and combine the
   results before dispatching. Do not fill slots with duplicate broad inquiries
   or work that adds no value.
4. Before dispatching, tell the user: `route | reason | time limit | stop
   condition`. For work that lasts beyond one interaction, update in the form
   `phase | elapsed | route | artifact/evidence | next decision`.
5. When a time limit expires, synthesize or report the bounded uncertainty.
   Do not silently retry, expand, or fan out; a further attempt needs a changed
   hypothesis, scope, or method.

Urgent work defaults to the smallest route that can answer the question: a
short read-only worker pass or a narrowly verified edit. It never defaults to
a broad full-context audit or agent fan-out. Search-shaped numerical work goes
to Torch and follows the smoke/checkpoint/health procedure below.

## When to reach for each agent

- **`explorer_fast`** — "where is X", "which files touch Y", "does convention Z
  exist", or "reduce these logs to a table". Default route for broad read-only
  grounding. Require conclusions and file/line references, not file dumps.
- **`mechanic_fast`** — a tiny rename, formatting change, or deterministic
  extraction. Its speed comes with a shallow workflow, so never assign model
  specification, calibration interpretation, or numerically sensitive code.
- **`worker_fast`** — a well-boxed implementation, an independent diagnosis, or
  a second implementation. Give it exclusive files, an explicit acceptance
  check, and `--write` only when the task actually requires edits.
- **`reviewer_strong`** — a read-only adversarial pass on a proposed diff,
  numerical argument, or research artifact. The lead adjudicates the result.
- **Built-in subagents** — convenient when live steering, shared thread context,
  or direct collaboration matters. They are not the default cost-saving route
  because the spawn API does not expose a per-child model selector.
- **Cluster** — see below. Anything long or search-shaped.
- **Lead keeps** — the model spec, the objective, identification, and the final
  read of any calibration result.

## Model-selected Codex workers

Use `ops/codex-workers/scripts/codex-worker.sh` whenever explicit model and
reasoning selection is part of the routing decision. Current model IDs and
reasoning levels live in `ops/codex-workers/config/models.env`; do not repeat
them in prompts or durable agent instructions. The wrapper is ephemeral,
read-only by default, rejects workdirs outside the repository, and never commits
or pushes.

`explorer_fast` and `mechanic_fast` use the context-minimal exception in
`AGENTS.md`; this prevents a tiny task from spending most of its context on the
calibration memory stack. The exception is unavailable for economics,
calibration, targets, model results, or model-code interpretation. Add
`--full-context` whenever a minimal-profile task crosses that boundary.

Start every worker prompt from
`ops/codex-workers/prompts/task-template.md` and specify the goal, file scope,
files not to touch, compressed return format, verification, and stop condition.
The wrapper checks the route and permissions; the caller remains responsible
for keeping the prompt bounded and complete.
Use the available concurrency when independent tasks or distinct verification
roles benefit from it; do not impose a fixed one-, two-, or three-worker cap.
Give each write-capable worker exclusive ownership of its files.

```bash
# Prepare one bounded task.
cp ops/codex-workers/prompts/task-template.md /tmp/codex-task.md

# Read-only search or grounding.
ops/codex-workers/scripts/codex-worker.sh \
  --profile explorer_fast \
  --prompt-file /tmp/codex-task.md

# Scoped implementation. The prompt must give exclusive file ownership.
ops/codex-workers/scripts/codex-worker.sh \
  --profile worker_fast \
  --prompt-file /tmp/codex-task.md \
  --write \
  --workdir code/model

# Inspect routing without making a model call.
ops/codex-workers/scripts/codex-worker.sh \
  --profile reviewer_strong \
  --prompt-file /tmp/codex-task.md \
  --dry-run
```

Do not use a faster service tier merely to reduce wall time when conserving
usage is the objective. Prefer a lower-cost worker profile and bounded prompt;
reserve the lead and deeper reasoning for judgment that actually needs them.

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

## Verify delegated output

Nothing delegated is "done" until the lead has verified it
(`verification-before-completion`): run the smallest check that proves the change,
read the diff on anything model-critical, and confirm calibration claims against
the full target-fit table and identification. Report what was and was not run.
