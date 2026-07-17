# Design: Standardize Delegation, Cluster, and Best-Practice Discipline

Date: 2026-07-02
Status: historical design note; implementation is maintained in
`docs/workflow/delegation_and_cluster_playbook.md` and `ops/codex-workers/`

The model and agent names below reflect the July 2 design discussion. They are
kept for historical context; use the current named worker profiles in the
canonical playbook for live delegation.

## Purpose

Make three things the *default* every session, not ad hoc choices:

1. **Delegation.** Even when the lead model is Opus 4.8 or Fable, route search,
   bulk/mechanical coding, independent investigations, and second-opinion
   diagnosis to cheaper/faster agents (Explore, Codex/GPT-5.4, Haiku/Sonnet/Fable
   subagents). The lead spends its budget on economics, identification,
   calibration judgment, spec, and review.
2. **Cluster.** Long or search-shaped compute goes to Torch by a single,
   recoverable procedure — not a hand-written 40-line prompt each time.
3. **Verification discipline.** Delegated output that touches the model spec or
   the SMM objective is verified against the math by the lead before it is
   trusted.

This is a workflow/documentation change. No model code changes, no calibration
changes.

## The routing standard

| Work | Owner |
|---|---|
| Economics, identification, calibration judgment, spec, final review, decisions | **Lead** (Opus 4.8 / Fable) |
| Codebase search & grounding | **Explore** agents — parallel, read-only |
| Heavy/independent coding, second-opinion diagnosis, "stuck" rescue | **Codex (GPT-5.4)** via `codex:rescue` |
| Bulk mechanical edits, log-scraping, table-building, many independent tasks | **Haiku/Sonnet/Fable** subagents (`model` override) |
| Any run > ~30 min, calibration sweeps, DE/grid searches | **Torch cluster** — smoke-test first, always |
| Model-critical numerics (objective, solver, DP `(1−phi)` threshold, KFE, targets) | **Lead writes, or lead verifies line-by-line against the math** |

The last row is the guardrail. Balanced routing means cheap agents do volume and
independent passes, but anything touching the model spec or the SMM objective is
checked against the equations by the lead. This repo's historical bugs
(`rho_hat` not recomputed, `psi_child` sign, unclamped golden-section kernel)
lived exactly there.

Design choice — balanced routing (not aggressive default-delegate): the lead
never surrenders identification and calibration judgment, because those are the
scarce skill here and the cost of a wrong-but-plausible delegated calibration
conclusion is high.

## Deliverables

### 1. Playbook doc — `docs/workflow/delegation_and_cluster_playbook.md`

The operational reference, loaded when delegation or cluster work starts. Contents:

- The routing standard (table above) with one line of "when to reach for each."
- **Torch SOP** consolidating the current hand-written prompt: SSH alias `torch`,
  the `BatchMode=yes` auth probe, account rule `torch_pr_570_general` (never
  `571`), the working-dir pattern (`$SCRATCH/projects/Fertility_Spring26/code/cluster`
  or the dated snapshot named by `CALIBRATION_STATUS.md`), smoke-test-first,
  explicit run-tag / walltime / array / seed / result-dir, `best.json` +
  `cases.jsonl` health check, and the collect step.
- **Codex handoff**: when to call `codex:rescue` (independent implementation,
  second-opinion diagnosis, stuck), how to scope the handoff, and the rule that
  Codex output touching model-critical numerics is verified by the lead.
- Cross-links to existing discipline: Long-Run Search Safety and Identification
  discipline in `CLAUDE.md`; `verification-before-completion`.

### 2. `CLAUDE.md` + `AGENTS.md` — "Delegation & Cost Discipline" section

One concise section, mirrored into both files in the same edit (repo rule).
States the balanced-routing default, names the primitives, points to the
playbook. No volatile state (no job IDs, tags, losses).

### 3. Global `~/.claude/CLAUDE.md` — cross-project routing block

A distilled, project-agnostic version: lead reasons and reviews; delegate search
to read-only agents, heavy/mechanical coding to Codex or cheap-model subagents,
long runs to a cluster or background; verify delegated output before trusting it.
No Fertility-specific content.

### 4. Lightweight automation

- **SessionStart hook** in `.claude/settings.json`: injects a 3-line reminder —
  the routing default plus "read the playbook before any cluster or Codex work."
  A durable nudge that survives `CLAUDE.md` growing long. Must be robust (never
  block session start).
- **`code/cluster/torch.sh`** wrapper (runs locally, shells over `ssh torch`)
  with subcommands: `status`, `smoke <family> <tag>`, `submit <family> <tag>`,
  `logs [pattern]`, `collect <family> <tag>`. Encodes the account rule, the
  working-dir pattern, and the exact submit/collect commands so the ceremony is
  one consistent call for any agent (especially Fable overnight). Fails loudly
  if the `torch` SSH probe fails (never silently works around auth).

## Out of scope today

- The menu-size / segmentation research plan (already in motion). The new
  machinery becomes its vehicle next session.
- Any model, solver, objective, or calibration change.

## Implementation checklist

1. Explore agent inventories exact current cluster commands (submit scripts, env
   vars, smoke tests, collectors, scratch/snapshot path). — in progress
2. Write `docs/workflow/delegation_and_cluster_playbook.md`.
3. Write `code/cluster/torch.sh` from the verified inventory; `chmod +x`; smoke
   the local arg-parsing (`torch.sh` with no args prints usage) without hitting
   the cluster.
4. Add the mirrored "Delegation & Cost Discipline" section to `CLAUDE.md` and
   `AGENTS.md`.
5. Add / create the cross-project block in `~/.claude/CLAUDE.md`.
6. Add the SessionStart hook to `.claude/settings.json`; validate JSON.
7. Update `code/cluster/README.md` to point at `torch.sh` and the playbook.

## Verification

- `torch.sh` prints usage and rejects unknown subcommands without touching SSH.
- `.claude/settings.json` parses as valid JSON.
- `CLAUDE.md` and `AGENTS.md` "Delegation & Cost Discipline" sections are
  byte-identical (repo mirror rule).
- Playbook's quoted cluster commands match the Explore inventory (which is read
  from the live scripts).
- No model/calibration files touched.
