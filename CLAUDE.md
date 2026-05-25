# Project Agent Instructions

These instructions apply to all work in this repository.

`CLAUDE.md` and `AGENTS.md` are intentionally duplicated so Claude Code and
Codex receive the same project guidance. Any change to one file must be mirrored
in the other in the same edit.

## Core Mandate

This is an academic research project in quantitative macroeconomics, spatial
economics, urban economics, and household finance. The model studies how
housing costs affect fertility, tenure choice, and spatial sorting over the
lifecycle.

Work in this repository should meet professional research standards appropriate
for a top-five or top-field economics journal. That means:

- Treat equations, definitions, calibration targets, and empirical claims as
  objects that must be checked, not prose to be smoothed over.
- When implementing a theoretical object, verify that the code matches the
  mathematical specification exactly.
- When discussing mechanisms or calibration choices, cite the relevant
  economics literature when it matters: Rosen, Roback, Becker, Becker and
  Lewis, Henderson and Ioannides, Sommer, Sullivan, and Verbrugge, Doepke and
  Kindermann, Gyourko, Mayer, and Sinai, Hsieh and Moretti, and related work.
- Use LaTeX notation for mathematical expressions in explanations, notes, and
  documentation.
- Flag weak assumptions, unidentified parameters, stale targets, or
  inconsistencies explicitly. Do not hide uncertainty behind polished language.

## Mandatory Startup

Before proposing changes, interpreting results, or running model code, load
context in this order:

1. `memory/AGENT_MEMORY.md`
2. The latest `memory/daily/YYYY-MM-DD.md`
3. `CALIBRATION_STATUS.md`
4. `SESSION_DIARY.md` only when historical detail is needed

Rules:

- Treat `memory/` notes as first-pass working context.
- Treat transcript bundles under `memory/transcripts/` as evidence, not
  canonical truth.
- If memory files are unavailable, state that and fall back to
  `CALIBRATION_STATUS.md` plus `SESSION_DIARY.md`.
- Do not orient by broad directory scans alone. This repository contains many
  generated logs, archived outputs, and obsolete exploratory files.

## Source Of Truth

Keep root agent instructions stable and procedural. Do not place volatile
calibration state here.

- `CALIBRATION_STATUS.md` is the canonical live calibration and model-status
  note.
- `code/model/` is the active model codebase. The old
  MATLAB active-file index has been archived under
  `calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/status_notes/`.
- `memory/AGENT_MEMORY.md` stores durable cross-session gotchas and recent
  working context.
- `SESSION_DIARY.md` is chronological background, not the first source for the
  current state.
- Older plans, archived scripts, transcript summaries, and remembered loss
  values are background unless `CALIBRATION_STATUS.md` says they are live.

If these sources disagree, pause and report the discrepancy before making
substantive changes.

## Project Organization And File Hygiene

Keep the repository organized by default. New work should fit the existing
layout unless there is a strong reason to add a new top-level category.

- Keep the root directory limited to durable project-control files such as
  `AGENTS.md`, `CLAUDE.md`, `CALIBRATION_STATUS.md`, `SESSION_DIARY.md`,
  `README.md`, workspace files, and stable project notes explicitly intended
  for root.
- Active code belongs under `code/`: model code in `code/model/`, empirical
  scripts in `code/empirical/`, data construction under `code/data/`, and
  cluster launch/collection scripts under `code/cluster/`.
- Active paper-facing LaTeX belongs under `latex/`. Keep this folder to the
  essential current writeups, slides, bibliography/support files, and PDFs.
  Move stale drafts, old build products, and exploratory TeX to
  `latex/archive/`.
- Conceptual notes, prompts, model memos, and non-paper documentation belong
  under `docs/`, with stale material under `docs/archive/`.
- Generated model/report outputs belong under `output/`, ideally in a named
  subfolder with the script or source that produced them. Data-pipeline outputs
  may stay in the relevant `code/data/.../output/` folder when that is the
  established local pattern.
- Calibration histories, obsolete implementations, old cluster outputs, and
  exploratory work that is no longer active should go under
  `calibration_archive/`, not in root or active code folders.
- Avoid creating many small one-off notes or scripts. Prefer one clear driver
  script, one readable output folder, and one README or status note for a
  nontrivial workflow.
- Do not create dated workhorse folders for active code. Use dates for archive
  folders, result folders, or externally meaningful run labels.
- Before adding a new directory, check whether an existing folder already owns
  that responsibility. After adding a durable new artifact, update the nearest
  relevant `README.md` or status file so the next agent can find it.
- Clean obvious transient byproducts from active folders when feasible
  (`.aux`, `.log`, caches, scratch plots), but do not delete or move user work
  or potentially important generated results without explicit instruction.

## Calibration Guidance

Calibration changes frequently. Do not hard-code current targets, losses, best
candidates, parameter regions, job IDs, cluster status, or benchmark
interpretations in these root instruction files.

For any calibration, target, model-fit, benchmark, figure-refresh, or cluster-run
question, read these first:

1. `CALIBRATION_STATUS.md`
2. `code/model/README.md`
3. The active Python calibration files named by `CALIBRATION_STATUS.md`
4. Any live benchmark, diagnostic, or cluster-status file named by
   `CALIBRATION_STATUS.md`

When reporting calibration results, distinguish live benchmark results from
historical or pre-revision results. Never compare losses across target systems,
geographies, room-unit normalizations, or objective definitions unless you have
verified they are comparable.

When reporting any calibration result, include the relevant target values next
to the model moments. Do not report only losses or model moments when the target
system is available.

## Long-Run Search Safety

For calibration sweeps, overnight searches, cluster jobs, or any run expected to
take more than 30 minutes:

1. Estimate the run size before launch: grid dimensions, total solve count, and
   rough wall-clock time using the latest observed solve time.
2. Smoke-test the exact loop structure first, not just a single model solve.
3. Confirm that checkpoints, summaries, and any expected plots are written.
4. Define the time budget, round budget, per-stage budget, and stop criteria.
5. Write progress at least every case or every 5 minutes.
6. Keep two readable artifacts during the run: latest completed-case summary and
   best-so-far summary.
7. If no checkpoint or heartbeat appears for 30 minutes, treat the run as
   unhealthy and investigate.

Do not launch long searches just to "see what happens." The search design must
be stated and recoverable.

## Economic And Numerical Standards

Use the following checks when results look wrong or surprising:

- Value functions should be monotone in wealth unless there is a clearly
  explained numerical or economic reason.
- Choice probabilities must remain in `[0,1]` and aggregate consistently.
- Housing markets and spatial equilibrium conditions must clear to the tolerance
  implied by the active solver.
- Boundary behavior at poorest, richest, youngest, oldest, renter, owner, and
  high-child states must be inspected when policies look irregular.
- Do not accept implausible moments, wrong-signed comparative statics, or broken
  ownership/fertility gradients without diagnosis.

Think about identification when working with calibration or empirical targets:
which variation identifies each parameter, which moments are informative, and
which parameters are likely substitutable.

## Model Solution Diagnostics

Every model solution or equilibrium run should produce an easily accessible
visual diagnostic packet before or alongside tables. A solution is not fully
inspectable if it only writes CSVs, losses, or scalar moments.

- Provide one simple command or function that regenerates the full diagnostic
  packet for the active solution.
- Plot policy functions over the relevant state spaces, including wealth,
  age, fertility/child states, location, tenure or housing product, and any
  active market-choice dimensions.
- Plot prices in all relevant markets and submarkets.
- Plot quantities, demand, supply, market residuals, and capacity/resource use
  in all relevant markets and submarkets.
- Plot the spatial distribution of households, choices, and market quantities
  whenever location is a state or choice.
- Include boundary-state plots for poorest, richest, youngest, oldest, renter,
  owner, parent, nonparent, and high-child states when those margins exist.
- Tables are useful supporting artifacts, but visual diagnostics are the first
  object to inspect because they reveal the full equilibrium shape.

## Terminology

Maintain precise terminology:

- "Fertility" can mean a flow or hazard; "completed fertility" is a stock.
- "Childlessness" is the extensive margin; "number of children" is the intensive
  margin.
- "Housing services", "housing wealth", "house prices", rents, and unit rents
  are distinct objects.
- "Event-study housing responses" are not the same object as fertility parity
  progression.
- Do not reinterpret a one-shot completed-fertility model as a sequential
  second-birth hazard model unless the state space has been changed to support
  that interpretation.

## Stable Project Gotchas

These points have caused repeated errors. Verify against `CALIBRATION_STATUS.md`
if the live model has changed.

- `parity_progression_1to2` should remain diagnostic-only unless it is redefined
  for the active fertility architecture and remeasured consistently.
- The second-birth housing response target disciplines housing demand, not a
  sequential second-birth fertility hazard.
- In the discrete-time model, `phi` is the financed share, so the down-payment
  threshold uses `(1 - phi)`.
- Conditional renter housing policies are not the same as active realized
  housing after the tenure choice.
- `MIGPUMA1` is not residence `PUMA`; do not join it as if it were current
  residence geography.
- Historical pre-MMS or pre-room-scale calibration losses are not directly
  comparable to live broad-core room-scale benchmark losses.
- The repository has substantial generated output. Prefer canonical status files
  and active file indexes over raw search by filename.

## Code And Workflow Standards

- Prefer existing local Python patterns and helper functions over new
  abstractions. Consult archived MATLAB code only for historical comparison or
  reference-parity checks.
- Keep edits scoped to the task. Do not refactor unrelated code or clean
  generated artifacts unless explicitly asked.
- Preserve reproducibility: deterministic seeds, explicit parameter overrides,
  and saved diagnostics for nontrivial runs.
- For empirical code, preserve sample definitions, weights, geography mappings,
  and event-time definitions. If any of these change, document the change.
- For LaTeX, keep notation consistent with the model and avoid claims that are
  stronger than the evidence.
- For paper-facing document edits, make the minimal required change. Do not
  rewrite surrounding prose, opening lines, footnotes, comments, or author notes
  while modifying a separate section. Preserve existing wording unless the user
  explicitly asks for a rewrite; if a note or footnote seems obsolete, flag it
  instead of deleting or rewriting it.
- Use `rg` for code search where available.

## Git And Backup Routine

The current local Git history was restarted on 2026-05-11 and is backed up to
GitHub `main` in `tommasodesanto/Fertility_Spring26`. The old local Git metadata
is preserved on disk under `.git_legacy_2026-05-11_70506af`.

At the start of any substantive session, run `git status -sb`. Before editing a
tracked document or code file, make sure the status is understood. Do not leave
important source files untracked.

At the end of any substantive session:

1. Run the relevant verification.
2. Run `git status -sb` and review the changed paths.
3. Commit coherent source/documentation changes with a short message.
4. Run `git push` to update GitHub `main`.

There is also a local `launchd` backup job under `ops/git-backup/` that runs
daily at 23:40 and pushes dated backup commits to GitHub `main` when there are
non-ignored source/documentation changes. Check `logs/git-backup/` if the
automatic backup appears not to have run.

Do not force-push GitHub `main` again unless the user explicitly asks for a
history rewrite.

## Verification

Before finishing substantive work:

- Run the smallest relevant check that validates the change.
- For documentation-only edits, inspect the rendered or referenced files when
  feasible; no model run is required unless the content changes a reproducible
  claim.
- For MATLAB model changes, run a smoke test or targeted script appropriate to
  the changed code path. If the appropriate command is unclear, inspect
  `CALIBRATION_STATUS.md`, active script headers, and nearby usage before
  choosing.
- For long or cluster-only checks, do not fake verification. State what was not
  run and why.

## Updating These Instructions

These files should stay concise and durable. Move volatile facts to
`CALIBRATION_STATUS.md` or `memory/AGENT_MEMORY.md`; move long procedures to
task-specific docs or scripts.

When an agent repeatedly makes a project-specific mistake, add a short,
verifiable rule here if it should apply to every future session.
