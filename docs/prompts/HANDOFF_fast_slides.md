# Fast slide handoff — September 13, 2026

You are the execution-focused Terra-medium chat for Tommaso's September 14 mock
slides. Finish the two figure updates below first. Be extremely economical:
credits are exhausted and usage is now paid. Work directly on small edits;
use a cheap bounded worker only when it saves work. Do not supervise agents in
long polling loops, reopen the calibration, scan chat archives, or start new
scientific experiments. Report a concrete blocker promptly rather than spending
ten minutes silently investigating. Do not resume the paused automatic monitor.

Project: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.
Follow `AGENTS.md` startup, but keep it bounded. Memory's leading September 7
entry is stale and the September 12 daily note is a placeholder; the leading
`CALIBRATION_STATUS.md` entries and this handoff contain the relevant state.
Read the `fertility-paper-slides` skill and presentation style guide. The subtree
`latex/JMP_DS_draft/` is strictly read-only. The worktree has many unrelated
edits: preserve them and stage only your work.

## First deliverable

Refresh exactly these frames in
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/september_14_presentation.tex`:

1. **Lifecycle Fit in 2023**: retain the three panels—ownership, rooms capped
   at nine, and households with children at home—against the existing ACS data.
2. **Intergenerational Allocation: Model vs Data**: retain the six bars for
   ages 22–39, 40–59 and 60–85, each split by children at home. The denominator
   is all large owner-occupied homes (six or more rooms), not all households.

Replace the old patch model curves with the **2023 cross-section of the saved
one-permanent-shock experiment, iteration 3**. Preserve the graph formats and
empirical definitions. It is unconverged, so keep the short numerical-status
qualification already implemented in the new plotting route. Do not substitute
a stationary economy or the old historical fit.

## Recovery already running — do not duplicate it

At the handoff, Torch job **17707689** was RUNNING, 10/22 backward steps done.
Remote batch:
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/permanent_2023_profile_20260913_retry2`.

Check its `output/backward_progress.json`, `output/latest_date.json`,
`output/verification.json`, `output/failure.json`, and `logs/recovery_17707689.err`.
Use `ssh torch`. The job has its own manifest/deadline and runs independently.
Earlier base/retry attempts failed on a list-versus-array interface; use retry2.
The corrected source is already running—do not overwrite frozen cluster files.

This is one fixed-price recovery, not a new root or calibration. Because there
are 17 finite household age cells, 22 backward dates recover the continuation
needed for 2007–2023 exactly. Five dates are then carried forward. **All five
aggregate rows must reproduce the saved 100-date iteration to 2e-10.** A failed
reproduction means the new figures cannot be used.

It saves `model_2023.json`, `native_2023_snapshot.pkl.gz`,
`continuation_2027.pkl.gz`, `rows.json`, `fertility.json`, and `verification.json`.
Collect into:
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_original_queue_20260913a/terminal_restart_v1/fertility_replay_iter3/profile_2023`.

Worker source (backed up as work in progress; end-to-end validation remains pending):
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/cluster/recover_e5f_permanent_2023_profile.py`.
The profile excludes childless readiness states from “children at home.” Check
that housing and household totals reproduce the native 2023 aggregate row.

## Plotting is implemented; finish verification and deck insertion

The work-in-progress addition to
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/build_e5f_patch_readout.py`
accepts `--permanent-profile PATH/model_2023.json --figures-dir PATH/figures`.
Use `/Users/tommasodesanto/miniconda3/bin/python` locally; default Python lacks
NumPy. The new route requires the adjacent passing verification and generates
both PNG/PDF figures plus source/value checks. Existing routes are unchanged.

Verify the source hash, row-reproduction gap, age aggregation, six shares summing
to 100, and visual layout. As an additional check, directly measured completed
fertility in the model's age-42 cell should match the previously reconstructed
2023 value **1.8230999869371063**; investigate any material discrepancy rather
than silently changing the earlier graph.

Then replace only the two `includegraphics` paths (currently around lines
604 and 609). Compile twice, inspect these two rendered frames, and deliver
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/pdf/september_14_presentation.pdf`
plus the usual adjacent `latex/september_14_presentation.pdf`. Keep source
edits minimal. Record the provenance in the local profile folder README and
`CALIBRATION_STATUS.md`; commit coherent changes and push without unrelated edits.

## Transition graphs already ready

Folder:
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_original_queue_20260913a/terminal_restart_v1/fertility_replay_iter3/output`.

- `macro_transition_mock.png/pdf`: four panels through 2103.
- `macro_transition_full_horizon.png/pdf`: the entire saved path through 2403.
- `macro_transition_mock.csv` and `macro_transition_verification.json`: values
  and definitions, including reconstructed completed fertility.
- `fertility_data_model.png/pdf`: the simpler orange model curve versus data.

These are all **one permanent preference decline of 38.1%, not a fitted
four-shock history**. The initial stationary fertility is about 2.10. The saved
path is unconverged; the terminal steady state has independently passed fresh
reproduction and stationarity checks. Housing demand and supply are shown
separately because markets are not yet cleared. Dates retain the existing
start-of-window convention: 2007 data denotes 2008–2011; 2019 denotes 2020–2023.
Do not label household mass as resident-person population.

The announced four-shock alternative has not produced a usable path: its smoke
passed the stationary checks but failed a dated household budget check under
the test price guess. Do not wait for it or launch a replacement for this task.

## Technical reference chat

**Review quantitative model**, task ID
`01a06dd4-9a45-7ff1-bab9-cd22c98c2a29`, host `local`, is the predecessor chat.
Tommaso authorizes contacting it for a precise question about technical history
or a scientific decision. Prefer a bounded task read; send a message only if
needed, requesting a brief verdict and exact source. Do not override its model
or create a duplicate quantitative audit. This handoff should suffice for the
two slide updates. The prior subagent is stopped and is no longer editing.
