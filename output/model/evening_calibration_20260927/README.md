# Six-hour evening calibration, September 27

Author-authorized new Torch run. Revised gated job **18672459** is SEARCHING on
`cs669`, with **24 actual single-thread model workers** verified at17:28EDT.
All six smoke cases passed; lead inspected full14/31 tables and all17plots before
issuing explicit approval. Search and repeats use the same24CPU/128GB allocation.
See `lead_search_approval.json`, `cluster/gated_v2_launch.json` and `cluster/gated_launch.json` for
actual scheduler evidence and preserved rejected/cancelled submissions.

Starting primary score33.820603648867845; relative-identity0.16026598337629389;
block3.3496849785080824. These use different weights and are not comparable loss
scalars. Complete starting fits/parameters: `cluster/smoke_review_v2/smoke_0000_primary/`.
Exact pair/cross-lane comparison: `cluster/smoke_review_v2/full_table_comparison.json`.
Three contact sheets there represent all17plots, byte-identical across six cases.
Early fertility, ownership and wealth remain material misses. Full-grid policy
caveats remain; acceptance permits search, not certification of a final calibration.

The first gated run,18671834, timed out all six smoke objectives after reaching
the fertility target but before completing verification/export. Search never
started. Contractv4 changes only the numerical initial fertility-benefit guess
to0.14281100340255604 and initial bracket step to0.005; it still solves the same
normalization equation with unchanged tolerance and gates. Source/weight hashes
are unchanged. The previous attempts remain in `cluster/smoke_review_v1/`.
Contractv4 and its exact difference receipt are in `cluster/contract_v4/`.

Baseline source equivalence is accepted through the supplemental frozen/current
Torch comparison: all 64 solution arrays and 31 parameters agree exactly. The
original Mac/Torch bitwise-array assertion remains failed and preserved; the
supplement demonstrates platform differences. All 14 moments agree within
1e-13, scientific gates pass, and all 17 baseline plots were inspected. Read
`lead_baseline_review.json` and `cluster/crosshost_v3/` for complete evidence.

## Fixed scope and clock

- Window: September 27 16:12–22:12 EDT; absolute end epoch 1790561520.
- Search cutoff 21:27 EDT, epoch 1790558820. Setup/queue consume the window.
- At most 24 single-thread model workers on Torch; no local model computation.
- At most 384 objectives including smokes and repeats. Revised search cap360
  leaves room for six failed plus six new smokes, six repeats and two numerical
  prerequisite checks: at most380 in total.
- Each objective has a parent-owned 900-second cap and at most 23 stationary
  solves for child-benefit normalization. No restart or extension of old plans.

## Economic contract

Activate the author-adopted DUE existing-owner borrowing rule with separate death
solvency, retaining buyer origination limits. Search the nine existing coordinates
and tenure choice scale; normalize child benefit to completed fertility 2.1.
Tenure scale exploratory bounds [0.001,0.1], proposed in logarithms. Retain all
other bounds, 2% annual real interest, B15 earnings, entry distributions, grids,
housing preferences, rental cap, conception schedule and scientific gates.

All 14 target rows and 31 parameter rows remain visible. Ten scored moments plus
fertility normalization; first births at age 30+, family rooms and older wealth
dispersion become explicit untargeted validation. Counting restrictions does not
establish identification. Recent-parent ownership is a joint anchor for tenure
scale, ownership preference and child housing loading.

Three separately pinned weight lanes: inherited active weights (primary),
identity after fixed relative-gap scaling, and equal block averages of inherited
standardized errors. Exact formulas and caveats are in
`docs/model/e5f_evening_weight_review_20260927.md`. Rescore every lane winner under
the common primary weights; old loss scalars are not directly comparable.

## Required release gates

Remote source hashes and focused tests; authenticated baseline replay; six
exact-loop DUE smoke objectives; full 14-target/31-parameter comparisons and
standard 17-plot inspection. Search requires an explicit pinned lead approval.
Preserve failures; source corrections require new immutable versioned contracts.
No failed gate may be bypassed to fill the time window.

Check actual concurrency, checkpoints and memory; investigate progress stale for
30 minutes. Final per-lane winners require two fresh repeats and standard plot
exports. No claim of an optimum from a finite search. No automatic overnight
continuation: the author returns around 22:00 to decide the next plan.

## Open interpretation issues

Current first-birth housing response is a matched four-year model contrast, not
the exact empirical event-study estimator. Empirical 1.465 provenance is
`sa_rooms_first_birth_v2.do`, with calendar-year windows. The bequest data are
child-directed transfers while the model counts positive estates. These remain
explicit provisional approximations; reviews do not change targets during search.
Full frictionless equilibrium transition remains unsolved and outside this run.

See the three dated review notes in `docs/model/` for weighting, tenure target and
target measurement evidence. No Google ledger edits or outside messages.
