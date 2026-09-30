# Codex worker task

## Goal
Prepare (DO NOT SUBMIT) a tightly bounded Torch native compiled validation of the
isolated renter-taper removal. Lead will review then submit. Author says fix it
now and validate; positive new unsecured credit is NOT selected.

## Scope and ownership
Own only credit_no_taper_v1/native_validation_v1/ under this output folder.
Read existing patch_driver.py, test_renter_no_taper.py, verification_v1/, README.md,
and the parent run_fixed_price.py and credit_v1/run_credit.py for reusable runtime
authentication and checks. No edits to any existing sources or receipts.
Use task-template scope. Full mandatory startup, narrow sections of huge status.
Read docs/workflow/delegation_and_cluster_playbook.md. Cheaper worker task limit
20 minutes. Return best preparation if incomplete; no restart or broad audit.

## Mathematical specification
Reference “2007 stationary reference — block0506, September 28 verified export”.
Use existing Torch isolated overlay_v1 parameters.py, SHA
9c6def300f76b2d5ac55c392e8a595fca881b1d78b15c1cba96016a30a3b83b9.
Native solver and kernels remain hash-pinned, unchanged. Fresh-process overlay
registration MUST precede runtime package imports. Enable exactly
renter_no_taper_estate_bound=True and rebuild debt caps on an authenticated
deepcopy of reference P. lambda_d=0. For death-impossible choice j:
b' >= min(b,0); possible death or terminal: b' >= 0. Owners/buyers unchanged.
No preferences, psi, income, entry, grid, fiscal, supply, target, timing changes.
Retain original 160-node grid. Do not turn on natural credit or renormalize births.

## Required deliverable
One driver and launcher + plan and short README. Reuse parent helper functions
where sound; do not rewrite a solver or copy large repository/checkpoints.
Torch read-only source paths may be inspected. DO NOT submit or run numerical
work; lead owns launch. No heavy Mac imports/tests/rendering.
Controller exact order: baseline flag-off control first, taper-removal second,
one fresh process each, one lifecycle solve per case maximum. Authenticate all
reference checkpoint/source/target/parameter pins. Control must reproduce all
113 reference arrays and 14 fit/31 parameter rows plus standard 17 diagnostics.
Then changed case validates actual compiled household paths, budget, mass,
probabilities, occupied monotonicity, transaction accounting, nonnegative estates,
and operative renter lower bounds. Collect occupied-state incidence by age:
current negative renter assets, policy saving negative, old/new floor binding and
violations; use ACTUAL branch choices/distribution, not full-grid averages.
No assumption that all global loc_probs equal one; unused arrays are not choices.
Keep full 14 fit/31 parameter tables and 17 standard diagnostic plots.
Fixed prices/rents, fixed psi and fiscal rules. These are validation/conditional
responses, NOT cleared equilibrium, demographic stationary endpoint or transition.
Report actual fiscal residuals: do not change pension to clear them or silently
relax an existing gate; if a maintained gate fails preserve failure and stop.
No large checkpoint download; if runtime state useful save it only on Torch with
atomic temporary-file rename. Avoid saving when unnecessary for validation.

## Budget and provenance
Exactly 2 lifecycle solves, 360 seconds per entire case (including reporting),
1200 seconds total from launcher entry; one CPU, 24 GiB, 20-minute Slurm cap.
No retry. Fixed start/deadline in launch.json. Every phase writes progress, latest
completed and best-so-far. Existing old search/recovery budgets are closed; this
is a NEW validation-only budget, not a restart of GE or elasticity.
Pin helpers, overlay, all files and plan immutably before lead submission; no
mutations beneath jobs. New remote root /scratch/td2248/projects/fixed_reference_credit_no_taper_validation_20260929.
Original source root /scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
with original container mount /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26.
Account torch_pr_570_general; use existing anaconda/container pattern.

## Verification / stop
Prepare zero-solve exact-loop mock tests (success, control failure, timeout,
duplicate-output refusal and changed pin) for Torch, gated before numerical cases.
Return concise diff/files, exact command, limitations and any scientific ambiguity.
Stop if rule needs economic assumptions, targets/identity conflict, missing source,
unsafe branch or budget exceeded. No git commits/push; no other worker dispatch.
