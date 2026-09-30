# Single-market phase 1: full equilibrium equivalence passed

Torch job **18880497** passed against the existing 160×15 control: six lifecycle evaluations each, eight exact closure receipts, all 87 arrays in each of two saved bundles, every 14/31 table row and 17 actual PNG hashes at each final. Compact authenticated receipts are in [collected](collected/comparison.json). Both runs use the same experimental D=0.14. This verifies implementation equivalence for that specification, not complete location-axis removal or a measured speedup. The preparation and launch history below remains for reproducibility.

# Single-market phase-1 full-closure verification — submitted, not verified

This is one verification arm for the committed single-market specialization
(`1e24652d`). It compares against the **160 × 15 control** of Torch job 18879780,
which already runs the reviewed indexed refactor on `cs713`; no duplicate
control solve is prepared. Phase 1 keeps singleton location axes. This packet
cannot certify complete axis removal or a measured speedup.

## Economic and numerical contract

The comparison changes source implementation only. Both arms use the same
160-node wealth grid, 15-state Markov income process, retained entrant
conditional distribution, preferences, fixed child-benefit parameter, targets,
fiscal objects and housing supply. Both use the already authorized experimental
constant unsecured-credit limit **D = 0.14**, no renter taper, full raw repayment
when selling into renting, and nonnegative estates. Those credit rules differ
from the frozen September 28 reference; they are not a new adopted reference.
No recalibration, normalization of fertility, entrant change or transition is
included. Price clears actual birth renewal; population scale clears fixed
absolute physical housing supply. The default normalized-population lab GE
command is not used.

## Reviewed source and execution

`prepare.py` copies the four prior reviewed indexed-driver modules without
changing their bytes. It materializes the committed cleanup engine as
`source/source/small_credit_lab`, including its relative-import `contract.py`.
The unused joint, two-shock and fertility-nest modules stay absent. Runtime
imports are self-contained; frozen code is loaded only for the authenticated
observer contract and unchanged 14 target rows, 31 parameter rows and 17 plots.
`preparation.json` pins every copied origin, while `source/source.sha256` pins
all staged runtime files. `source.tar.gz` and `source_archive.sha256` are ready
for transfer to the new, exclusive remote directory:
`/scratch/td2248/projects/single_market_verification_v1`.

`run_arm.py` supplies the same fixed-D14 seed call as the grid control and calls
the unchanged full GE root loop. The global deadline is **1200 seconds from
launcher entry**, with at most **six lifecycle evaluations**, **300 seconds per
case**, one CPU and 24 GiB. The inherited root search reserves 700 seconds for
reporting and the selected exact repeat. A timeout/failure ends the attempt;
there is no retry, cap widening or economic fallback.

`launch_torch.sh` requests `cs713`, authenticates staged files and bundle, runs
an authenticated zero-solve smoke first, then runs production with a fresh
Numba cache. It uses the existing anaconda3/2025.06 Python and container, with
all compute thread limits one. The frozen repository and input bundle are
bound read-only. Outputs are never overwritten. The baseline path is
`/scratch/td2248/projects/grid_resolution_120x9_v1/results/full/control_160x15`.
If its completion receipt is absent, this launcher preserves the cleanup
results and writes `comparison_pending.json` instead of waiting or solving.

## Preparation verification

`preflight.json` records the completed local exact GE-loop smoke with mocked
model calls and **zero lifecycle evaluations**. It verifies seed/root/repeat
controller shape using a synthetic seed, six-case cap rejection and expired
deadline rejection. Local observer data are explicitly mocked; this is not a
numerical or frozen-observer certificate. The remote smoke authenticates the
frozen observer before any production solve.

All copied model bytes match commit `1e24652d`; all original source pins remain
unchanged; the deletion inventory is checked; launcher shell syntax passes.
The source package's 19 component tests were already independently reviewed
and passed before this packet. No SSH, transfer, submission, lifecycle or GE
calculation occurred during preparation. The lead subsequently deployed the
pinned archive and submitted replacement Torch **18880497** directly, with no `afterok`
dependency. `launch.json` records the submitted, unverified state and source
archive hash. The initial job **18880464** was cancelled while pending, with
zero numerical work, after discovery of an extra directory level during
archive extraction. The lead moved `initial_extraction/source` to the proper
remote `root/source`; all 25 staged source hashes passed. The replacement uses
the same archive, launcher and budget, with no source change or numerical
retry. The 18879780 control arm passed, but that job later failed in
the separate 120×9 proposal at its entry relocation guard; the passed control
arm is the comparison baseline. No cleanup verification result is available yet.

The lead has deployed the source archive to the exclusive remote directory
and submitted the reviewed launcher. Results will be collected here after the
job completes. No command in this packet submits itself.

## Comparison contract

`compare.py` requires passed full-closure receipts and six lifecycle cases in
both arms, equal input identity and D, exact closure dictionaries, every finite
native array in every `solution_arrays.npz`, all 14/31 rows, and all 17 actual
PNG hashes at both selected and selected-repeat finals. Additional grid-only
`common_support_policies.npz` files are excluded. It does not alter output or
use common-support interpolation to establish source equivalence. Run it on
Torch against completed outputs; copying native array archives locally is not
needed for this check. Runtime node/Python/NumPy/thread information is retained
for auditing; timing differences alone do not establish a matched speedup.
