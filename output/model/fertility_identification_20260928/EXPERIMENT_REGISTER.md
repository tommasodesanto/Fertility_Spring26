# Fertility experiment register

This is the authoritative index for experiments in this packet. A calibration
reference is a frozen benchmark; a selected experiment is not adopted unless
the register says so. Report each experiment against its named anchor, and keep
the one-at-a-time economic restriction distinct from any subsequent parameter
refit.

## Frozen reference

**R00 — 2007 stationary reference — block0506, September 28 verified export.**
Primary loss is `19.581310760138322`. Source is
[`resume_v1/selected_export/primary/`](resume_v1/selected_export/primary/);
the authenticated manifest and its exact identity are documented in the
[packet README](README.md). This is the comparison reference, not an
experimental result. Do not overwrite or promote an experimental export into
R00 without an explicit author decision and a recorded before/after identity.

## Completed and authorized experiments

All completed calibration experiments below use target fingerprint
`20e531855075b807886da8e936930dca9090b989b2ce381de39042cf206efdbd`, the same
ten scored rows and three validation rows, with completed fertility normalized
to `2.1` and the renewal gate retained. Their saved selected-case evidence
includes the complete 14-row target fit, 31-row parameter table, and 17
standard diagnostic plots. Consult the linked packet for the full fit and
parameters; scalar losses alone are not a fit comparison.

| ID / status | Question and sole economic change | Anchor and outcome | Run record and disposition |
|---|---|---|---|
| **E01 — completed; not adopted** | Same one-birth model and economic specification as R00; this is a re-estimation within the original ten-parameter bounds, not a new economic treatment. | [`numerical_pair_v1/`](numerical_pair_v1/), retained-start anchor `one_birth_024_gn1_0`; overnight selected result loss `7.826226594410982`. Two final repeats passed. | Two-stream overnight run, array `18766206` task 0; 36-evaluation lane cap, 1 CPU, 24 GB, seven-hour worker allocation, with two repeats reserved. Full evidence: [overnight results](two_stream_overnight_v1/morning_readout_v1/RESULTS.md). Not adopted; R00 remains the reference. |
| **E02 — completed; not adopted** | Adds a conditional extra taste opportunity and inclusive-value term, plus an independent conception draw after a successful first birth within the four-year cell; the event-age proxy remains the same. This is one bundled experimental change. | [`two_births_optimized_v2/run_v1/`](two_births_optimized_v2/run_v1/), initial anchor; overnight selected result `two_birth_024_gn1_0`, loss `7.8420175378092205`. Two final repeats passed. | Same array `18766206`, task 1 and lane budget as E01. Full evidence: [overnight results](two_stream_overnight_v1/morning_readout_v1/RESULTS.md). Not adopted; this remains an experimental alteration to the active sequential birth-choice model. |
| **E03 — completed; not adopted** | Sets the one-birth model's first-birth utility fixed cost to exactly zero; this is a restriction on utility, not a cash expense. | Positive-cost E01 candidate (`first_birth_fixed_cost=0.35270914196085973`, `psi_child=0.12184551693359474`). Mechanical zero-cost center loss `518.3061508777971`; bounded nine-parameter refit loss `37.07645939301989`, with normalized `psi_child=0.10890056006250005`. | Isolated [zero-cost packet](zero_first_birth_cost_v1/README.md); target fingerprint `20e531855075b807886da8e936930dca9090b989b2ce381de39042cf206efdbd`; configuration SHA256 `7035e2c15641a6cf3dd5ddecd3eac1094bbaabe54057a7cf788c72081d73ce19`. Search job `18817312`; saved-output readout `18830124` completed and verified. Budget: at most 16 evaluations, four-hour controller cap, eight stationary solves and 35 minutes per case, with 70 minutes reserved for final repeats. Original readout job `18817376` failed on missing matplotlib; this was repaired in the later readout job. Two refit repeats passed. The separately observed large mechanical loss and locally refitted result do **not** establish infeasibility; no adoption. |
| **E04 — authorized preparation; not launched** | Sets the exact-zero **first-birth** taste scale \(\kappa\); exact zero makes this choice deterministic while conception risk remains. No birth-opportunity or timing change. The fixed cost is held at the E01 value for the initial center, then is free within its original `[0,8]` bound in the nine-parameter refit. | E01 selected candidate: fixed cost `0.35270914196085973`, `psi_child=0.12184551693359474`. This positive-cost E01 anchor is distinct from E03's zero-cost refit. | Staged plan: replay the anchor; evaluate the zero-\(\kappa\) center holding the other nine coordinates fixed, with completed-fertility normalization; if gates permit, bounded nine-parameter refit and repeat checks. Same ten scored/three validation rows and target fingerprint above. Budget, job ID, exact run configuration, and source certification are **not yet established**. No result exists; do not launch from this register alone. |
| **E05 — authorized preparation; not launched** | Sets the exact-zero **later-birth** taste scale \(\kappa\), retaining the first-birth taste scale and all other E01 economics. Exact zero makes the later-birth choice deterministic; conception risk remains. The fixed cost is held at the E01 value for the initial center, then is free within its original `[0,8]` bound in the nine-parameter refit. | Same E01 positive-cost candidate as E04. Run one taste-scale restriction at a time; do not combine E04 and E05. | Same staged control/zero-center/refit/repeat structure and target fingerprint as E04. Budget, job ID, exact run configuration, and source certification are **not yet established**. No result exists; do not launch from this register alone. |

For E04 and E05, `kappa=0` lies outside the original `0.02` lower bound and
must be reported as an explicit experimental restriction. For each initial center the other E01 parameters, including the positive
first-birth fixed cost, are held at the anchor values. If the center passes the
existing gates, the subsequent nine-parameter refit frees the first-birth fixed
cost along with the other eight coordinates, each within its original bound
(the fixed-cost bound is `[0,8]`). A refit changes parameter estimates as
recalibration; it is not part of the one-object economic comparison. Preserve
the complete target fit, every free-parameter estimate and bound status, all 17
standard plots, repeat evidence, target/source fingerprints, budgets, and job
receipts before classifying either experiment as completed.

## Measurement diagnostics (no specification change)

**M01 — age clock and teenage births, completed.** Saved-data jobs 18833274 and
18833422 used zero model solves/imports/checkpoint reads. The matched age-26
comparison does not close E01's early-child-count gap. Births before model entry
are observed in CPS/NCHS, but their age-25 cohort contribution is not identified
by these cross-sectional stock/period-flow measures. No target was replaced.
See [definitions and full age comparisons](measurement_audit_v1/age_tail_v1/README.md).

## E04/E05 numerical preparation gates

The active solver divides by both fertility taste scales, so setting either to
zero without an isolated solver branch is invalid. At zero, the choice value
must be the maximum of the unchanged wait and conception-risk-adjusted attempt
values; choice probabilities select the maximizing action. The proposed exact
tie convention is wait, consistent with the existing deterministic argmax
ordering. This is an explicit numerical tie convention: the positive-scale
softmax limit instead splits exact ties equally. Dead states retain zero choice
probabilities. First-birth `eps_fert` must remain synchronized with its scale.

The isolated binder must disclose the external zero restriction, validate the
other nine native bounds and verify the actual bound and saved parameters.
Validation must cover both action rankings, ties, dead states, each separate
scale restriction, positive-scale replay, probability caches and KFE accounting.
All tests and model smokes belong on Torch. Source implementation, search method
and budgets remain outstanding; this register authorizes no automatic launch.

## Earlier diagnostic profiles

The September 28 early-fertility-weight and fixed continuation-scale profiles
are diagnostic candidates recorded in the [packet README](README.md) and its
saved-review evidence. They use altered weighting or fixed-scale profile
designs, are not R00, and were not adopted. Keep their original common-primary
comparison and profile labels; do not silently substitute them for E01/E02 or
the frozen reference.
