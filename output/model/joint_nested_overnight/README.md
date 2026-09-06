# Full simultaneous-choice experiment

Experimental source is isolated in `tmp/e5f_joint_nested_full_20260906a`,
branch `codex/joint-nested-full`. Main production code and calibration are
unchanged. Read `CALIBRATION_STATUS.md` for live job/results status.

The fixed-price core passed on Torch (job17070652). Initial full histories
(job17071299) reached all five dates but failed the marginal-probability audit
at 1.0000000000000002 and are not certified. Complementary marginal construction,
parameter metadata and policy diagnostic fixes are included in the corrected
snapshot. The default-off ten-array reproduction in job17074777 passed exactly.

Corrected snapshot:
`/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260906b`.
`frozen_contract.json` mirrors the remote
`output/model/joint_nested_overnight/contract.json`. The latter's SHA is
`7558feaf55ddae7058f481569cda72a8dba5e09a0d094300aa171b66806adb1d`;
scientific bundle `5020a3e77ec8a0f7ee2deb6cd4b642c67dcaa115f26b68df1f14a06e780f9766`.

Job17074777 runs the full-loop smoke: two exact calibration repetitions,
two perturbations of all eleven estimated parameters, and four two-date
policy loops. A long run must first pass all of these. Its controller searches
all eleven coordinates against the original twelve targets and weights.
No rejected proposal receives a reported calibrated loss. The immutable contract
bounds the search, reserves identification and exact-repeat checks, and pins
the final equilibrium-path/PDF driver.

The complete case directories contain `target_fit_long.csv`,
`parameter_table.csv`, `case_receipt.json`, seventeen standard PNGs and a
compressed checkpoint. `best_so_far.json` and `latest_completed_case.json`
provide compact progress. Rejected proposals are recorded separately.
The final local Jacobian reports its own anchor and missing/one-sided columns;
parameter count alone is not identification evidence.

Post-2023 policies retain the closed-national finite-horizon benchmark:
M=0, retention1, inherited four-slot queue, births divided by2.1, supply
elasticity0.63, current prices treated as permanent each date. The four paths
are baseline,20% supply,95% dependent-child LTV and doubled property tax
without rebates. No perfect-foresight, stationary-endpoint or welfare claim.

The cluster generates a provisional review PDF. Before delivering the final
morning PDF under `output/pdf/`, collect the verified winner and policy ledgers,
verify every displayed number, render every page and inspect it. Use the PDF
skill and available document tools for the final user-facing PDF. The existing
`build_e5f_independent_audit_pdf.py` uses Pandoc/XeLaTeX, not ReportLab.
The protected author manuscript remains read-only.

Review evidence: `integration_review.md`, `implementation_review.md`,
`core_review.md`, `controller_review.md`. Delegated changes were reviewed and
corrected by the lead; worker outputs alone are not certification.

## Queued overnight run

Source commit `4b4ba8e` is pushed on `codex/joint-nested-full`. Long job
`17075663` waits on successful completion of smoke `17074777`, with
cancellation on invalid dependency and independent receipt checks before
searching. The existing finite monitor is active every fifteen minutes.
Both anchor histories now pass and reproduce exactly; local `smoke_anchor/`
contains the full target fit, all parameter bounds, receipts and seventeen
standard graphs. The anchor is not a new searched calibration. Seventeen
local reference/plumbing/controller tests pass. The other smoke stages
and the long search still need to complete.

Scheduler status: `ssh torch 'squeue -j 17074777,17075663'`.
Remote search progress will be in `output/model/joint_nested_overnight/search/`.
Remote final policy results will be in `output/model/joint_nested_overnight/equilibrium_path/`.
The final PDF operation marker was run once on September 6 in the launching
turn; do not repeat it for the same logical artifact. No final PDF has yet
been produced for this full-calibration extension.
