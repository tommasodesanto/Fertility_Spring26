# Completed bounded morning refinement

The author authorized this September 5 same-model search while working on theory.
All 39 evaluations are complete. The best candidate is `round2/task_009`, with
loss 23.15344472757269, 24.0446% below retained production 30.482966707698903.
Two fresh runs exactly reproduce all 12 fit rows, all estimated parameters,
preference normalizations and 253 numeric history entries. Production is not
promoted and no policy response is recalculated. The report contains every
production and candidate target and parameter, including bounds and fixed objects.

## Read first

- Discussion source: `docs/model/e5f_bounded_calibration_refinement_review.md`.
- PDF: `output/pdf/e5f_bounded_calibration_refinement_review.pdf`.
- `completion_receipt.json`: counts, stages, losses, stop reason and CPU allocation.
- `best_so_far.json`: cross-round incumbent with original summary and plan hashes.
- Each stage's `report/all_target_fits.csv` and `report/all_parameters.csv`:
  complete evidence for every evaluated point, including unsuccessful moves.
- The selected directory's `target_fit_long.csv`, `parameter_table.csv`,
  `case_receipt.json`, `dated_state.pkl.gz`, and `standard_diagnostics/`:
  selected solution with its unchanged 17-graph set. The PDF includes that full set.
- `report/coordinate_identification.json`: all 11 smaller-step columns, singular
  values, weakest loadings and one-sided agreement. This did not generate proposals.
- `round1/lead_diagnostic_review.json`: retained nonwinning occupied value flags.
  `round2/task_009/supplemental_spike_exposure.json` shows zero current mass at
  the selected age-30 first-birth spike. Off-support optimality is not certified.

## Immutable stages and budgets

| Stage | Cases | Torch array | Plan SHA-256 |
|---|---:|---|---|
| Smoke |2|17001372|8ed820afb20c6c5851a4b333e2fba81e45cac3069fc12b79d57c0c9957dc666e|
| Small coordinate |23|17002396|04c7bc99b58edeeee74a0887119c40e04e353b93850cb35b5bebdbf584cebbf4|
| Joint patterns |12|17002846|c709cdd0607c68ac5d26a463e1264df510fba5b4f30afbc06bada79e328d8065|
| Exact repeats |2|17003313|40cb7c93e7c53ff55489045838ba0e94fb6b84e22cca5cf3724c151c8b63cffe|

`slurm_final_accounting.psv` confirms 39 completed allocations and 5.8611 allocated
CPU-hours. Six-way maximum concurrency and 30 minutes per case were retained.
The four-hour session budget started 15:27 UTC; all computation completed around
17:05 UTC. No third round or further search is scheduled.

The immutable plans identify every generator argument, input, observer/helper hash,
scientific fingerprint and deadline. New plan generation additionally verifies
its explicitly supplied planner source hash; older plans retain their original
planner lineage. `code_review_resolution.json` documents that correction and
retrieval of the initial checkpoint receipts. No scientific bundle file changed.
The original ridge planner's stricter rejection remains rejected.

## Reproduction and storage

Frozen scientific revision: `ac676c2a2c6b27f92319625f77affde2fd6b4b74`.
Fifty-file scientific bundle:
`630ba20bca6a1b54eb4c46aca904c4a087afb8c808b9c7f4660d5fcd316a970e`.
Target fingerprint:
`3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497`.
Original input fingerprint:
`0afcb82d4735bd15aaa143ea04e3105a5d43df152122d02b983372102f20eef6`.

The production snapshot and all raw case outputs are also retained on Torch at
`/scratch/td2248/projects/Fertility_Spring26_independent_audit_20260905/`.
Local snapshot: `tmp/e5f_overnight_20260905a/frozen_production/code/model/`.
Use the frozen 50-file code for numerical reproduction; the active checkout's
later default-off target extension has a different manifest. A raw case rerun
must use a fresh output plan/directory; the adapter refuses overwriting completed
results. Keep the original generator, center and source/target fingerprints.
The original source JSON is in the immutable plan. Local full checkpoints and
all 819 recorded artifacts have been hash-verified after retrieval.

Regenerate the report and all standard graph pages, with no new solve:

```bash
python3 code/model/tools/build_e5f_bounded_refinement_report.py --run-dir output/model/e5f_morning_refinement_20260905a --source docs/model/e5f_bounded_calibration_refinement_review.md
python3 code/model/tools/build_e5f_independent_audit_pdf.py --source docs/model/e5f_bounded_calibration_refinement_review.md --output output/pdf/e5f_bounded_calibration_refinement_review.pdf --date '5 September 2026' --heading 'FERTILITY / BOUNDED CALIBRATION REFINEMENT' --no-source-index
```

The separate code review is `docs/model/e5f_full_code_correctness_efficiency_review_20260905.md`.
Its incomplete-policy-bundle finding does not establish an affected historical
candidate. Dated elasticity remains 0.63; the inherited initialization still uses
1.75. Neither value or ordering was changed during this search.
