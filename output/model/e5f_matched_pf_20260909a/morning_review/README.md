# Initial equilibrium: refreshable working review

`review.pdf` is a four-page advisor review of the revised parenthood-only housing utility and balanced initial pensions. It reports the mapped initial point, not a re-estimated calibration or a certified perfect-foresight/policy result.

- Page 1: verified initial loop and fiscal checks; economic interpretation and pending scope.
- Page 2: all 13 restrictions (one normalization plus twelve proposed scored rows), two validations and both CPS stock projections. Every actual weight and loss is unavailable.
- Page 3: all 17 parameter/restriction rows, including nine structural starting values with supplied bounds and near-bound flags.
- Page 4: measurement limitations, the feasible passive recent-parent observer, next gates and source boundary.

## Inputs and refresh

The builder reads the sibling `initial_fit_readout` CSV tables, summary and preserved observation/smoke receipts. It verifies their source-snapshot hashes, common checkpoint, row counts, null actual weights/losses and moment-gap arithmetic. All printed target/model/gap values and parameter values are checked against extracted PDF text. The report also explains that partial sensitivity cases are not its baseline table source.

`report_status.json` is an explicitly dated operational checkpoint taken from canonical `CALIBRATION_STATUS.md`; update it only from verified status evidence. The builder does not inspect the cluster or infer completion. Rebuilding updates the generated timestamp, not the evidence-as-of timestamp.

Run from the project root:

```bash
/Users/tommasodesanto/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 output/model/e5f_matched_pf_20260909a/morning_review/build_review.py
pdftoppm -r 105 -png output/model/e5f_matched_pf_20260909a/morning_review/review.pdf output/model/e5f_matched_pf_20260909a/morning_review/page
```

The builder deliberately fails if the initial report's null-weight, no-calibration contract changes. Adopting actual weights, changing targets or introducing a new calibrated result requires a corresponding reviewed report revision; simply rerunning this diagnostic builder is insufficient.

Review all four `page-N.png` renders after any substantive refresh. `qa.json` contains output/input hashes and numeric checks; its visual-review field resets to pending each build. Update it only after inspecting the latest rendered PDF. This folder owns all report source, outputs and intermediates. No model solve, source-code change, empirical re-estimation, new figure or policy run is performed.

## Qualifications preserved

The 2.1 normalization is the model completed-fertility restriction, not female period TFR. The recent-parent ACS group is unavailable from a static count mask; a passive realized-birth/empty-home observer is feasible and under development, with its timing/residence mapping still requiring tests. The older all-dependent versus lifetime-childless contrast is not substituted. Family rooms remain an explicit dependent-count proxy; CPS exactly-one fit changes sign with the two age projections.
