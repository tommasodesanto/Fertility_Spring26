# Simultaneous tenure-nested choice experiment

Author authorized this isolated experiment September 6, 2026. The production
model, calibration, target fingerprint and numerical gates are unchanged.
This is a fixed-price, fixed-population, fixed-baseline-continuation one-date
diagnostic. It does not evaluate the twelve-moment calibration objective.

## Reproduction

The final exclusive Torch snapshot is
`/scratch/td2248/projects/Fertility_Spring26_joint_nested_experiment_20260906c`.
The corresponding local prepared snapshot is
`tmp/e5f_joint_nested_experiment_20260906c/`.

1. Reconstruct the scientific files using `git archive ac676c2a2c6b27f92319625f77affde2fd6b4b74 code/model` into a new empty snapshot. Copy the independent numerical-audit helper from git revision `83bb064` (SHA ba0f94f52705f43fd0bf7a18ea8c2053f9897674abeb20b030eaf37c57de8623).
2. Copy the three joint-choice source/test files and `code/cluster/submit_e5f_joint_nested_experiment.sh` from this committed implementation. Verify every `driver_hashes` entry in `final_contract.json`, then copy that contract as the snapshot's `contract.json`. The scientific bundle must hash to 630ba20bca6a1b54eb4c46aca904c4a087afb8c808b9c7f4660d5fcd316a970e.
3. Use the previously verified checkpoint `output/model/e5f_overnight_independent_verification_20260905a/numerical_smoke/dated_state.pkl`, SHA bbe10a21a843facaf2bceed56e89281e00992d632cddab640ff0c572d3eb494f. Its existing Torch copy is `/scratch/td2248/projects/Fertility_Spring26_independent_audit_20260905/output/model/independent_numerical_smoke/dated_state.pkl`.
4. Run the pure tests: `python3 -m unittest discover -s code/model/tools -p test_e5f_joint_nested_choice.py -v`.
5. Submit from the frozen `code/cluster`, setting `E5F_JOINT_STAGE=smoke`, `E5F_JOINT_CHECKPOINT` to the absolute checkpoint path, and `E5F_JOINT_CONTRACT_SHA256=b92a8e2a14634213995acd7748ddc65ed2b92604d178d371c040b9f8352876a0`. Use `sbatch --time=00:32:00 submit_e5f_joint_nested_experiment.sh`.
6. Only after checking its completed exact-loop smoke, set `E5F_JOINT_STAGE=panel` and `E5F_JOINT_SMOKE_SHA256` to that run's smoke-summary hash and submit with `--time=00:20:00`. No additional Bellman calls occur in the panel.
7. Collect the snapshot's `output/` into this result directory's `final_run/`. Run `python3 code/model/tools/collect_e5f_joint_nested_experiment.py --root output/model/e5f_joint_nested_experiment_20260906a`. This command verifies the saved states and regenerates the complete comparison figure; each case already contains its supplemental age/wealth/choice plot. The smoke retains two copies of the original seventeen standard reference graphs.

The original source specification is `docs/model/e5f_joint_nested_experiment.md`.
`independent_review.md` is a mathematical review of the first written spec,
not certification of the subsequently written code. `review_resolution.md`
records the lead's adjudication.

## Development receipts

- First smoke 17067221: stopped after 43 seconds on an absolute budget-excess gate. Its output is preserved in `cluster/smoke/`, with the original contract in `contract.json`.
- All seven occupied budget exceptions, including each state's mass and expenditure gap, exactly match the benchmark's previous audit. The revised gate examines positive state-by-state increases in violating mass, without offsets. Absolute violations are still reported.
- Revised smoke 17068080: completed in 73 seconds. Its contract and compact output are `revised_contract.json` and `revised_smoke/`. A final reporting correction then separated explicit births from topcode-adjusted child units in the reference aggregate receipt, matching the case reports. It changed no choice equation.
- Final smoke 17068184: completed, 54.71 seconds in the driver. All twelve reference arrays reproduced exactly twice; all four saved case hashes checked locally before the panel.
- Panel 17068310: completed all 26 cases in 208.53 driver seconds (215 allocated-job seconds), with no additional Bellman solve. All four smoke cases reproduce exactly, all 26 state hashes verify, and the flat-logit limits agree on full arrays. The final receipts are `final_run/panel/summary.json` and `collection_checks.json`.

Both birth measures refer to the four-year model period. Ownership in the
panel covers the entire inherited adult population; it is not the calibration's
ages-30-55 ownership row. The birth response to a first child in the empirical
event study is not measured by these one-date flows.

`reference_target_fits.csv` and `reference_parameters.csv` preserve all twelve
reference fit rows and the full parameter record, including external and
derived objects. They are copied without re-estimation from retained September
4 task_010; the newer unpromoted overnight candidate is not substituted.

## Readout and final checks

Five-page source: `docs/model/e5f_joint_nested_readout.md`; PDF:
`output/pdf/e5f_joint_nested_readout.pdf`. Page 1 is the short advisor readout.
Render with:

```bash
python3 code/model/tools/build_e5f_independent_audit_pdf.py --source docs/model/e5f_joint_nested_readout.md --output output/pdf/e5f_joint_nested_readout.pdf --date '6 September 2026' --heading 'FERTILITY / SIMULTANEOUS CHOICE EXPERIMENT' --no-source-index
```

All 17 standard reference graphs and 26 supplemental case plots were visually
inspected, as were all five final PDF pages. The numerical PDF verification
receipt is `document_verification.json`. No experimental run estimates a new
loss. Maximum joint/control difference is about one birth unit per million
households over the four-year period; scales move levels much more. Full
lifecycle continuation, equilibrium and original-moment identification remain
outstanding. Total elapsed allocation across all four jobs was 395 seconds.
