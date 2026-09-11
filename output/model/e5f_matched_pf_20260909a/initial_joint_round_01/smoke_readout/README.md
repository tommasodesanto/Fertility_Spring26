# Completed joint smoke readout

Job **17360699** completed two fresh normalization loops, eight GE solves, in 1,150.19 seconds (scheduler elapsed 19m27s; peak resident memory 2,805,060 KB). This is a diagnostic proposal readout. Actual SMM weights and loss contributions remain unavailable; no total loss or selected candidate is reported. The unchanged source is 7e872053 and the executed smoke contract retains its own identity.

All saved numerical and parameter checks pass locally, including exact early-measurement equality across repetitions. The completed-fertility normalization is 2.100276518, within the unchanged 0.0005 tolerance. The final market residual reconstructed from saved quantities is 2.18313e-7. The initial remote launch check stopped because the driver writes the market quantity report only for the final repetition; this receipt-path failure is preserved in `initial_gate_failure.json`. The corrected remote gate now passes all 634 source pins and 31 artifact hashes, including both checkpoints and the original 17 graphs. A second receipt correction excludes only elapsed solve seconds from exact normalization equality; all other normalization fields, prices, legacy moments and early measurements remain exactly equal. Both corrections preserve the scientific contracts and numerical gates. The reviewed plan SHA256 is `c466736866e19f4b6115272b5e1ebc1d490c985e5fc11b900b68de8ea62a6edd`; the executed smoke contract remains unchanged.

## Complete target fit

The primary projection assumes births are uniformly timed within each four-year cell. Gap is model minus target. Baseline comparisons use exactly the same saved observation definitions.

| Restriction | Target | Baseline | Smoke | Smoke gap | Actual weight | Contribution |
|---|---:|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.10003286 | 2.10027652 | 0.000276517533 | — | — |
| Childless women, ages 40–44 | 0.198278751 | 0.180585019 | 0.205084872 | 0.00680612061 | — | — |
| Exactly one child among mothers, ages 40–44 | 0.213655325 | 0.233236067 | 0.20266911 | -0.0109862148 | — | — |
| Period mean first-birth age | 25.9762639 | 26.5851408 | 25.9787003 | 0.00243640258 | — | — |
| First births at age 30+ | 0.249278013 | 0.265167566 | 0.227355543 | -0.0219224698 | — | — |
| Wealth / annual gross labor earnings | 6.14586139 | 6.09687696 | 5.81349612 | -0.332365276 | — | — |
| Annual bequests / aggregate wealth | 0.0088 | 0.012368035 | 0.0106133475 | 0.00181334755 | — | — |
| Old wealth/income p90 / median, ages 76–84 | 3.51593509 | 3.72975572 | 4.07963706 | 0.563701977 | — | — |
| Mean occupied rooms, capped at 9 | 5.56109738 | 6.24638234 | 5.60932678 | 0.0482294051 | — | — |
| Ownership, heads 30–55 | 0.648334034 | 0.574768359 | 0.631009107 | -0.0173249268 | — | — |
| First-birth room response, −1 to +3 | 0.720246262 | 0.438846663 | 0.691587925 | -0.028658337 | — | — |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347066932 | 0.166215341 | 0.173255385 | -0.173811547 | — | — |
| Recent-parent ownership gap | 0.162895509 | — | — | — | — | — |

Recent-parent ownership remains unavailable. The family-room row retains its dependent-count proxy qualification; the existing age, residence, income, first-birth measurement and uncertainty limitations remain unchanged. Full provenance and reference precision fields are retained in both 13-row target CSVs and `summary.json`. The complete 20-row precision-options table is copied byte for byte and does not activate weights.

## Both CPS timing projections

| Restriction | Projection | Target | Smoke | Gap |
|---|---|---:|---:|---:|
| Childless women, ages 40–44 | uniform_birth_time | 0.198278751 | 0.205084872 | 0.00680612061 |
| Childless women, ages 40–44 | constant_post_cell | 0.198278751 | 0.194598426 | -0.0036803248 |
| Exactly one child among mothers, ages 40–44 | uniform_birth_time | 0.213655325 | 0.20266911 | -0.0109862148 |
| Exactly one child among mothers, ages 40–44 | constant_post_cell | 0.213655325 | 0.180990226 | -0.0326650995 |

## Two validation observations

| Observation | Target | Baseline | Smoke | Smoke gap | Weight / contribution |
|---|---:|---:|---:|---:|---|
| Ownership, heads 25–34 | 0.431158377 | 0.31360335 | 0.350573081 | -0.0805852956 | — / — |
| Old wealth/income median, ages 76–84 | 7.285793 | 6.40231025 | 5.33770884 | -1.94808417 | — / — |

## All 17 parameters and restrictions

| Parameter | Value | Lower | Upper | Near bound | Status |
|---|---:|---:|---:|---|---|
| beta_annual | 0.995156639 | 0.94 | 0.9995 | False | diagnostic candidate; not a certified estimate |
| kappa_fert | 1.48094958 | 0.02 | 50 | False | diagnostic candidate; not a certified estimate |
| kappa_fert_continuation | 1.47461895 | 0.02 | 50 | False | diagnostic candidate; not a certified estimate |
| chi | 1.09752107 | 0.1 | 5 | False | diagnostic candidate; not a certified estimate |
| H0 | 6.61784826 | 0.2 | 80 | False | diagnostic candidate; not a certified estimate |
| theta0 | 0.345934695 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| theta1 | 0.137971713 | 0.02 | 16 | True | diagnostic candidate; not a certified estimate |
| first_birth_fixed_cost | 4.2030072 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| h_P | 1.18188586 | 0.1 | 2.3 | False | diagnostic candidate; not a certified estimate |
| hbar_child_rooms | 0 | — | — | False | zero restriction |
| psi_child | 0.309860746 | — | — | False | normalized to 2.1 |
| payroll_tax | 0.179 | — | — | False | externally fixed |
| pension_period | 2.04636139 | — | — | False | budget derived |
| housing_supply_elasticity | 0.63 | — | — | False | externally fixed |
| tenure_choice_kappa | 0.005 | — | — | False | externally fixed |
| alpha_cons | 0.733 | — | — | False | externally fixed |
| sigma | 2 | — | — | False | externally fixed |

`parameters.csv` is an exact copy of the second repetition receipt. `comparison_parameters.csv` retains the starting value and change for every row. The nine structural coordinates are trial inputs; the initial preference level is normalized, and the remaining rows retain their original restrictions.

## Recovery and subsequent collection

The prepared `collect_joint.py` reads 23 fresh cases plus proposal 21 reused from smoke repetition 02. It calls only `launch_gate.py check-smoke`, verifies all source/seed/contract identities, applies the pinned panel validator, hashes checkpoints and the original 17 graph names per case, and transfers only JSON/CSV observations and receipts. It preserves separate proposal and actual smoke contract hashes. It never submits, retries, ranks candidates, or downloads PNGs. Supply the reviewed plan SHA recorded above; do not use a guessed fingerprint. The collector rejects mixed plan fingerprints during ingestion.

Run its `--remote --plan-sha256 SHA --array-job-id JOB` mode on Torch from outside model source, saving JSON output locally; ingest with `--ingest SAVED_JSON`. All ingestion stays below this readout folder (`array_collected/`, `raw_case_summary.json/csv`, `array_collection_receipt.json`). No array collection or submission was performed by this task. The initial raw light archive and scheduler receipt are retained.
