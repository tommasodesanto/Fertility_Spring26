# Paired overnight calibration

Submitted once on September 25, 2026 and now complete. Smoke array `18490899`, 40-worker search `18490902`, repetitions `18490903` and export `18490906` all completed with exit 0:0. The final artifacts are in `final_results/`; the launch narrative below remains as the reproducibility record.

## Reviewed results

All 728 planned objectives scored, including four seeds and four selected repetitions. The 17.9% arm completed 2,024 stationary solves and the 8.751% arm 2,217: 4,241 started/completed, zero incomplete. There were no rejected or unrun objectives. Workers finished at 06:12:39 EDT, repetitions at 06:25:16 and exports at 06:25:52, before the amended eight-hour ceiling. Exhausting this finite bank is not evidence of convergence.

| Moment | Target | Tax 17.9% | Tax 8.751% |
| --- | ---: | ---: | ---: |
| Ownership, heads 30–55 (%) | 67.626 | 64.524 | 75.676 |
| Mean occupied rooms | 5.608 | 6.711 | 6.544 |
| First-birth room response | 1.465 | 0.924 | 0.789 |
| Aggregate wealth / annual gross earnings | 6.927 | 5.357 | 6.045 |
| Mean first-birth age | 25.976 | 26.542 | 26.342 |
| Weighted loss, identical objective | — | 291.732 | 280.411 |

The lower-tax selection is closer on wealth and birth timing, but farther from ownership and first-birth housing. The actual balanced pensions are 51.159% and 25.011% of mean gross working earnings. These selected comparisons include recalibration; they are not a tax-only effect. The common structural-seed comparisons, with separate fertility normalization, are saved under `final_results/matched_smokes/`.

The selected points pass the retained numerical gates and match both exact repetitions. Five other high-tax points have tiny budget-excess mass (maximum 1.58e-23). The unchanged frozen rule counts spending gaps above 1e-9 and rejects mass above 2e-10; both selected points have exactly zero excess mass. Do not confuse the spending-gap threshold with the mass threshold or with the strict-zero rule used by an older workflow.

Independent direct reads of all six saved selected/repeated checkpoints on Torch confirm identical native prices, values and distributions, evaluated values/current distributions, pre-choice distributions, fertility utility scales, and byte-identical full target/parameter tables. See `final_results/collection_scientific_repeat_audit.json`. All native moment and normalization equality is additionally verified by the completed exporter; it did not persist separate per-moment hashes. The final selected checkpoint hashes are `de1da882335328f9c0a6ace673c85e10522c8e68c09555f6a7da1717d93c33d7` and `d5ef71bdaf9960273035c722a2428a55f14bab160e0596881c8a981e67b8ead1`; the field `selected_checkpoint_sha256` in inherited receipts names the original seed, not these final cases.

Neither tax rate is adopted by this experiment. A small weighted-loss difference cannot determine the desired pension coverage and benefit level. Both fits still overpredict mean rooms and underpredict the chosen first-birth response; the empirical/model observer differences below remain unresolved. Exact repetition does not prove grid convergence or validate a transition.

- **17.9%:** [reviewed PDF, full tables and 17 figures](final_results/greaney_179/paired_greaney_179_review.pdf), [all target fits](final_results/greaney_179/target_fit.csv), [all parameters and restrictions](final_results/greaney_179/parameters.csv).
- **8.751%:** [reviewed PDF, full tables and 17 figures](final_results/oasi_087510/paired_oasi_087510_review.pdf), [all target fits](final_results/oasi_087510/target_fit.csv), [all parameters and restrictions](final_results/oasi_087510/parameters.csv).
- `final_results/pdf_qa/` records remote rendering of all 42 pages and visual inspection. Some original income-state legends remain crowded; the agreed figure set is preserved. Full-size table inspection found no overlapping columns.

Every free parameter and its actual bound appears in the CSVs. The two taste-shock estimates are inside the inherited 1%-of-full-range lower-bound screen; neither equals the lower bound of 0.02. The bounds are broad (0.02–50), so these flags do not establish binding constraints.

## Common specification and fiscal comparison

The author selected **1.465 rooms**, exactly, for both arms. Rates are 0.179 and 0.08751017424959717 under the same existing pension-only PAYGO rule. The higher rate is a historical DUE-table comparator, not a replication of DUE's broader government budget. Neither rate is adopted as the paper baseline.

Relative to the September 24 commute calibration, the new source enables the adopted birth-based adult-entry accounting and enforces the intended annual beta upper bound of 0.99. It retains B15 persistent earnings with the age profile, inherited heterogeneous entry wealth, floor utility, current-income purchase eligibility, existing transfer rules, the accepted target/cost updates, and existing numerical gates. National housing values use 2005–2006 ACS data for the 2007 reference economy. All 12 positive weights remain fixed. Eight structural coordinates are searched; theta1 is externally fixed and psi is normalized separately to completed fertility 2.1.

The dependency probability remains 2/9 per four-year period. The separate adult-entry queue splits births equally at 16 and 20 years and converts births to households once by 1/2.1, preserving top-bin weights. At constant births this changes the accounting and dated entry clock, not the stationary household decisions. No outside-entry valve or additional parental-death entry was introduced.

The native first-birth observer remains an approximation to the PSID event-study estimator. The model bequest observer includes positive estates of childless households while the SCF target is child-directed. These experimental comparisons do not resolve those measurement differences or establish grid convergence or a validated transition.

## Search and verification

Both arms use the same two structural seeds and the same 360-point bank. Each of 20 workers has nine joint three-coordinate and nine one-coordinate proposals, with both kinds around each seed. There are at most 728 normalized objectives across both arms: four seed evaluations, 720 search evaluations and four selected repetitions. At roughly six native stationary solves per objective this plans for about 4,368 calls, not a promised completion count. The observed 817-second median measures a whole normalized objective.

The run was frozen and submitted with a six-hour clock. The author then explicitly requested about eight hours. Only the live runtime clock was extended: the original deadline is preserved in `deadline_original_6h.json`, and `deadline_amendment_8h_receipt.json` records the authorization and changes. The unchanged frozen lock still describes the original six-hour launch contract. From the first smoke at 00:49:12 EDT, the amended search cutoff is 07:34:12, repetition cutoff 08:34:12, and hard export deadline 08:49:12. Worker/repeat/export stages read this amended clock when starting. Slurm rejected the worker walltime extension; its 6h10 allowance from worker start remains, so some search time may be lost. All job IDs remain unchanged.

A single objective is capped at 3,100 seconds. Queue delay and failures may leave proposals unrun. Unknown failures stop a worker and remain in its records; no retries or gate changes are authorized. Initial progress confirmed two complete stationary solves in the high-tax smoke and one in the low-tax smoke, with each proceeding to its next normalization solve; no scored new objective was yet available at that check.

Eight focused unit/caller tests passed on Torch. Native zero-solve preflight bound both seeds at both rates and checked the generated runtime sources. A three-page report-only fixture verified the table layout for all 13 target rows and 25 parameter/restriction rows. Independent and lead source reviews are in the [maturation review packet](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/demographic_accounting_review/maturation_implementation_20260925/README.md).

The final exporter compared the original selected checkpoint against both exact repetitions for native prices, V, g, policy values, current/pre distributions, all native moments, psi, normalization, source/target contracts, tables and loss. Only stationary_solve_seconds is excluded from normalization comparison. Both numerical comparison tolerances are zero; all four repetitions passed. Final all-page visual QA is complete.

## Files and responsibility

- [Complete target definitions and weights](contract_review/target_decision_table.csv).
- [Targets, bounds, external restrictions and provenance](contract_review/proposed_common_contract.json).
- `runner/`: frozen runner, launcher, input builders, report exporter and compact validation/launch receipts.
- Remote source and results: `/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/nightpair_20260925_v1/`; results in `results/run_001/`.
- Ready lock SHA256: `6443195fa3f7de0dce5cc8a4c05e2709b99d586421c96a35dba3c07e93e061a1`.
- Objective SHA256: `cce32a5d8208603e9237607da8253995a5f5517ef39fd5a1ee75bee90718ed58`.
- Source inventory SHA256: `237904131d159f775c7ae89d1bbf1e8d1dd70c79ad012c80a36d948658d6f9c6`.

The existing cluster integration agent owns initial smoke monitoring and compact collection. The lead owns final numerical/economic review, canonical status and delivery. All heavy execution, hashing and PDF/plot rendering stay on Torch. No giant checkpoints are to be downloaded. Deliver only the paired report: the ancestor PDF and receipt are quarantined under `INTERNAL_NOT_FOR_DELIVERY` because their old front-page labels do not describe this experiment.
