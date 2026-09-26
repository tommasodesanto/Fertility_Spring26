# Collected four-arm utility comparison

**Collection is complete; the planned recalibration was not achieved.** All four searches stopped during their initial populations, and no differential-evolution generation ran. Normal cluster access returned on September 26, permitting collection of the existing share outputs. Each selected original has two independently verified exact repetitions. No model solve was added during collection.

The floor variants increase required housing with children; the share variants instead change the consumption–housing utility shares. Each is paired with linear or diminishing benefits from children at home. These preferences, the reference-rent normalization and the 0.86 curvature remain experimental. The common pension update is author-adopted.

At the retained points, the share variants fit the old mean-room and first-birth housing-response targets better, but fit childlessness, fertility timing, overall ownership, old-age wealth dispersion and recent-parent ownership worse. Their family-room gap also overshoots. Diminishing child benefits worsen childlessness and recent-parent ownership within both families; they improve the family-room gap only in the floor family. The four weighted losses are 285.640, 302.722, 406.642 and 434.646 in the table order below. Unequal search coverage and different selected structural points across families prevent a calibrated or causal specification ranking.

**Input-contract correction for the next calibration:** the author had already chosen AHS mean rooms 5.729. This frozen experiment instead used ACS mean rooms 5.608 and a model observer capped at nine rooms. The results and weights below remain exactly those of the frozen experiment; they do not implement or supersede that earlier author decision. The next contract needs the accepted AHS target, a matching model observer and a declared weighting rule. This repair/collection pass does not choose that rule.

## Every target and retained model moment

All thirteen rows appear below. Each linked original table retains the full-precision target, model moment, gap, weight and contribution; the completed-fertility row is the separate unscored normalization and twelve rows are weighted.

| Moment | Frozen target | Floor, linear | Floor, concave | Shares, linear | Shares, concave |
| --- | ---: | ---: | ---: | ---: | ---: |
| Completed fertility (normalization) | 2.100 | 2.100 | 2.100 | 2.100 | 2.100 |
| Childlessness, ages 40–44 (%) | 19.828 | 18.279 | 17.737 | 17.131 | 16.609 |
| Exactly one child among mothers, ages 40–44 (%) | 21.366 | 23.126 | 23.564 | 23.055 | 23.482 |
| Mean age at first birth (years) | 25.976 | 26.324 | 26.316 | 25.440 | 25.426 |
| First births at age 30 or older (%) | 24.928 | 24.600 | 24.597 | 19.415 | 19.385 |
| Wealth / annual gross earnings | 6.927 | 6.215 | 6.216 | 7.641 | 7.641 |
| Annual bequest flow / wealth (%) | 0.729 | 0.688 | 0.688 | 0.739 | 0.739 |
| Wealth/income p90 / median, ages 76–84 | 3.516 | 2.910 | 2.909 | 2.763 | 2.763 |
| Mean occupied rooms, capped at nine | 5.608 | 6.543 | 6.544 | 5.979 | 5.980 |
| Homeownership, ages 30–55 (%) | 67.626 | 76.181 | 76.161 | 78.792 | 78.771 |
| First-birth housing response (rooms) | 1.465 | 0.785 | 0.785 | 1.561 | 1.552 |
| Room gap: 3+ versus 1–2 children at home | 0.385 | 0.220 | 0.228 | 0.610 | 0.622 |
| Recent-parent ownership gap (percentage points) | 12.761 | 8.781 | 8.482 | 3.774 | 3.596 |

Original full target-fit tables: [Floor, linear](floor_linear/export/target_fit.csv), [Floor, concave](floor_concave/export/target_fit.csv), [Shares, linear](shares_linear/export/target_fit.csv), [Shares, concave](shares_concave/export/target_fit.csv).

Both floor variants retain initial point 0001; both share variants retain initial point 0038. Structural coordinates agree within each linear/concave pair, while the separately normalized child-benefit coefficient differs. Coordinates also differ across the floor/share families, so those comparisons combine specification and selected-point differences.

## Every parameter and restriction

All free, fixed and derived parameters are included. The four original parameter tables retain exact values and original status text. The near-bound flags use the existing one-percent-of-physical-interval screen; they are not evidence of a bound optimum. The two fertility-shock parameters are flagged near the lower bound in every arm.

| Saved parameter | Floor linear | Floor concave | Shares linear | Shares concave | Bounds / restriction |
| --- | ---: | ---: | ---: | ---: | --- |
| `H0` | 8.646 | 8.646 | 6.727 | 6.727 | Free; [0.200, 80.000] |
| `beta_annual` | 0.962 | 0.962 | 0.975 | 0.975 | Free; [0.940, 0.990] |
| `chi` | 1.127 | 1.127 | 1.125 | 1.125 | Free; [0.100, 5.000] |
| `first_birth_fixed_cost` | 0.507 | 0.507 | 0.719 | 0.719 | Free; [0.000, 8.000] |
| `kappa_fert` | 0.209 | 0.209 | 0.222 | 0.222 | Free; [0.020, 50.000]; near lower bound |
| `kappa_fert_continuation` | 0.482 | 0.482 | 0.460 | 0.460 | Free; [0.020, 50.000]; near lower bound |
| `theta0` | 0.088 | 0.088 | 0.162 | 0.162 | Free; [0.000, 8.000] |
| `h_P` | 1.890 | 1.890 | — | — | Free; [0.100, 2.300] |
| `theta1` | 0.008 | 0.008 | 0.008 | 0.008 | externally fixed B15 |
| `psi_child` | 0.134 | 0.144 | 0.104 | 0.110 | normalized to 2.1 |
| `payroll_tax` | 0.080 | 0.080 | 0.080 | 0.080 | derived from adopted pension ratio and baseline demographics |
| `pension_period` | 0.918 | 0.918 | 0.918 | 0.918 | endogenous balanced PAYGO |
| `housing_supply_elasticity` | 0.630 | 0.630 | 0.630 | 0.630 | retained external setting |
| `tenure_choice_kappa` | 0.005 | 0.005 | 0.005 | 0.005 | retained external setting |
| `alpha_cons` | 0.733 | 0.733 | 0.733 | 0.733 | retained external setting |
| `sigma` | 2.000 | 2.000 | 2.000 | 2.000 | retained external setting |
| `selling_cost` | 0.060 | 0.060 | 0.060 | 0.060 | retained external setting |
| `financed_share` | 0.800 | 0.800 | 0.800 | 0.800 | retained external setting |
| `annual_depreciation` | 0.014 | 0.014 | 0.014 | 0.014 | author adopted input |
| `period_depreciation` | 0.055 | 0.055 | 0.055 | 0.055 | four-year compounded |
| `annual_property_tax` | 0.011 | 0.011 | 0.011 | 0.011 | author adopted input |
| `period_property_tax` | 0.042 | 0.042 | 0.042 | 0.042 | four-year linear source convention |
| `income_process` | 15.000 | 15.000 | 15.000 | 15.000 | retained B15 persistent-state count |
| `entrant_conversion_factor` | 0.500 | 0.500 | 0.500 | 0.500 | legacy child-departure diagnostic; inactive in split-birth entry |
| `adult_entry_birth_to_household_conversion` | 0.476 | 0.476 | 0.476 | 0.476 | effective closed stationary birth conversion |
| `child_benefit_exponent` | 1.000 | 0.860 | 1.000 | 0.860 | fixed sensitivity; not estimated |
| `utility_reference_rent` | 0.110 | 0.110 | 0.110 | 0.110 | fixed experimental common normalization |
| `pension_to_gross_worker_earnings` | 0.229 | 0.229 | 0.229 | 0.229 | adopted CPS2007 target |
| `delta_alpha_jump` | — | — | 0.078 | 0.078 | Free; [0.000, 0.250] |
| `delta_alpha` | — | — | 0.039 | 0.039 | Free; [0.000, 0.250] |

Original full parameter tables: [Floor, linear](floor_linear/export/parameters.csv), [Floor, concave](floor_concave/export/parameters.csv), [Shares, linear](shares_linear/export/parameters.csv), [Shares, concave](shares_concave/export/parameters.csv).

Each floor table has 28 rows (eight free coordinates); each share table has 29 rows (nine free coordinates). A missing cell above denotes an inactive specification coordinate, not an estimated zero. Earnings, entry wealth and rank coupling, split entry at ages 16/20, births divided by 2.1 once, dependency departure, transfers, mortality and estate allocation retain the frozen reference. No new mortality/estate/wealth-grid proposal is adopted.

## Complete attempt accounting

The original finite design allowed 652 new objectives. Exactly 114 were attempted: 59 succeeded, 48 were inadmissible, five were recorded as failed, and two as timed out. All attempted objectives now have terminal records. Another 538 planned search slots were unrun. Four initial slots reuse the first verified smoke and are not additional attempts. Raw statuses are preserved; a traceback diagnosis does not rewrite an original failed record.

| Arm | New initial attempts | Initial success | Inadmissible | Raw failed | Raw timed out | Initial reuse | Initial unrun | DE attempted / unrun | Successful smokes | Successful repeats |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Floor, linear | 10 | 2 | 7 | 1 | 0 | 1 | 29 | 0 / 120 | 2 | 2 |
| Floor, concave | 10 | 2 | 7 | 1 | 0 | 1 | 29 | 0 / 120 | 2 | 2 |
| Shares, linear | 39 | 19 | 17 | 1 | 2 | 1 | 0 | 0 / 120 | 2 | 2 |
| Shares, concave | 39 | 20 | 17 | 2 | 0 | 1 | 0 | 0 / 120 | 2 | 2 |

The selected floor points come from only three distinct successful points per arm, versus twenty and twenty-one for the two share arms. All four Slurm jobs exited with code 2 because the intended search was incomplete, despite successful repetitions and exports. The original [bounded failure readout](readout.md) and [deadline postmortem](share_deadline_postmortem.json) preserve the stopping evidence. The remaining concave-share `initial_0022` traceback now also establishes an inner-alarm TimeoutError followed by SystemError wrapping; its [separate causal receipt](shares_concave/export/shares_concave_initial_0022_failure_review_receipt.json) preserves the raw failed classification.

## Scientific and visual evidence

Each selected original was compared separately against both repeats at zero absolute and relative tolerance, including native value/distribution arrays, prices, normalization, moments, full tables and scientific receipt fields. The original exporter performed those comparisons on Torch. Collection authenticated those saved proofs and rehashed the share checkpoints remotely; it did not rerun the solver or unpickle native models locally. Large checkpoints remain remote.

All selected cases pass the inherited normalization, market, fiscal, entry/operator, budget, purchase and occupied-value/probability gates. These are saved-case checks, not grid-convergence evidence or validation of underwater-owner transitions. The unmatched PSID first-birth observer and the SCF child-directed target versus all-positive-estates observer remain explicit limitations.

Each report includes the stable 17-figure diagnostic set and 22 pages. Every actual page and standard figure was reviewed. The saved source figures retain crowded legends on pages 14–17, overlapping income-state ticks on page 18 and consumption-legend/tick overlaps on pages 19/21. No new table or page clipping was found. Tiny gaps use scientific notation in the preserved PDFs.

**Unresolved policy-boundary feature:** in both selected share packets, housing after tenure choice falls from ten to six rooms near the upper wealth boundary for age-30 childless renters; owner-entry probability declines there too. The selected floor packets do not show that housing drop. The selected structural points differ, so this is not causal evidence about the share preference. Saved compact receipts do not establish occupied mass at those specific nodes; occupancy and the economic/numerical explanation remain unknown. No claim of harmlessness or grid robustness is made.

Actual reports and complete QA receipts:

- [Floor, linear report](floor_linear/export/utility_floor_linear_review.pdf) · [visual review](floor_linear/export/visual_review_receipt.json) · [scientific collection](floor_linear/export/collection_receipt.json)
- [Floor, concave report](floor_concave/export/utility_floor_concave_review.pdf) · [visual review](floor_concave/export/visual_review_receipt.json) · [scientific collection](floor_concave/export/collection_receipt.json)
- [Shares, linear report](shares_linear/export/utility_shares_linear_review.pdf) · [visual review](shares_linear/export/visual_review_receipt.json) · [scientific collection](shares_linear/export/share_collection_review_receipt.json)
- [Shares, concave report](shares_concave/export/utility_shares_concave_review.pdf) · [visual review](shares_concave/export/visual_review_receipt.json) · [scientific collection](shares_concave/export/share_collection_review_receipt.json)

The supplemental birth-versus-wait value-gap diagnostic remains [unavailable](supplemental_birth_wait/availability.json): the required branch values and eligibility flags were not saved. No new solve or inversion of saved probabilities was performed.

Original source commit: `d5dbf04d68ff000e1a1b8e66994cdffd31a01580`. Original launch-contract SHA256: `c3fa4d5b7925a54a747c72511538b8182029e1493e4d02e5ef703a30fcecc28a`. Common target/weight/measurement fingerprint: `8fad0155053df30fc0ccf968c733fcd86e916278ad37cb4b1be62462274557d6`. Arm-specific objective fingerprints and all source/input hashes remain in the linked collection receipts. The heartbeat stays paused; collection does not authorize another run.
