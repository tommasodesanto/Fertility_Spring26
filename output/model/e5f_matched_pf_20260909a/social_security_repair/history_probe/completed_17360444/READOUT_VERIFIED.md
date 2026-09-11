# Verified historical conditional replay — job 17360444

The numerical mapping and its exact full-path replay pass. This is not a historical equilibrium, fitted history, calibrated SMM result, or adequate-horizon certificate. The path prescribed prices and pensions; it did not take a root-update step.

Independent verification rehashed all 639 remote source files, seven input pins, 19 collected JSON/CSV files, the saved checkpoint, and all 17 graph files against their receipts. Both mappings agree exactly on 60 policy-array signatures, 18 distribution signatures, prices, pensions, and actual fiscal ledgers. All 12 dated household-audit packets pass. Independently recomputed market/fiscal residuals and CSV-to-ledger comparisons differ by zero. The largest occupied budget excess is 7.11e-15; the largest recorded mass/head identity gap is 1.43e-15. No graph visual review or independent pickle reload was performed by this collection task; the frozen driver performed and certified checkpoint reload.

## Dated residuals

Market residual is (housing demand minus supply)/supply; pension residual is (payroll revenue minus pension outlays)/max(abs(revenue),abs(outlays),1e-12). Positive pension residual here means a revenue surplus. Pensions are four-year period amounts. The last column is an accounting-balanced pension evaluated at this particular population, not a new solved path.

| Year | Asset price | Pension | Housing residual | Pension residual | Implied balanced pension |
|---|---:|---:|---:|---:|---:|
| 2007 | 0.633257 | 2.046361 | -5.830159% | +41.742268% | 3.512601 |
| 2011 | 0.604819 | 1.923736 | -0.068538% | +40.812299% | 3.250230 |
| 2015 | 0.576380 | 1.801111 | +7.667710% | +36.729690% | 2.846692 |
| 2019 | 0.547942 | 1.678486 | +16.949960% | +32.743451% | 2.495646 |
| 2023 | 0.519504 | 1.555860 | +35.543129% | +33.224734% | 2.329995 |
| 2027 | 0.519504 | 1.555860 | +42.640908% | +25.416839% | 2.086075 |

The old stationary pension of 2.046361 is therefore not balanced at the empirically reweighted 2007 population: the conditional ledger implies 3.512601. This is the documented reason to initialize dated pensions from the reweighted dated age/earnings population, subject to the lead’s independent prefix-pension helper validation. From 2023 onward the person/head distribution evolves endogenously, so these implied amounts are starting guesses to be recomputed inside the root, not externally fixed pensions. Property-tax revenue is intentionally unrebated under the maintained historical contract; a nonzero property-tax government surplus is not the PAYGO residual being solved.

## Next bounded actual root

Use the reviewed balanced-history adapter unchanged in a separately pinned actual-root driver. The current probe driver deliberately requires two identical mappings and cannot be relabelled an equilibrium root. Retain the same initial and accepted terminal packets, original supply curve, all demographic arrays, announced preference path, payroll tax 0.179 and zero property-tax rebates.

First solve six dates with eight total evaluations: one fresh initial path, at most six update/trial mappings (rejected trials also count), and one reserved fresh replay. Use explicit price/pension guesses; the current measured accounting-balanced pensions can seed the tail, while the independently reviewed prefix helper should provide the historical guesses. Keep the current price bounds [0.05,5], pension bounds [0.05,10], positive-rent projection, initial block slopes 1.63 and 1, max log step 0.2, damping 1, and exact household/mapping gates unless the lead explicitly pins another reviewed numerical control. Housing tolerance remains 2e-4, fiscal tolerance 1e-6 and fresh reproduction 2e-10. No dense numerical Jacobian is needed for this bounded first round. Report best/final payloads and dated ledgers after every mapping; stop at budget or any invalid household mapping, and never promote an incomplete result.

Measured mappings took 173.855 and 168.472 seconds (mean 171.163); all 24 Bellman calls, checkpoint/reload, and graphs finished in 360.080 seconds. Peak recorded RSS was 2.23 GB. Eight six-date mappings require at most 96 Bellman calls, about 22.82 minutes at the measured mean: use a 27-minute root budget, 30-minute driver watchdog and 32-minute Slurm limit, one CPU, 8 GB with margin; omit manual partition selection. This is a prospective request, not a submission.

After the actual six-date loop has demonstrated valid root updates and honest fresh-final association, the smallest informative longer stage is 28 dates (2007–2115, arriving at a 2119 terminal state), six total mappings including initial and fresh replay. That is 336 Bellman calls, linearly estimated at 79.88 minutes. An explicit 100-minute root budget, 105-minute driver cap and 110-minute Slurm limit leave margin; use one CPU and 16 GB until long-path memory is measured. Eight mappings instead mean 448 calls and 106.50 estimated minutes, needing a larger pinned budget. These long timings are extrapolations, not observed costs. The implementation must preserve all terminal-distance diagnostics and horizon_verified=False; convergence of a finite path is not horizon validation.

The present six-date terminal state has resident-person mass 58.03% and household-head mass 54.13% above the supplied terminal endpoint; normalized household-distribution L1 distance is 0.36336. Thus the current short horizon explicitly fails the retained terminal gates even though terminal price, rent, preference and last pension match by construction.

## Age masses for the lead’s pension-prediction comparison

`observed_age_masses.csv` contains 85 rows (17 model ages × 2007/2011/2015/2019/2023). The 2007 vector is directly in preflight.json → initial_bridge → age_reweight → reweighted_initial_age_mass; its sum is the inherited unit household mass. The original stationary vector is alongside it as stationary_age_mass. For later historical dates, each previous-date transition_path.csv row stores historical_bridge_audit; parse that field with ast.literal_eval and read groups[*].target_mass, mass_after_bridge, and model_age. These are actual pinned empirical bridge receipts, not the age-only illustrative projection in preflight.age_reweight.path.

To obtain the full joint 2007 age/income marginal, load the pinned original initial checkpoint using the frozen source; take stationary_g_pre and apply the exact age factors reweighted_initial_age_mass/stationary_age_mass along axis3. This is precisely e5f_approved_initial_state.py:48–64 and run_e5f_open_population_transition.py:514–534. The resulting NormalizedOldState.initial_state.g_pre is the announcement-state distribution. Summing axes (0,1,2,4,5,6) yields its age marginal. Preserve all other dimensions when computing the wage-type/retirement exposure; do not replace the distribution with age shares alone. The current dated_2023.pkl.gz does not contain the full g2007 array, so no such extraction was falsely inferred from that checkpoint.

The actual 2007 payroll base 3.4237665554 and retiree exposure 0.1744730749 are already in root_receipt.json → final.payload.fiscal_accounts[0], in both dated_reproduction signatures, and in transition_path.csv. These provide the direct numerical prediction check without solving households. Historical later-date bridge function: run_e5f_open_population_transition.py:567–605. All source references are relative to tmp/e5f_matched_pf/code/model/tools/.

Collection and verification only; no source edits, model solves or job launches.
