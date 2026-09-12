# September 12 autonomous cluster campaign

## 13:15 UTC: matched first window, consolidated PDF and final horizon checks

Warm follow-up arm 0 passed its short native root: model fertility **1.974855548** against **1.974875**, gap **-0.000019452**. Housing/PAYGO maxima are 3.22250e-5/3.97512e-7. Complete rows, root receipt, provenance and 17 graphs are collected in `matched_short/`. The -0.0175 arm failed the unchanged mass gate (1.318e-8 versus 1e-8). No accepted four-shock sequence or new policy result exists.

The older -0.045 28-date root finally converged (housing 3.60843e-5, PAYGO 9.19031e-7) but failed terminal-distance checks: household mass gap 12.56% and pension gap 7.78%. It stopped on the historical fitting budget; no accepted state was passed to policy.

Final parallel horizon array **17489644_0–1** checks **28 and 56 dates** at the matching preference decline -0.01414, each with **9,000 seconds, one CPU and 32 GiB**. It reuses the pinned verified short smoke and nearby converged 28-date prices/pensions as numerical guesses; only the 28-date case reuses the dimension-matched Jacobian. Both tasks passed preflight and are running. Root: `/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/night_matched_horizon_20260912`; plan hashes: `c83b770c04c8b7a5820b4e9da7494b2ef0cd614863c65fd8d0f5fb4587dc3386` and `7d2048336848379ca83349946b3c12974c151a54aeecc49fe775aff82b540b28`. These are diagnostic-only checks, not fitted histories or policy launchers.

The consolidated **14-page review** is [e5f_overnight_review.pdf](../../../../pdf/e5f_overnight_review.pdf). It includes every target and parameter, all 17 standard initial graphs, the actual forecast shapes, fiscal receipts and unresolved issues. All 14 pages were rendered and inspected. `report_support/verification.json` records the source fingerprints, recomputed loss, native residuals and reviewed PDF hash. Existing graph legends are preserved, including their original crowding. The PDF is not a certificate of a completed historical model.

Rebuild with `code/model/tools/build_e5f_overnight_review.py` using Python with matplotlib and reportlab. In this session: `PYTHONPATH=/tmp/e5f-pdf-deps MPLCONFIGDIR=/tmp/e5f-mpl /opt/anaconda3/bin/python -B code/model/tools/build_e5f_overnight_review.py`; the isolated `reportlab` package is linked from the Codex bundled Python runtime. No model solve is performed. Final-stage monitoring now checks every 30 minutes, stays quiet while healthy, collects on completion and pauses no later than 16:00 UTC.

## 11:00 UTC: two smaller-shock equilibria and a numerical-start test

Follow-up: **both warm terminal probes passed** (`passed_terminal_root_diagnostic`), and both dependent forecasts 17488254_0–1 are now running. Verified terminal equilibria exist at the newly tested smaller declines. Earlier failures at initial guesses must not be interpreted as general equilibrium nonexistence.

Both frontier jobs completed: first-window fertility is **1.90752953** for -0.0225 and **1.86726183** for -0.0275, against data **1.974875**. Both clear finite six-date housing and PAYGO checks and have 17 standard graphs. [frontier/README.md](frontier/README.md) contains the comparison and points to their complete receipts and graphs. No horizon or historical-fit certificate follows from these short roots. The original long arm is in its third 28-date continuation; the previous continuation's maximum housing/PAYGO residuals were 8.21433e-4/1.77627e-5, still above tolerance.

Failures of smaller shocks at the old initial price/pension guess do **not** establish nonexistence of a terminal equilibrium. Array **17488021_0–1** tests -0.01414/-0.0175 from the nearest verified terminal's numerical starting coordinates, at most eight mappings/30 minutes each. Root: `/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/night_warm_terminal_20260912`. The generated entrypoint's full diff was reviewed: only numerical starts, schema/provenance and explicit source-root location differ; all economic source pins remain unchanged. Generated-driver SHA256: `e820d378a1de532f53f37772ebd4f8b80daba6093be2bfa235b6d8ca7bf7f4b7`. The starting receipt SHA256 is `460a02e584fd8d25b4d85ceb056e7f03d38d8a8de3915b322ed6f3ffc0419841`. Both probes have produced valid native mappings.

Conditional array **17488254_0–1** waits for those probes. A passed terminal is used as a numerical start and independently reproduced before its six-date forecast; a failed terminal skips the dependent forecast. Each forecast has a 90-minute budget and cannot promote a history or launch policies. Root: `/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/night_warm_followup_20260912`; inspect `plan_*.json`, `results/arm_*/summary.json`, `failure.json` and `skipped_*.json`. All 19 existing loop tests still pass locally. Original long-run source and output directories were not changed.

## 09:00 UTC: completed calibration and first native equilibrium

Calibration job 17440306 finished all three bounded batches and 252 search trials. Final loss **158.5411910626344**, a **0.438%** improvement from launch, reproduces twice. Beta is still at its estimated 0.99 cap. [verified_final/README.md](verified_final/README.md) contains all targets, weights, contributions and parameters/bounds; the full 17-graph packet is collected. This candidate does not replace the historical workers' pinned initial condition.

Recovery arm 1 passed the six-date native root: maximum housing residual 1.13348e-4 and pension residual 7.92154e-7, within unchanged gates. Terminal-distance checks fail; its 28-date extension continues with valid mapping progress. The first forecast window gives fertility 1.72574 versus 1.974875 in data. This is an unfitted constant-preference forecast, not the realized historical sequence. [Native forecast PDF](native_smoke/native_forecast.pdf) and [PNG](native_smoke/native_forecast.png) are supplemental; their source rows, fiscal/horizon receipt and plot reproduction script are alongside them.

Recovery arm 0 stopped at a 1.075e-8 relative cohort-mass discrepancy against the unchanged 1e-8 gate (absolute loss 2.75e-10). The old hypothesis of float32 tenure probability normalization has not been verified; no scientific code or gate changed.

Independent diagnostic array **17486113_0–1** tests initial shocks -0.0225/-0.0275. Each gets 90 minutes, one CPU and 16 GiB; both were observed running with all 19 tests passed. These tasks stop after the six-date loop, retain diagnostic/horizon flags and cannot launch a history or policy. This is a bounded early-fertility/terminal-existence check, not another full calibration. Root: `/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/night_surprise_frontier_20260912`; plan SHA256 `9a804837c6ae5b10d4b001833efc0ae96a8fbf441c2c85eabf13c4c9a56b7571`; collector **17486114**. Inspect `results/arm_*/summary.json`, failure receipts and `native_graphs/standard_diagnostics/`. The long recovery arm and its conditional policy logic remain unchanged.

## 06:55 UTC recovery

Arms 1 and 2 of the original array failed at the second native mapping with `KeyError: pension_relative_gap`. The native terminal report returns its global tolerance dictionary by reference; extending that returned dictionary polluted the next call. The adapter now deep-copies the report before adding the pension tail test. A repeated-root regression exercises the shared-dictionary case. All 19 tests pass locally, on Torch and inside the new workers. No equations or thresholds changed.

Recovery array **17473943_0–1** is running with distinct initial preference changes -0.03 and -0.045, seven-hour budgets and the original policy reserve. Original failed outputs and source pins are preserved. New snapshot directory: `/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/night_surprises_recovery_20260912`; plan SHA256 `47e0b7212a9b440144fc7b6bd81760293c852d9407b869c6b70371be2727c7ed`. Read its `results/arm_*/` receipts for live transition progress. Recovery collector **17473944** generates its own `report/` after these workers and calibration end. Native forecast convergence remains pending.

Calibration remains active in batch 2. Completed batch 1 is collected in [verified_round_1/README.md](verified_round_1/README.md), with all target/parameter tables, exact-repeat receipt and standard graphs. The selected loss 158.92995942909278 is an intermediate result; it does not change the historical workers' initial calibration. Weekly remaining at this wakeup: 29%.

The model-task author explicitly authorized overnight calibration, successive unexpected permanent preference changes, and conditional policy experiments, retaining the 0.7202462623815278 first-birth rooms target and its weight. This authorization is separate from the cancelled accidental data-task delegation. No target, utility, empirical estimator, or numerical tolerance is changed.

## Submitted work

| Job | Work | Resources and budget |
|---|---|---|
| 17440306 | Initial calibration refinement; 18 simultaneous coordinate trials followed by joint proposals | 18 CPUs, 96 GiB, 10-hour scheduler cap; up to three bounded search batches, at most 270 new cases |
| 17445263_0–2 | Three independent successive-surprise searches, starting from preference changes -0.005, -0.015, -0.030 | Each 3 CPUs, 48 GiB, 10 hours; history uses one worker, then three policy subprocesses concurrently |
| 17445799 | Graphs and evidence collector, after all above jobs terminate, including failures | 1 CPU, 4 GiB, 10 minutes |

All four numerical tasks were observed RUNNING. Calibration seed reproduced twice. The 18 surprise/controller tests passed locally, on Torch, and in each array task. Native terminal and six-date exact-loop smoke stages still determine whether a long transition is permitted. A small preference decline may fail the retained finite-level demographic endpoint; that rejection is preserved and the next bounded trial is attempted.

At the first follow-up, arm 0 had exhausted its four small-shock terminal trials: all failed finite demographic endpoint existence (last renewal ratio 1.00206). Arms 1 and 2 and calibration continued. The graph collector was smoke-tested successfully on these partial receipts. The failed arm does not terminate its siblings.

The historical workers use the previously twice-reproduced capped calibration, not a changing live search incumbent. Beta remains estimated over [0.94, 0.99]. New initial candidates are reported separately; they cannot silently replace a historical chain's initial distribution or parameters.

## Search and gates

Each surprise search has at most four candidate preference levels for each of four birth windows ending 2011, 2015, 2019, 2023. Every candidate has its own constant-preference terminal solve and forecast. Only the first period is implemented before the next unexpected shock. The inherited distribution and birth queues are retained. Preferences remain constant after the last historical shock.

Each terminal is limited to eight root mappings and 30 minutes. The native six-date loop permits three bounded continuations, then 28- and 56-date forecasts permit up to three each. Each continuation retains the Jacobian and prices/pensions; no tolerance is relaxed. The history reserves two hours for policies within a 35,400-second total worker budget. Maximum counts are ceilings, not promises that every trial fits overnight. Earlier observed 28-date mappings took roughly 12–17 minutes: a four-shock, eight-mapping chain alone can take 6.4–9.1 hours before scalar-search repeats and terminal solves. Independent starts improve the chance of progress; difficult convergence can exhaust the budget.

Acceptance requires housing residual <= 2e-4, PAYGO relative residual <= 1e-6, fresh replay <= 2e-10, all household/population gates, and terminal-distance checks including pension distance <= 1%. The shock-fit threshold is absolute fertility error <= 0.005 for the retained four-year target. This is a declared fitting threshold, not an empirical standard error. Full horizon sensitivity remains unverified. The 0.169 outside-origin entry share remains explicitly diagnostic/outstanding, not a production normalization. The existing household-rate analogue of female TFR is retained and visibly labeled.

After a complete accepted history, each arm automatically attempts baseline 1% property tax without rebate, equal rebate at 1%, and equal rebate at 2%. PAYGO and property-tax rebates have separate budget equations. Failures in one policy subprocess do not stop its siblings. Results remain diagnostic pending the outstanding closure and horizon checks.

## Evidence and monitoring

Calibration root: `/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/night_20260912/search`.

Transition root: `/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/night_surprises_20260912`. Its pinned `plan.json` SHA256 is `470ae06a40f5f4ba39050fd6ff833dbe81836c0cdba3edc5170ec2768c68f37a`; 644 original source entries are pinned, plus separate wrapper and input hashes. No active scientific source was overwritten.

Read `results/arm_*/latest_stage.json`, `heartbeat.json`, trial root receipts, `best_so_far.json`, `realized_fit.json`, and `failure.json`. The collector writes `report/overnight_diagnostics.pdf`, `fertility_model_vs_data.png`, all completed fits and failures, the full historical initial score, and any new selected calibration tables. Stable 17-graph packets remain beside accepted solutions. Blank model coverage is explicitly reported as failure to obtain accepted windows, never shown as a fit.

The existing task heartbeat is active every two hours, quiet while healthy. It checks usage first and pauses agent work at 25% weekly remaining, leaving cluster computation autonomous. The last observed weekly remaining was 32%. No reset redemption is authorized. Morning collection deadline is September 12 at 16:00 UTC; the heartbeat then pauses. Authentication failures are recorded once per wakeup. A 30-minute absence of numerical progress requires investigation.
