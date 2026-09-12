# September 12 autonomous cluster campaign

The model-task author explicitly authorized overnight calibration, successive unexpected permanent preference changes, and conditional policy experiments, retaining the 0.7202462623815278 first-birth rooms target and its weight. This authorization is separate from the cancelled accidental data-task delegation. No target, utility, empirical estimator, or numerical tolerance is changed.

## Submitted work

| Job | Work | Resources and budget |
|---|---|---|
| 17440306 | Initial calibration refinement; 18 simultaneous coordinate trials followed by joint proposals | 18 CPUs, 96 GiB, 10-hour scheduler cap; up to three bounded search batches, at most 270 new cases |
| 17445263_0–2 | Three independent successive-surprise searches, starting from preference changes -0.005, -0.015, -0.030 | Each 3 CPUs, 48 GiB, 10 hours; history uses one worker, then three policy subprocesses concurrently |
| 17445799 | Graphs and evidence collector, after all above jobs terminate, including failures | 1 CPU, 4 GiB, 10 minutes |

All four numerical tasks were observed RUNNING. Calibration seed reproduced twice. The 18 surprise/controller tests passed locally, on Torch, and in each array task. Native terminal and six-date exact-loop smoke stages still determine whether a long transition is permitted. A small preference decline may fail the retained finite-level demographic endpoint; that rejection is preserved and the next bounded trial is attempted.

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
