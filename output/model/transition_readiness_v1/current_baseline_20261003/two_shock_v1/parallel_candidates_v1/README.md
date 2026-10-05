# Parallel first-shock candidate solves

Author-directed priority: obtain the two-shock 2023 numbers using at least 12 concurrent Torch nodes. This fresh wave has six fixed first-shock preference values, each at both 24 and 32 four-year periods: 12 separate jobs. It preserves the healthy fixed-preference diagnostic 19193055. No earlier job is restarted.

## Contract and verified facts

Reference: one-birth Estate-A case `20261003T212605039706Z_a739edc3`, chain 13, baseline preference 0.17892072066041628. Post-interest transactions, soft financing, net estates, fixed-H0 elastic housing supply, actual inherited population, both birth queues and fixed payroll tax with endogenous pension remain unchanged. Economic changes relative to the authorized two-shock contract: **none**. No concurrent calibration changes are adopted.

Each job freshly reconstructs the reference, measures five 12-date/date-5 native seed maps, and solves its own stationary endpoint. It evaluates one fixed preference at one horizon using the original joint price/pension solver and its 12-map allowance including fresh replay. No old reference, endpoint, measured seed, root or optimizer output is reused. The only numerical path change from frozen v5 is the cold price guess: with A = R_gross + delta + tau_H, fixed initial rent r = 0.13689881249028354, and freshly solved endpoint qT, q[t] = r/(A-1) + (qT-r/(A-1))*A**(-(H-t)). Pension starts at the newly solved endpoint. Subsequent solver updates are unchanged.

The six preference values are 0.13, 0.14, 0.145, 0.14736308634876963, 0.15 and 0.16. The original absolute bounds remain [0.001789207206604163, 0.35784144132083257]. The informative first-stage target is 2012–2015 fertility 1.861, measured at local fertility index 1. No preference estimate is accepted yet.

| Birth window | Target | Weight | Model/gap/loss |
|---|---:|---:|---|
| 2008–2011 | 1.974875 | 0 | Pending |
| 2012–2015 | 1.861000 | 1 | Pending |
| 2016–2019 | 1.755375 | 0 | Pending second shock |
| 2020–2023 | 1.645750 | 1 | Pending second shock |

Each job requests one CPU and 24 GiB on a distinct `cs` host, for at most 10,800 external seconds and 10,700 internal seconds. Original inner caps remain: candidate 7,200 seconds, path 6,000, endpoint/seed/mapping 1,800 each, 48 endpoint evaluations and 12 path evaluations. Native calls are capped at 20,000 per job; actual expected counts are far lower. Based on observed 24-date map time 416 seconds plus fresh preparation, approximate job time is 6,500–9,500 seconds, with uncertainty and all caps governing. The paired 32-date paths may reach the unchanged path cap; no extension is implicit.

The five focused fake-native-loop tests execute the copied adapter and cover original gates, full checkpoints, plot hashes, failure receipts, queue/source rejection and pin tampering. Lead AST checks verify exact original endpoint, mapping, Jacobian extension, solver arguments and terminal-check calls; only the cold initializer and reporting differ. These cheap tests are execution checks, not fresh native convergence evidence. The launcher performs full immutable-v5 authentication on host and inside the container, enforces independent persistent claims, and saves actual terminal receipts.

## Collection and continuation

`submission.json` preserves the first scheduling requests; these were all cancelled before allocation (zero elapsed time and no claims/outputs) because their initially selected hosts were reserved. `pending_node_correction.json` records Torch rejecting the in-place host update. `unstarted_cancellation_verified.json` proves all twelve ended unstarted. **`replacement_submission.json` is the authoritative active submission record**, with each replacement attempt saved before sbatch. No numerical attempt was repeated. `deployment_receipt.json` pins transferred files; `lead_review.json` and `config_check.json` record review and the zero-native-call input check. Never duplicate a task with any known or unknown submission outcome. Each remote claim is `jobs/TASK_ID.claim`; outputs are `results/TASK_ID/run/`.

Monitor native counts in `progress.json`, five-map seed evidence, and `candidate/map_*/native_record.json`, `candidate/latest_completed_full.json`, `candidate/best_so_far_full.json` and `candidate/root_progress_full.json`. Preserve exact failure tracebacks and completed native artifacts. `result.json` distinguishes a root-gate pass from acceptance: a single horizon always has `accepted_candidate=false` and `pending_horizon_comparison=true`.

Require both 24/32 roots, fresh replay, source/target fingerprints, the original historical-horizon comparison and fertility gap at most 0.005 before any first-stage selection. Target identification and scalar-fit/final-reproduction evidence remain required. Exported 2015 packets retain the exact initial and inherited states, both queues, V0/V1/V2, terminal value, full prices/pensions/preferences and genuine endpoint packet. They support the unchanged two-period prefix replay. Only after that accepted own-vintage prefix may the unexpected 2015 second shock be fitted from its actual state and the [2:14] forecast slice. The stage-one 2023 continuation exported here is **not** the requested two-shock 2023 state.

The 17 standard plot names are preserved when original root gates pass. Diagnostic 24/32 horizons do not certify 104/128 production paths; estate closure and terminal/projection limitations remain explicit. The frozen September 14 presentation used 100 periods after four announced-shock periods (104 total), and its retained record labels that curve provisional and unconverged. It has a different model and information contract.

## Active scheduling records

**Verified October 5 at 01:18 New York:** five fresh 24-period candidates (preferences 0.14, 0.145, 0.14736308634876963, 0.15 and 0.16) completed with exit 0 and all original root/accounting/fresh-replay gates passing. Four used 455 native calls and ten root maps; the 0.16 case used 413 calls and nine maps. All 1,176 dated audits pass, projected mass is zero, and fresh replay discrepancy is zero. Terminal diagnostics remain false and nongating. The closest fixed preference produces 2012–2015 fertility **1.861092771 versus target 1.861**, but no preference estimate is accepted. [Full tables, bounds and lead review](completion_review.json).

The 32-period candidate **c02_h32 / 19195638 failed at 01:17:41**, exit 1:0, after 527 actual native calls: its original 6,000-second path deadline expired during the ninth mapping. Eight completed mappings account for 501 cumulative calls including preparation; 26 additional entries belong to the unfinished map. Last maximum housing/fiscal residuals were **6.1279227e-5 / 5.1678709e-5**: housing passed, fiscal did not, and no final replay occurred. The underlying TimeoutError was wrapped by Numba; this was neither OOM nor an external scheduler timeout. No job was extended or restarted. [Failure evidence](failure_review.json).

Six jobs remain allocated at this snapshot: c00_h24, c00_h32, c01_h32, c03_h32, c04_h32 and c05_h32. Five show native progress; c03_h32 remains stalled during container verification, with process state unverified. Do not repeat unchanged node-access checks. Original deadlines remain October 5 02:12:43, except c00_h32 at 02:14:17 and c00_h24 at 02:53:41. Original paired 24/32, scalar-fit/final-reproduction and exact own-2015 prefix gates remain required before the second surprise. No two-shock 2023 state is accepted.

| Task | Preference | Horizon | Active job | Requested host |
|---|---:|---:|---:|---|
| c00_h24 | 0.13 | 24 | 19195622 | cs602 |
| c00_h32 | 0.13 | 32 | 19195626 | cs622 |
| c01_h24 | 0.14 | 24 | 19195630 | cs784 |
| c01_h32 | 0.14 | 32 | 19195632 | cs747 |
| c02_h24 | 0.145 | 24 | 19195636 | cs669 |
| c02_h32 | 0.145 | 32 | 19195638 | cs641 |
| c03_h24 | 0.14736308634876963 | 24 | 19195645 | cs680 |
| c03_h32 | 0.14736308634876963 | 32 | 19195648 | cs689 |
| c04_h24 | 0.15 | 24 | 19195660 | cs756 |
| c04_h32 | 0.15 | 32 | 19195663 | cs717 |
| c05_h24 | 0.16 | 24 | 19195664 | cs748 |
| c05_h32 | 0.16 | 32 | 19195665 | cs718 |


## Completed fixed-preference readout

| Fixed preference | 2008–2011 model (target 1.974875; weight 0) | 2012–2015 model (target 1.861; weight 1) | Fitted-window gap | Loss contribution |
|---:|---:|---:|---:|---:|
| 0.14000000000000001 | 1.7865183985 | 1.8027095598 | -0.0582904402 | 0.003397775415 |
| 0.14499999999999999 | 1.8269099482 | 1.8424294022 | -0.0185705978 | 0.000344867104 |
| 0.14736308634876963 | 1.8460211906 | 1.8610927710 | +0.0000927710 | 8.60645808e-09 |
| 0.14999999999999999 | 1.8673521968 | 1.8818053598 | +0.0208053598 | 0.0004328629976 |
| 0.16 | 1.9482310634 | 1.9591896869 | +0.0981896869 | 0.009641214609 |

The remaining 2016–2019 and 2020–2023 rows retain targets 1.755375 and 1.64575 and weights 0 and 1; model values, gaps and losses await the second surprise. Both preference parameters remain unestimated, with unchanged bounds [0.001789207206604163, 0.35784144132083257]. The full four-row table for every point, exact 2015 checkpoint pins and residual trajectories are in the lead review. Saved 2023 continuations belong to the first shock alone.
