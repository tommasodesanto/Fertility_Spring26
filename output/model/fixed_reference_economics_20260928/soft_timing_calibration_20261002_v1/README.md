# Matched soft purchase timing calibration

This packet prepares four paired starts for each of two transaction clocks. The
first arm continues the original interest timing; the second estimates the
isolated post-interest transaction timing. Both retain the author's adopted soft
purchase financing rule, the physical parent room floor, constant child utility
share, no age-varying utility term, nonnegative-mean entrant distribution with
no unsecured debt, financed share `phi = 0.8`, normalized population `N0 = 1`,
and the same ten free coordinates. The complete 14-row target system and weights
are pinned before a solve. The post-interest arm is experimental until selected
native checks and author review; no calibration is promoted by this driver.

`driver_plan.json` records the exact four starts, parameter bounds, source and
target fingerprints, and budgets. Every chain runs in a separate fresh Python
process. The original arm never imports the alternative household overlay. The
alternative arm uses the authenticated sandbox sources, calendar forward map and
independent purchase-accounting audit already exercised by the fixed-coordinate
experiment. Both use the normalized price/renewal root and derived physical
housing supply coefficient at `N0 = 1`.

The driver is
`code/model/experiments/purchase_timing_sandbox/calibrate.py`. Its interface is:

```bash
code/model/.venv/bin/python code/model/experiments/purchase_timing_sandbox/calibrate.py \
  --arm original --chain 0 --out /absolute/new/chain/output \
  --deadline-epoch UNIX_EPOCH_SECONDS
```

Use `--mock-smoke` for two exact optimizer-loop iterations and checkpoint writes
with zero model solves. Use `--preflight-evaluator` to construct and authenticate
the exploratory search evaluator with zero solves. Use `--smoke` for two
exploratory objective calls followed by a full native selected-point verification.
The smoke checks all 14 target rows and 31 parameter rows against the native
report at the existing absolute tolerance of `1e-10`. Search and final
verification use the same evaluator modes in smoke and production. The full
native selected-point evaluator runs in a fresh child interpreter under the
same absolute deadline; the parent admits its receipt only after checking the
search receipt hash, target and weight fingerprints, tables, plots and exact
repeat. A production chain allows 250 objective
calls, at most 32 lifecycle solves per case under the native evaluator, six
hours from actual start, and 30 minutes reserved for native verification. Each
case writes `latest_completed.json`, `best_so_far.json`, and `heartbeat.json`;
the provisional search writes `search_completed.json`. Only a full selected
postcheck with 14 fit rows, 31 parameter rows, 17 standard plots and exact
repeat can write `selected_numerically_verified` to `completed.json`. The
controller must investigate a missing heartbeat after 30 minutes and must not
automatically retry, change targets, relax gates, or promote a provisional
search result.

The common first start is the saved soft chain 16, case 0046 point with loss
23.07830929416 under the original clock. The prior fixed-coordinate
alternative timing loss is 57.18633285456. These are starting evidence, not
the results of this matched calibration. The three additional starts are small
deterministic relative perturbations, paired exactly across clocks. The plan
records their full numerical values and unchanged search bounds, including the
author-requested `h_P <= 2.6`.

## Expanded matched starts

The optional `expanded_start_plan.json` supplies 24 deterministic starts for
**each** timing arm. Its first four vectors exactly match the eight already
running legacy chains across both arms. The 20 additional vectors per arm are
seven other saved soft-checkpoint best points (ranks two through eight), eight
modest perturbations around the selected point, and five broader stratified
points. The last two groups use independently shuffled Latin-hypercube
midpoints per coordinate with seed `20261003`. All 24 vectors are unique,
inside the unchanged ten-parameter search bounds, and paired across timing
arms. Per-start provenance and broad sampling ranges are in the JSON; those
sampling ranges do not change optimization bounds or economics. The historical
anchors retain the same 14 targets and weights.

Pass `--starts-file /absolute/path/expanded_start_plan.json` and
`--starts-file-sha256 383d3e6833082a0796985d0bf65c017eafcda98d4ad4fa9c401c3a6412706270`
with `--chain 0` through `--chain 23`. The table hash is checked before model
setup and passed unchanged to the native postcheck child. Without these flags,
the four-start behavior remains. Deploy only indices 4–23 in each arm (40 new
chains); preserve the eight running earlier chains.

At existing caps, all 48 chains represent at most 288 CPU-hours and 12,000
objective calls; the 40 additional chains account for 240 CPU-hours and
10,000 calls. Running all 48 simultaneously would reserve 1,152 GiB at
24 GiB each. The formal 32-lifecycle-solve cap per objective gives a loose
384,000-solve ceiling; seven solves per objective would be 84,000. These are
arithmetic budgets, not forecasts: the six-hour chain deadline can stop
searches earlier. The staged launch controls actual concurrency. An earlier
120-start proposal was superseded before submission and retained only as
`driver_verification/superseded_120_start_plan.json`.
