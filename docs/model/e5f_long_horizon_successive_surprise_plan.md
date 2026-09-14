# Long-horizon refit with successive preference surprises

## Author decision and scope

The author confirmed successive surprises, each solved with a long forecast,
on September 13. Retain the calibrated 2007 structural parameters and initial
household distribution. Re-estimate four fertility-preference levels; the
currently running announced sequence is a separate comparison. The author
subsequently authorized execution. Exact-loop smoke job 17732511 was submitted
for batch `long_successive_refit_20260914a`. It passed and automatically
dispatched long-stage job 17732824; receipts are indexed below.

At each shock date households learn the current preference level and believe
it will remain permanent. They do not anticipate later surprises. Solve their
entire forecast, execute its first period, and carry the resulting households
and birth vintages into the next shock. Never initialize an intermediate date
from a stationary distribution or from a different fitted history.

## Targets and identification

Four free preference levels are matched to four dated fertility observations;
all other parameters remain fixed. The source is
`output/model/e5f_matched_pf_20260909a/path_pilot_20260910/fertility_data/empirical_blocks.csv`,
SHA256 `8945aa2427e26157e2e01b44b327e374daae5078daa458727f025938569a2f74`.
Its adjacent README and source manifest document published NCHS annual rates.

| Decision year | Birth years | Data | Model | Gap | Weight | Squared-gap contribution |
|---|---|---:|---|---|---:|---|
| 2007 | 2008–2011 | 1.974875 | Pending | Pending | 1 | Pending |
| 2011 | 2012–2015 | 1.861000 | Pending | Pending | 1 | Pending |
| 2015 | 2016–2019 | 1.755375 | Pending | Pending | 1 | Pending |
| 2019 | 2020–2023 | 1.645750 | Pending | Pending | 1 | Pending |

The scalar fitter uses the signed model-minus-data gap; unit weights above
describe a reporting loss, not empirical precision weights. Retain the existing
absolute fertility tolerance 0.005 and numerical search bounds
\([\psi_{SS}-0.20,\psi_{SS}+0.02]\). Do not impose an additional monotonicity
restriction on the four estimated levels. Save every estimate, bound distance,
trial and rejection. Reproduce the complete fixed-2007 target contract in the
final packet alongside these four new rows.

The empirical statistic is the arithmetic mean of four published female TFRs.
The retained model counterpart is `period_tfr_topcode_adjusted`, the existing
household-based age-specific fertility analogue. Preserve and label that
measurement approximation; it is not completed fertility or births per total
household. Pin its builder and normalization in the launch manifest.

## Solve and fit

1. Freeze the current original-queue numerical source and target fingerprint.
   Preserve the fixed housing supply curve, no immigration or rescaling,
   births/2.1 entry conversion, 1% annual tax with equal household rebates, and
   balanced PAYGO at payroll tax 0.179.
2. Start each date from its inherited household state. Use the earlier fitted
   level only as a search seed. For each candidate preference, solve and verify
   the stationary endpoint under that candidate's own permanent preference.
   Never use the final shock's endpoint for an earlier, different preference.
3. Solve the dated housing, pension and rebate equations, with household
   choices solved backward and the distribution carried forward. Use measured
   Jacobian information and the tested direction-preserving step only after
   source-matched verification on the long problem. Preserve all economic and
   numerical gates. The frozen eight-mapping per-call limit can be continued
   through explicitly budgeted rounds with saved best coordinates/Jacobian.
4. Use a bounded scalar bracket/secant search on first-period fertility,
   restarting from nearby converged numerical solutions. An invalid or
   unconverged equilibrium is not a usable observation for fitting. When the
   fertility gap passes, replay and checkpoint the actual next household state.
5. Repeat for all four dates. Historical stages are sequential. Trial candidates
   conditional on the same inherited state, endpoint solves, and sensitivity
   checks may run on separate cluster jobs without mixing their states.

Proposed common terminal date is 2423: 104 decision dates from 2007, 103 from
2011, 102 from 2015, and 101 from 2019. This leaves exactly 100 decisions from
the inherited 2023 economy for both policy and its already fitted control.
Preferences stay at the final fitted level after the 2019 surprise.

## Policy and numerical acceptance

After all four shock stages pass, the 2019 forecast generates the inherited
2023 distribution. Freeze that state before 2023 choices. Reuse the remaining
1% baseline forecast and solve an unexpected permanent 2% tax with equal
rebates over the same dates. Re-estimate neither preferences nor structural
parameters under policy. Recompute the policy endpoint at the newly fitted
final preference; today's 2% endpoint need not remain applicable.

Maintain three distinct tests: finite-horizon market/fiscal convergence;
distance of the carried terminal distribution, prices and birth queues from
the verified stationary endpoint; and stability of fitted shocks and policy
effects to a longer forecast. Passing the first does not prove the other two.
If endpoint distance fails, retain the finite result as provisional and extend
the horizon before claiming long-horizon adequacy. No convergence guarantee
is asserted. A final longer-horizon comparison must inspect early outcomes,
not merely whether the last plotted price reaches the imposed boundary.

## Readiness, cost and remaining implementation

The isolated driver is `code/cluster/run_e5f_long_successive_refit.py`; its
preparer is `code/cluster/prepare_e5f_long_successive_refit.py`. It retains the
old scalar search pattern but uses the current original-queue terminal and
three-block forecast adapters. Its policy mode reads the final fitted
preference and saved 2023 state directly and verifies the checkpoint hash.

The measured-Jacobian/scaled-step method passed a ten-period finite equilibrium
test (13 root mappings across two rounds); that does not certify 100 dates.
First run an exact-loop smoke and a bounded long-horizon first-stage candidate.
Verify target observation, replayed next state, checkpoint/resume, rejected
candidate handling, and the stable diagnostic packet before full dispatch.

Recent long mappings cost roughly 28–37 minutes. At 13 mappings, one cold
candidate would take roughly 6–8 hours, before endpoint and reporting work.
Four stages with three candidate values each would therefore be roughly
72–96 hours sequentially at that illustrative rate; this is a planning scenario,
not a measured forecast. Warm starts and parallel candidates may reduce it.
Do not promise an overnight fit from short-horizon timings.

Submission pins conservative caps: 30 minutes for smoke, one hour per terminal,
ten hours per candidate, four eight-mapping rounds, 24 hours per historical
stage, 12 hours for policy, and a seven-day absolute calendar expiry including
queueing. These are stop limits, not completion forecasts. The first long
candidate supplies the long-solver benchmark. The cap is six trial values per date,
with latest-completed and best-so-far files, per-mapping plots, minute
heartbeats and independent rejection receipts. A 30-minute period without
checkpoint or heartbeat triggers investigation. Smoke dispatches stage 0 only
after native stationary-root, changed-preference root-loop, target measurement
and checkpoint/replay checks pass. Each fitted stage dispatches the next;
the final one dispatches policy. Failed candidate evaluations are recorded,
but a failed fit cannot pass its state onward. Local unit checks passed;
compute-node smoke passed and the first full fitting stage is dispatched.

Pinned contract and receipts:
`output/model/e5f_original_queue_20260913a/long_successive_refit/`.
