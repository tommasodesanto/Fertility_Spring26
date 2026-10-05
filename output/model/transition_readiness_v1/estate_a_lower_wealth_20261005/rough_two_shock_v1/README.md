# Rough sequential transition from the lower-wealth Estate-A estimate

## Current facts and scope

Status: submitted October 5 at 14:21:58 New York as job **19238662**, now
**RUNNING on cs604**. The first fixed-price bridge solve completed in
27.854843792 seconds with one lifecycle call and zero entry-censored mass;
the bridge's downstream accounting/array checks are still in progress, so no
completed bridge or two-shock fitted result is claimed. See
`deployment/submission.json`, `deployment/launch_status.json` (initial pending
state) and `deployment/running_status.json`. The user queue was empty before submission.
The persistent submission receipt and task claim prevent an unknown-outcome retry.

Remote job root:
`/scratch/td2248/projects/estate_a_lower_wealth_rough_transition_20261005_v1`.
Package SHA-256:
`ad68ed1b46c82eab45745226152b75c3aea329adb2f8cbb07a36859c12072d0f`.
Launcher SHA-256:
`e0dca9d99a5b43af1ae522337a69d709895e05214d490bd44a35aeded9a70072`.
The zero-call Torch reference preflight passed, including selected report and
repeat-array pins. Six focused local fake-loop/contract tests passed; the lead
review is `lead_execution_review.json`. The native loop check runs inside the job.

Tommaso's October 5 messages in `Reconcile yesterday’s threads`
(`01a10804-6aad-7072-bd35-a65fd4f080e3`) prioritize estimating the transition,
permit a rough first path with gentle convergence, and request the new estimate
using the lower wealth target. These human messages were read directly. This
is a separate rough estimation experiment, not approval of the old six-job H32
budget proposal.

Reference: Estate-A overnight array 19194495, chain 11, case `0149_nm`, saved
loss 15.021541735257825. It remains provisional pending the calibration's final
selected-point postcheck. The complete 14-row target/fit table and 31-row
parameter table are retained at:

- `output/model/overnight_dual_20261004_v1/terminal_collection/a_provisional/target_fit_new_contract.csv`
- `output/model/overnight_dual_20261004_v1/terminal_collection/a_provisional/parameters_estate_a.csv`
- `output/model/overnight_dual_20261004_v1/recovery/inputs/a11_best_so_far.json`

The actual remote selected report is
`/scratch/td2248/projects/estate_birth_overnight_20261004_v1/results/production_binary_chain_11/run/0149_nm/phase_b_ge/selected_root`.
The source target fingerprint is
`c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70`;
the weight fingerprint is
`f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4`.

The new calibration has population normalized to one, price
0.7551066122345627 and derived housing-supply coefficient
\(H_0=6.299624572680468\). Freeze that coefficient for the transition and use
the candidate's own physical household distribution. The checkpoint's
`fixed_h0_population_scale=1.016837356224955` instead refers to the inherited
housing coefficient and must not be combined with the newly derived coefficient.

Relative to the previous chain-13 transition, the selected parameter point
changes all ten estimated coordinates, derived housing supply and the derived
benefit coefficient. The new point was estimated under the lower wealth target
4.45838713455674. The old Estate-A report already *rescored* its old point using
that same lower target; the four transition fertility targets do not change.
Retain one birth per period, post-interest timing, soft financing, net estates,
fixed payroll tax and endogenous pension. No additional earnings, credit or
birth-cap changes are authorized here.

## Separate rough numerical controls

Use 16 four-year periods (64 years) per stage and a freshly measured 12-date
seed, perturbed at date 5 with log step 1e-5 and five native mappings. Rough
housing tolerance is 0.005; scaled fiscal tolerance is 0.001. Keep native replay
1e-10, inherited infeasible-mass tolerance 1e-12, stationary renewal 1e-6,
nonnegative population and original accounting checks. These controls do not
alter the existing production tolerances or certify horizon stability.

Fit the unexpected permanent 2007 preference to 2012–2015 fertility 1.861.
Implement exactly the first two periods, using their own 2015 boundary price,
pension, value function, physical population and both birth queues. Then reveal
the unexpected permanent 2015 preference and fit 2020–2023 fertility 1.64575.
Initialize the second stage from the accepted first-stage forecast slice
`[2:14]`, without padding. Require informative scalar responses, fitting gap
at most 0.005, and fresh selected reproduction.

Both absolute preference bounds remain
`[0.001789207206604163, 0.35784144132083257]`, with starts
`0.14736308634876963` and `0.12`. They are explicit bounds, not newly computed
multiples of the new reference preference `0.176948250201189`.
The four fertility rows are `[1.974875, 1.861, 1.755375, 1.64575]`, with weights
`[0, 1, 0, 1]`; the historical path uses the first two rows of each stage.

Initial execution envelope: one CPU, 24 GiB, four hours externally and 14,280
seconds internally across input preparation, a small execution check and both
fits. This is below the previously authorized six-hour empirical envelope.
Keep 20,000 actual calls, 12 scalar evaluations per stage, 12 path evaluations
including replay, 48 endpoint evaluations, 1,800-second seed/endpoint/map caps,
6,000-second path cap and 7,200-second candidate cap. No automatic extensions
or unknown-outcome retries. A stopped fit remains incomplete.

Rebuild the required initial state at the saved price if the array-only search
archive cannot provide the full native solution object; no new stationary price
search is needed. This bridge and the small execution check belong to the same
bounded job. Exact deployment pins and submission evidence will be added under
`deployment/` once prepared and reviewed.

Save compact parameter/loss/status records for every trial, plus full best,
latest and selected states. Generate the stable 17-plot packet for selected
results. Superseded root-map diagnostic packets are removed after a newer map
completes; superseded candidate diagnostic packets are removed after a newer
candidate completes, retaining best/latest/selected. Numeric records remain.
The job checks for at least 32 GiB scratch quota headroom initially and monitors
its own output against 16 GiB every 60 seconds; this sampled limit can overshoot
between checks. It stops itself on an observed excess or write failure.
Preserve all preexisting transition packages and evidence. Never use
the old H24/H32 measured seeds, endpoints or optimizer matrices as inputs here.
Final rough output must disclose provisional calibration, untested strict
horizon comparison and terminal diagnostics; neither a close first-stage grid
point nor its 2023 continuation is a completed two-shock estimate.

Storage prerequisite was checked at 13:55 New York on October 5: a 4 KiB
write/fsync/remove probe succeeded under the project scratch root; quota reported
4.63 TB used of 5 TB. Cleanup belongs to its separate owner and is not repeated
by this experiment.
