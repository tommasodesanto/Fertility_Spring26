# Bounded numerical comparison — complete

**Both real evaluations passed, Torch 18753562.** The proposed joint step
reduced loss from 19.5813 to 13.7745/13.7789. The predicted normalization start
used three stationary solves versus seven and saved 50.56% of total time at
this point. All paired moment/price/psi screens passed. Early fertility remains
about 0.533 versus 0.810. No candidate is promoted. See [complete results,
targeted/untargeted fits and parameter bounds](RESULTS.md).

The preparation/first-arm entries below retain the launch history; RESULTS.md
is the completed comparison.

Reference: **2007 stationary reference — block0506, September 28 verified export**.
This folder owns only a numerical comparison on the existing one-birth model.
It does not run the separate two-birth experiment, resume an expired search,
change the frozen reference, or promote a candidate.

## What changes

The existing full-precision `measurement_audit_v1/bounded_experiment_plan.json`
supplies one damped Gauss–Newton proposal for all ten estimated coordinates.
Those parameter changes and the corresponding normalized child benefit are
**diagnostic recalibration**, not adopted estimates. Earnings, entry wealth and
income, timing, preferences' functional form, transfers and floors, targets,
weights, parameter bounds and scientific gates are unchanged.

Both arms evaluate exactly that proposal. The retained arm begins the fertility
normalization at psi `0.14281100340255604`; the other begins at the derivative
prediction `0.13052783066857948`. Both retain bracket step `0.005`, the same
initial prices, the existing within-objective price warm start, the completed
fertility target of 2.1 and demographic renewal. The adapted initial value is
only a numerical guess; it is not the final child-benefit value.

`adapter.py` is opt-in and changes only `normalization.initial_psi` in a deep
copy, **after** the original runtime factory authenticates its unchanged source
and original contract. It changes no model or shared runtime code. The new pair
contract and diagnostic identity pin that adapter independently; the inherited
contract is ancestry, not authorization of the new initial guess.

## Budget and stopping

Two fresh evaluator processes, two workers, at most 23 stationary solves and
1,800 seconds per objective, at most 46 solves total, with a 4,200-second global
budget. A controller-created process group receives SIGKILL at its deadline,
including descendants. No retries or extensions. An unknown or integrity
failure cancels any still-running sibling. Expected scientific rejection is
recorded, and the other arm can finish within its original budget.

The saved reference used six stationary solves and 970.676 seconds. A pair
running concurrently could therefore take around 16 minutes if the proposed
point behaves similarly; no speed gain is assumed. The global cap is 70 minutes.
`heartbeat.json`, `latest_completed.json`, `best_so_far.json`, `arm_failures.json`
and each arm's stationary-solve ledger preserve live progress. A comparison
point called best-so-far is not a promoted reference.

## Verification and comparison

The original native evaluator writes each successful case's full 14-row fit
table, 31-row parameter/bound table, checkpoint and unchanged 17 standard plots.
Checkpoints stay on Torch. The original validator verifies table arithmetic,
the complete target fingerprint, parameter bounds, solve accounting and hashes.
The pair wrapper adds its new source/adapter/plan identities, explicit economic
change labels and normalization/renewal/market cross-checks.

Additional paired screens are fixed before results:

- Loss difference at most 0.05.
- Relative price difference at most 0.0001 (symmetric denominator: the larger
  absolute arm price).
- Normalized psi difference at most 0.0001.
- Maximum scored-moment difference, multiplied by the square root of its
  working target weight, at most 0.01. This states the plan's “working scale”
  operationally without changing any weight.
- Each validation-moment difference at most
  `0.01 * max(1, abs(reference model moment), abs(target))`.

Passing these numerical screens does not establish statistical equivalence.
Failure means sensitivity to the numerical start; no equivalent-speedup claim
or promotion follows. `predictions_vs_actual.csv` keeps all 14 rows, with roles
separating scored, validation and normalization moments. `paired_comparison.json`
reports actual versus predicted losses, elapsed times and stationary solve
counts. Failed and censored arms retain their partial solve counts separately.

## Launch protocol

All model imports, verification, hashing, synthetic tests and rendering run on
Torch. `run.sh prepare` writes a new contract and refuses an existing one.
`run.sh smoke` requires its explicit SHA through
`EXPECTED_NUMERICAL_PAIR_SHA256`; it does not evaluate a model. The same finite
two-process dispatch, authenticated receipt path and process-group timeout are
tested with conspicuously synthetic fixtures. Placeholder checkpoints/PNG files
inside the smoke subtree are not economic results and must not be collected as
such. Tests include success, expected rejection, owned timeout with a child
process, fatal sibling cancellation, receipt corruption, all paired sensitivity
screens, target-fingerprint rejection and the copied one-field adapter.

Smoke job **18753074 completed successfully in 19 seconds** with two CPUs,
24 GiB and a six-minute cap. All 13 tests pass, with **zero model solves**.
The compact [smoke receipt](smoke_receipt.json) and [launcher log](smoke_launcher.log)
are retained locally. Complete synthetic artifacts stay on Torch.
**No main run was submitted.**

Pair contract SHA256:
`7950be909086744f25390211e5a7b93d7ea65f25d73698d68a90fa51ce303303`.
Original remote smoke receipt SHA256:
`c3fc2326b60bb464ef6382cbb5f0f134283e21982de5175a3eb641173e7d329c`.
The local compact receipt is a byte-for-byte copy of
`numerical_pair_v1/smoke_v1/smoke_receipt.json` in the Torch stage.

Main execution additionally requires a separate reviewed approval JSON whose
hash is supplied in `EXPECTED_NUMERICAL_PAIR_APPROVAL_SHA256`. Its fields are:

```json
{
  "status": "approved_bounded_numerical_pair",
  "pair_contract_sha256": "the reviewed contract hash",
  "smoke_receipt": {"path": "absolute smoke receipt path", "sha256": "its hash"},
  "launch_not_after_epoch": 0
}
```

Use an explicit future launch deadline after review, not the example zero.
Set `NUMERICAL_PAIR_APPROVAL` to that file and submit `sbatch run.sh run` in
the staged packet directory. The wrapper reads the supplied pins; it does not
mint approval. `run_v1` must not already exist. The pinned smoke is mandatory,
and approval expires rather than silently permitting a late restart. Further
exact-repeat verification and adoption are outside this pair's budget.

## Current limitations

After lead review of the sources, pinned contract and passed smoke, main job
**18753562** was submitted. It was pending cluster priority at the last check.
The separate `lead_review_approval.json` records a 30-minute launch window;
its SHA is `9491b4a72561dd2ce3886c7b2e4c4e49cfca3fd71241595531534b4f4f77e5a3`.
No automatic retry, window extension, final repeat or promotion is authorized.

The synthetic smoke tests orchestration and artifact validation, not economic
predictions, the native opt-in factory at a real candidate, or a speed gain.
No claim of improved calibration or faster normalization is supported until
the bounded real pair has run and its full outputs have been reviewed.

## First completed arm (paired comparison pending)

The predicted-start arm passed the native scientific gates and independent
controller/table/receipt checks. It used three stationary solves (513.486s
solving, 579.629s total), with loss 13.778878 versus reference 19.581311 and
linear prediction 13.335253. Final child benefit is 0.1311737861; completed
fertility is 2.100139 within the unchanged 0.0005 tolerance. This is a
diagnostic candidate, not an adopted reference or an exact-repeat certificate.

Early fertility is 0.533010 against 0.809528, slightly below the reference
0.535426. The calibration improvement comes from other rows. The older-age
wealth-dispersion validation moment worsens to 2.998638 against 3.515935;
this zero-weight row is retained. All ten coordinates satisfy their original
bounds; both fertility choice scales have the inherited near-lower-bound flag.

- [Complete 14-row fit, with scored/validation/normalization roles, gaps,
  weights and contributions](run_v1/predicted_start/case/target_fit.csv)
- [All 31 parameter rows, bounds, restrictions and near-bound flags](run_v1/predicted_start/case/parameters.csv)
- [Scientific and diagnostic receipt](run_v1/predicted_start/case/receipt.json)
- [Standard diagnostic packet](run_v1/predicted_start/case/standard_diagnostics/summary.json)

All 17 standard figures were collected and visually inspected. Housing-market
residual is 9.86e-7. Consumption/fertility wealth policies are regular at the
plotted scale; the previously noted high-wealth age-30 housing drop, small
ownership downturn, and retirement wealth kink remain. These observations do
not certify grid convergence or establish a new economic explanation. No
checkpoint was downloaded. Speed and numerical-start consistency await the
other arm's completion.
