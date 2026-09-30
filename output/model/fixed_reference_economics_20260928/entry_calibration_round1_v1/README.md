# Four-hour entrant-wealth calibration, first round

## Verified specification and source facts

This packet is an author-authorized experimental continuation of the verified
`entry_calibration_pilot_v1` results, not a production calibration or adoption.
The run has three parallel lanes:

| Lane | Entrant wealth | Renter borrowing floor | Grid |
|---|---|---:|---|
| `empirical_credit_120x9` | Existing five PSID ratio-bin means times current annual gross entrant earnings | −0.25 | 120 assets, 9 income states |
| `nonnegative_mean_120x9` | Negative ratio-bin means set to zero; positive nodes scaled to preserve the original mean | 0 | 120 assets, 9 income states |
| `nonnegative_mean_160x15` | Same nonnegative five-node rule | 0 | 160 assets, 15 income states |

Credit is measured in economy-wide mean annual gross working earnings units.
Positive-node scaling is 0.3632385158888715. This transforms the existing five
bin means; it is not censoring individual raw survey records. Entry projection
is the existing native linear grid projection and must preserve the mean with
zero clipping in these lanes. No entrant mass is discarded or moved to repair
feasibility.

All lanes retain the common 2% annual real saving/borrowing rate, reviewed
collateralized homeowner financing and mortgage-sale repayment, mortality
repayment rules, earnings/tax/pension objects, preference functional form, targets and weights.
The experimental economic changes relative to September 28 block0506 are the
corrected no-taper/full-sale credit rules, entry-wealth law and borrowing floor,
and the birth-renewal price/population closure already used by the pilots.
The coarse/fine distinction is numerical income/wealth discretization.

Price clears actual birth renewal; population clears absolute housing supply.
Housing supply scale H0 and child-benefit level psi_child remain fixed. Nine
historical coordinates are searched against ten scored rows; all 14 target rows
and 31 parameter rows/bounds are reported. There is no benefit normalization.
The fine lane loads the actual authenticated original 160×15 bundle. It does
not interpolate the 120×9 input distribution or merely change dimension labels.

Each lane starts at its corresponding verified pilot's selected parameters.
Both nonnegative lanes use exactly the same nine-coordinate starting vector.
`plan.json` freezes those values, full selected parameter tables, pilot labels,
losses and source/repeat receipt hashes. The original checkpoint constructor
defaults are never substituted. Source and bundle hashes are authenticated.

## Bounded search

The total lane clock is 14,400 seconds, including launcher preflight and final
verification. `--deadline-epoch` passes the launcher clock into the native run.
One CPU and one computational thread are required on Torch. The hard evaluation
caps are 80 full GE evaluations for coarse lanes and 36 for the fine lane,
including baseline/repeat and the two selected-point full-GE repeats.

After repeated native baseline GE, each round forms nine fresh one-sided finite
differences around the current best admissible point. The finite-difference
size and coordinate transformation, ridge rule and componentwise trust radius
are exactly those of the pilot: ridge is the largest singular value squared
times 1e-4 (minimum 1e-10), trust is three finite-difference steps. Bounded
Gauss–Newton proposals use dampings 0.5, 0.2 and 1.0. Adding the full damping and
refreshing the center/Jacobian are numerical search changes; no solver or
objective edits accompany them. There are at most eight coarse/four fine rounds.

A partial Jacobian near the time/evaluation cap is explicitly incomplete;
missing columns are never replaced by zero. A full rank-nine Jacobian is required
for a Gauss–Newton proposal. Rank deficiency or no local improvement stops the
search. Successful derivative probes may improve the selected point. This is a
local preliminary search; completion is not a convergence or identification
certificate. The latest identification receipt is at its recorded round center,
not necessarily at the selected point; prior complete receipts remain in
`rounds.json`.

Reserve is the maximum of 1,200 coarse/1,800 fine seconds, three observed GE
seconds, and 710 plus 2.6 observed GE seconds. Conservative observed floors are
210 coarse/600 fine seconds. Exploratory evaluators receive an earlier deadline
that protects this reserve; selected repeats receive the global deadline. Exact
native time-reserve exits stop exploration and retain the best completed point.
Accounting, feasibility, parameter, source and per-case timeout failures halt the
lane. Only explicit unbracketed exploratory renewal roots may be rejected as
inadmissible. No retries, relaxed gates, extensions or specification fallbacks
are introduced.

Native GE retains its 20-lifecycle cap and 300-second per-lifecycle/observer
clock, original safeguarded bracket and native repeat gates. Each GE includes
its selected-price native repeat; the final selected point additionally needs
two independent full-GE repeats, exact 14/31 tables, closure and 17 PNG hashes.
Insufficient verification yields `provisional_budget_exhausted`, never a verified
calibration. Latest/best/case receipts are written at every evaluation and native
GE writes its own point-progress receipts.

## Verification and commands

The 13 pure/mock tests pass with zero model solves. They check true lane grids,
seed provenance, all nine bindings and fixed objects, target/dimension contract,
multiple current-best updates, partial/failed/rank-deficient Jacobians, baseline
and selected-repeat failures, fatal accounting errors, guarded search-budget
stops, reserved repeat deadlines and native phase-B 31-field stale-binding tests
for both grids.

```sh
cd output/model/fixed_reference_economics_20260928/entry_calibration_round1_v1
../../../../code/model/.venv/bin/python -m unittest -v test_round1
python runner.py --mode preflight --lane nonnegative_mean_160x15 --out /new/path --deadline-seconds 300
```

Exact mocked launcher preflight performs the actual controller sequence with
zero lifecycle calls: 80/80/36 mock GE calls for empirical coarse, nonnegative
coarse and nonnegative fine. Caps leave a partial last Jacobian and preserve both
selected repeats. Native baseline/repeat gates still execute in the actual jobs;
mock preflight does not certify a native equilibrium or prediction of fit.

Results and launch receipts will be indexed here after dispatch. Use the existing
standard diagnostic packet; do not generate unsolicited PDFs or copy large model
arrays during collection.

## Launch — September 30 evening

Torch array **18899847** submitted for tasks0/1/2 in the table order above.
Preparation array18899803 passed all three exact CLI preflights (80/80/36
mockGE, two repeats each, zero lifecycle solves) and13 tests per task. Lead and
independent review passed source/binding/budget checks;86pins authenticated.
The source archive has88files,355268bytes, SHA
`4657194f2dd8cad599aa491b35b97a7ace28583a7fa2ecf10b8bef9a1331937e`.
Native repeated-baseline gates execute inside each four-hour job.
See `launch.json`, `lead_verification.json`, and `deployment/preflight_receipts/`.
The existing thread monitor follows this array every15minutes, notifying only
meaningful failures/actions or completion; no automatic restart or extension.

All three tasks started at18:38:17NewYork, with a22:38:17hard end.
All built-in preflights and native entrant-income/mass checks passed; each is
solving its baseline. This is launch verification, not a completed calibration.
