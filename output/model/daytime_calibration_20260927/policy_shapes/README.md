# Supplemental policy-shape inspection

Saved selected calibration checkpoint only; zero Bellman, stationary or forward
kernel solves. `inspect_e5f_policy_shapes.py` uses the authenticated household
loader and preserves the standard 17 figures. `source_receipt.json` identifies
the exact checkpoint and target contract. This is a bounded first inspection,
not certification of all policies.

The frozen diagnostic implementation mixes conditional renter consumption,
modal-tenure housing, and probabilistic ownership. These are not a single
realized policy. The supplemental plots distinguish them and construct expected
housing/consumption by the saved six tenure probabilities and transaction-map
interpolation, conditional on the childless incoming-renter state (one location,
child state0). Fertility is not integrated out. CSV includes all15 income states;
plots select actual ascending income ranks1,8,15. Horizontal range is population
post-fertility/pre-tenure wealth p0.1–p99. This is a different timing from the
post-tenure financial-wealth distribution.

## Findings

- High-wealth ownership decline is off occupied support at extreme wealth:
  both examined slices have zero mass above300. At wealth2673.753 the middle
  income ownership probability is .833(age30) and .832(age42). Renter probability
  is close to1/6 and owner product probabilities become less concentrated,
  consistent with fixed logit shocks dominating shrinking value differences.
  They are not exactly equal shares; no asymptotic theorem is claimed.
- The large standard age30 housing downturn is the highest-income modal choice
  changing from10rooms to renting6rooms at the top wealth node3000. That node
  has zero mass in these slices. The smaller occupied-region features remain.
- Ownership is not globally increasing with wealth even inside the zoom. With
  adjacent-node mass each >1e-10, decreasing-pair upper nodes contain0.086% of
  the examined age30 slice and8.826% of the age42 slice. Their whole-population
  masses are0.00135% and0.0729%. These are descriptive grid-node metrics, not
  causal or lifetime incidence; they do not establish a numerical error.
- Expected housing is considerably smoother than modal housing. Decreasing-pair
  upper nodes contain about0.570% and0.541% of the respective slices. Expected
  consumption is mostly smooth; age30 has a small .033 decline at negative
  wealth with meaningful mass. Savings can decline as tenure lotteries change.
  Housing discreteness, saving-grid switches, and endogenous tenure mixtures
  must be distinguished before classifying these features as bugs.

`policy_slices.csv` contains4800 rows; `shape_audit.json` includes all15 income
slices and tail tenure probability vectors. `occupied_declines.csv` conditions
on both adjacent nodes having mass >1e-10 and reports the denominator and worst
local decline. The latter is a deterministic CSV reduction of policy_slices.
No conclusion about all ages, parents, or owners follows from these two renter
slices. All-grid probability sums include zero-filled infeasible states; the
script checks finite probabilities everywhere and unit sums only on occupied
post-fertility states. An initial unconditional unit-sum assertion failed for
that reason; no numerical model or gates were altered.

## Lead verification and remaining diagnosis

Lead reviewed the extraction against the native transaction-map kernel: lower
node weight1−w and upper-node weightw, with saved tenure probabilities. All10
occupied-decline summaries independently reproduced from CSV to1e-12; their
reducer is now included in the reproducible script. Both supplemental figures
visually inspected. This is saved-array extraction, not fresh solution accuracy
verification. Conditional renter controls outside a feasible renter branch
should not be interpreted as realized consumption.

Concrete age42 income-rank8 interval: wealth1.419→4.209, ownership.484→.361,
expected housing4.545→5.593. Four-room owner probability falls.479→.205, while
six-room owner probability rises.005→.154 and renter probability rises.516→.639.
This identifies the choice substitution underlying the aggregate decline; it
does not yet determine whether its magnitude is economic or grid-driven.
No change to tenure scale, housing grid or borrowing rules has been adopted.

Reproduce (one process, zero solves):
```sh
EXPECTED_UTILITY_OVERNIGHT_SHA256=3b770d8c8c22d2b0449b34a575d6353b063bc015d74ce11016dad7e22ed7ca5e E5F_LOCAL_EXECUTION_AUTHORIZATION=tommaso_authorized_20260927_local_primary_continuation_v1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MPLBACKEND=Agg code/model/.venv/bin/python code/model/tools/inspect_e5f_policy_shapes.py --contract tmp/e5f_overnight_local_20260927/portable/night_launch_v4/primary_continuation/production_contract.json --case output/model/daytime_calibration_20260927/local_continuation/run_v1/combined_scale_curvature/case --output output/model/daytime_calibration_20260927/policy_shapes
```
