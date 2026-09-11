# Social Security correction — September 11, 2026

The author requires Social Security to balance at every date. Existing matched
PF runs fix both payroll tax and the pension derived from a reference age ratio.
Their household, population and property-tax checks did not include a Social
Security balance condition. Existing quantitative results remain provisional.

## Economic decision still outstanding

The lead recommends fixed payroll tax with a common pension adjusting to actual
tax revenue. The author was also offered fixed pensions with the payroll tax
adjusting; no answer has yet been received. Neither rule has been adopted for
production. In period units, with actual household distribution
`g[wealth,tenure,location,age,z,parity,child]`, the identities are

\[
B_t=\Delta\sum_{i,j<J_R,z} w_i a_j z\,m_{tijz},\qquad
E_t=\sum_{i,j\ge J_R,z}[1+s(z-1)]m_{tijz},\qquad
\tau_t B_t=b_t E_t.
\]

Here `m` sums over wealth, tenure, parity and dependent-child states; `Delta`
is the model's flow period scale; `s` is its retirement-income dispersion
coefficient (zero in the common-benefit specification). Resident nonheads,
reference income-state weights and property-tax rebates do not enter payroll.
Both the initial stationary economy and the announced ACS-reweighted 2007 state
need their own balance checks; the latter has a different age distribution.

## Implemented and verified locally

Source commit **a654219c**, pushed branch `codex/balanced-social-security`, in
`tmp/e5f_matched_pf`. That worktree is clean. The production checkout and old
numerical snapshots were not modified by this implementation.

- `code/model/tools/e5f_social_security.py`: actual-household budget integration
  and idempotent dated income binding in period units.
- Optional dated pension/payroll-tax paths enter both backward values and
  forward policies, including the split historical/person-tail evaluator.
- `code/model/tools/e5f_social_security_root.py`: explicitly chosen fixed-tax or
  fixed-pension joint housing/budget root, separate gates and fresh replay.
- `code/model/tools/audit_e5f_social_security.py`: pinned, read-only accounting
  on the old initial and terminal distributions; no model solves.

**70 pure tests pass** (2.250 seconds of test execution): fiscal arithmetic,
income-state and period units, invalid inputs, separate residual gates, root
safeguards, cached-value mismatch rejection, and future-income anticipation.
The core root retains its previous behavior when no reset matrix is supplied.
The existing household evaluator retains its previous income when no fiscal
paths are supplied. All new production fiscal calls must explicitly supply a
path and validate the endpoint and source contracts.

## Cluster checks and actual budget gaps

| Job | Scope | Budget | Latest result |
|---|---|---|---|
| 17352552 | Read-only old initial and terminal budget audit; zero model solves | 1 CPU, 6 GB, 10 minutes | Completed in 13 seconds; collected |
| 17352615 | Real compiled household smoke: six two-date cases, 24 dated Bellman calls plus one tiny-grid fixture setup | 1 CPU, 8 GB, 20 minutes | Failed overall: two occupied-policy-response assertions; other checks pass |
| 17353361 | Same six conditional paths; inspect saving corners, age-specific response and all versus occupied state support | 1 CPU, 8 GB, 20 minutes | Completed in 11 seconds; collected |

The read-only audit verifies all 508 inherited source files and both endpoint
checkpoint hashes. All amounts below are in model period units; every original
state uses payroll tax 17.9% and pension 1.7184. Implied adjustments hold the
distribution fixed and are **not new equilibria**.

| State | Payroll revenue | Pension outlays | Revenue less outlays | Implied pension at fixed tax | Implied payroll tax at fixed pension |
|---|---:|---:|---:|---:|---:|
| Pre-announcement stationary | 0.530414 | 0.445407 | +0.085007 | 2.046361 | 15.0312% |
| Announced 2007 inherited state | 0.612854 | 0.299815 | +0.313040 | 3.512601 | 8.7569% |
| Terminal person/household endpoint | 0.342634 | 0.379912 | -0.037278 | 1.549786 | 19.8475% |

The 2007 state spends 48.9210% of payroll revenue on pensions; the terminal
budget falls short by 9.8122% of pension outlays. At fixed tax, its frozen-state
balancing pension would rise 104.4111% in 2007 and fall 9.8122% at the terminal
endpoint. These ratios use different explicitly stated denominators.

The compiled smoke executed all six tests in 12.244 seconds. It passed actual
first-date balance under both instruments, household budgets, mass, queues,
baseline reproduction and backward/forward policy replay. Future income did
change occupied initial values. Its saving-response assertions failed: occupied
current consumption/saving did not move (maximum saving/consumption change zero
for future pension and roundoff for future payroll tax). The original `.err`
and `.out` files are retained. No model gate was weakened.

An independent source check confirms that next-date, next-age values enter the
actual compiled optimization. However, a policy response is not universal:
the fixture forbids borrowing, has only one retirement age, and the two-date
pension change affects only its last working cohort. Diagnostic 17353361 is
collected in `anticipation_diagnostic.json`: conditional saving changes by up to
0.436831 under future pensions and 0.840777 under future payroll tax. Those
states have zero mass in this tiny fixture; occupied values respond while
occupied controls remain unchanged. This establishes that anticipated fiscal
income reaches actual policy optimization. It does not establish borrowing
corners as the sole explanation: the affected pre-retirement cohort has positive
saving. The original assertion requiring an occupied-control response is not
universal. The original failed suite remains preserved and has not been rerun
with a revised test; do not mark it as passed.

Both checks are short validations, not calibration searches. The compiled test
examines actual household budgets, first-date balance under each instrument,
both income paths' anticipation effects, policy replay, and birth-entry queues.
Changed cases hold prices and the terminal value fixed and are **conditional**.
The fixture has 20 wealth nodes, 6 ages, 2 income states and 2 owner products.
Expected runtime is several minutes including compilation; the hard cap is
20 minutes. No empirical target or loss is fitted.

Audit directory:
`/scratch/td2248/projects/Fertility_Spring26_paygo_audit_20260911a`.
Compiled snapshot:
`/scratch/td2248/projects/Fertility_Spring26_paygo_smoke_a654219c`.
The submitted scripts are retained alongside this README. The 4.3 MB source
archive was generated with `git archive a654219c code/model`; local and remote
SHA-256 both equal
`1d89372212a07af9114f9da29c3656245dd4ff70933d97ec2e225c6eeb84ba1b`.
The audit checks all 508 inherited source fingerprints, both checkpoint hashes,
and helper hash `f7bb7db2774be7ef7d2ca2a635eafa3cb27ba79a503fe69b6f00c1eef9658455`.

All three checks are collected and the monitor is paused. No model run,
production root or calibration is active in this repair. Original failures
remain preserved.

## Required before claiming the model is repaired

1. Settle the economic adjustment instrument with the author.
2. Complete the compiled test and inspect the actual saved-state budget gaps.
3. Rebuild and reproduce the initial normalization and terminal stationary
   endpoint under the chosen fiscal rule, with new source/input contracts.
4. Jointly solve dated housing prices and the pension or tax path, so both
   fiscal balance and household expectations use the same actual population.
5. Recheck horizon convergence and all existing accounting/numerical gates,
   then reassess historical preference fit and policy outcomes.

No corrected full equilibrium, new calibration or policy result is available
from this repair yet. Property-tax surplus/rebate rules remain a separate
fiscal contract; this correction does not finance pensions from that budget.
