# Social Security correction — September 11, 2026

The author requires Social Security to balance at every date. Existing matched
PF runs fix both payroll tax and the pension derived from a reference age ratio.
Their household, population and property-tax checks did not include a Social
Security balance condition. Existing quantitative results remain provisional.

## Historical reconciliation and restored repair direction

The author's recollection has direct support. On July 10 at 23:05:47 UTC the
assistant said pensions balanced within the stationary age distribution. At
23:08:33 and 23:12:39 UTC it explicitly described an external 17.9% tax and a
pension determined internally by stationary budget balance. The author asked
about precisely this distinction and then approved the source-table wording.
See `memory/transcripts/2026-07-10/combined_user_assistant.md:17586` and
`:17627`; the original response records are also retained in that day's raw
transcript folder. No instruction abandoning that condition was found in the
reviewed exchanges. These exchanges establish what was communicated; they are
not themselves numerical budget checks.

The July 15 mortality experiment supplies the concrete numerical link. Its
saved `report/lifecycle_decomposition.csv` has 12 equal-sized working cohorts
and five equal-sized retirement cohorts in the no-mortality control, M0.
The pension shortcut uses the ratio 12/5. Adding retirement mortality in M1
reduced the mass of retirees without updating that shortcut:

| Saved age distribution | Worker mass | Retiree mass | Demographic outlays/revenue factor, (12/5) R/W |
|---|---:|---:|---:|
| M0, no mortality before terminal exit | 0.705882353 | 0.294117647 | 1.000000000 |
| M1, post-retirement mortality | 0.740801474 | 0.259198526 | 0.839734374 |

The factor assumes actual mean worker earnings equal the pension formula's
reference mean. It is an accounting reconstruction from saved age masses,
not a recovered historical fiscal-residual receipt. The independent audit of
today's actual pre-announcement stationary payroll gives outlays/revenue
0.839734375: the mortality reconstruction differs by only 1.24e-9. This
identifies the stale no-mortality age ratio as the source of the current
pre-announcement stationary gap. Subsequent reweighting to the 2007 population
and changing transition demographics introduce further discrepancies.

The unchanged pension formula is visible in
`code/model/intergen_eqscale_seq_optimized/parameters.py:742`. Commit
`411616d22ce64634ca2ddd0d31a9a4dc8c27aa8c` records the July 10–16 machinery,
including mortality, without replacing this formula. The July 15 experiment
contract explicitly was not a production promotion; the exact date of its
later production adoption is not established here. Older MATLAB branches
also differed: January had an actual-distribution pension update, whereas a
March branch explicitly retained a fixed reference benefit. Do not conclude
that every historical version either balanced or failed.

The fiscal issue was subsequently flagged in August transition reviews and
the September 4 audit, but remained unresolved; see
`memory/transcripts/2026-08-11/combined_user_assistant.md:5281` and
`memory/transcripts/2026-09-04/combined_user_assistant.md:3740`.
The lead's initial claim
blurred the current audited failure with the older stationary specification;
the author was right to challenge that history. Arithmetic and source hashes
are preserved in `historical_reconciliation.json`; no model was solved or
changed for this reconciliation.

The repair direction is therefore to retain the external 17.9% tax and
determine pensions from actual payroll revenue and retiree exposure, consistent
with the previously communicated stationary specification and the author's
current requirement to balance Social Security. Extend the condition to every
transition date. Do not keep presenting the instrument as a wholly new,
unanswered choice. This is the repair direction, not a claim that the production
launcher or a corrected equilibrium has already adopted it.

In period units, with actual household distribution
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

1. Carry the recovered fixed-tax, endogenous-pension specification into the
   explicit endpoint and transition contracts; retain 17.9% as the external tax.
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
