# First-child housing loading without Stone-Geary

September 26, 2026. **Mechanism diagnostic completed; no new calibration or
production specification adopted.** The author requested a conventional
utility formulation with an explicit first-child housing loading. The tested
loading applies when at least one child is at home and is the same for one,
two or three children at home. The equivalence scale and child benefit still
vary with the number of children.

Primary readable report: [PDF](../../pdf/first_child_loading_probe.pdf).
It contains the specification, summary, incentives, all target and parameter
tables, and the unchanged 17-figure diagnostic set for each of five cases.

| Specification | First-child housing loading | First-birth housing proxy (rooms) | Fertility | Childlessness (%) | Descriptive loss |
|---|---:|---:|---:|---:|---:|
| Original floor, linear child benefit | Floor | 0.785 | 2.100 | 18.279 | 285.640 |
| Floor, concave child benefit | Floor | 0.783 | 2.042 | 19.206 | 296.264 |
| No floor, concave benefit | 0 points | 0.120 | 2.633 | 3.267 | 3707.254 |
| No floor, concave benefit | 10 points | 0.941 | 2.605 | 3.737 | 2871.895 |
| No floor, concave benefit | 20 points | 1.538 | 2.540 | 4.972 | 2084.010 |

The first-birth housing target in the retained comparison is 1.465 rooms.
An explicit first-child loading can generate that magnitude without a minimum
housing constraint. The zero-loading diagnostic is weak at this fixed point;
positive loadings materially increase the model housing response. This is
evidence for the proposed mechanism, not empirical validation of the particular
loading or the full model. The statistic is still an unmatched stationary proxy
for the empirical event study.

The other moments do not jointly fit. For the 20-point loading, fertility is
2.540 versus 2.100, childlessness is 4.972% versus 19.828%, first-birth age is
24.210 versus 25.976, and the recent-parent ownership gap is 0.404 percentage
points versus 12.761. Mean rooms is 6.769 versus the frozen old ACS target 5.608.
The adopted AHS target is deliberately not substituted into this controlled
comparison. A larger housing response cannot stand in for these failures.
The original floor/concave point also misses the housing-response target.

## Tested preferences

Let c be nonhousing consumption, s housing services including the retained owner
premium, and m children currently at home. Material consumption is C. Utility is

\[
u(C,m)=\frac{C^{1-\sigma}}{1-\sigma}
       +\psi\frac{m^{1-\kappa}}{1-\kappa},\qquad
C=A(m)\frac{c^{\alpha(m)}s^{1-\alpha(m)}}{e(m)},
\]
\[
\alpha(m)=0.733-\lambda\mathbf1\{m>0\},\qquad
e(m)=\left(\frac{2+0.7m}{2}\right)^{0.7}.
\]

Sigma is material curvature, retained at 2. Kappa is child-benefit curvature,
fixed at 0.140 as a literature-guided experimental restriction in every new
case. It is distinct from the fertility taste-shock scales and is not estimated.
The original control retains kappa=0. The code coefficient b is held at the
original first-child benefit; therefore psi=(1-kappa)*b in the displayed CRRA
form. This preserves the benefit at m=1 while changing benefit curvature.

Lambda is a first-child housing-weight increment, fixed at 0, 0.100 or 0.200.
The positive cases raise housing's weight from 0.267 to 0.367 or 0.467 whenever
children are at home. Subsequent children do not add another direct housing
loading. In a future calibration, this parameter can replace the old housing
floor parameter; neither these trial values nor the curvature are estimates.

A(m) is the inherited substantive fixed normalization, not an estimated
parameter. Define K(a,r)=a^a*((1-a)/r)^(1-a). The normalization is
A(m)=K(0.733,r*)/K(alpha(m),r*), where r* is the fixed reference rent in the
parameter table. It equalizes compensated expenditure for uncapped renters
at that reference rent, before applying the family needs scale. It does not
equalize every owner's or capped renter's utility. The unchanged owner premium
also receives a different weight when alpha changes; that is a consequence of
the proposed preferences, not a separate altered primitive.

## Changes, closure and scope

Experimental changes relative to the named floor/linear control are explicit:
child-benefit curvature is fixed at 0.140 in all new points; the share cases
also remove the housing floor and impose the normalized first-child share
loading. The floor/concave case isolates the benefit change. The zero-loading
share case helps distinguish floor removal from adding the housing loading.

Earnings, initial wealth/income, timing, mortality, bequests, credit, fertility
and tenure shocks, housing products, supply, fiscal rule, targets, weights and
numerical gates are unchanged. There is no transfer floor or public assistance.
Child-benefit levels are not renormalized and no structural parameter search is
performed. Housing and pension markets/accounts clear under normalized entry.
The demographic replacement gaps are 2.767%, 20.251%, 19.401% and 17.327% across
the four new cases. These are not closed-renewal production equilibria.

Lead interpretation: the conventional first-child share specification is a
credible candidate for joint calibration. This diagnostic rejects neither the
floor model nor the share model as a calibrated architecture. Further utility
complexity is not justified by these results; the outstanding task is fitting
the housing and fertility evidence together with consistent measurement.

## Full evidence and reproduction

- Each case folder has `target_fit.csv` (all 13 rows: target, model, gap, weight,
  loss contribution), `parameters.csv` (all 32 parameters/restrictions, bounds
  and near-bound flags), incentive summaries, and scientific receipts.
- [Design and authenticated source pins](design.json), [case-loop smoke](smoke.json),
  [complete results](complete.json), [final verification](final_verification.json).
- Standard figures are embedded unchanged in the PDF; originals and all four
  large solution checkpoints are retained on Torch, not copied to the laptop.

Torch job18604498 completed four new stationary housing equilibria and the
authenticated saved control in 631 seconds after setup. Every case passed the
unchanged scientific gates. The control target table is identical to its saved
reference. All 65 target rows, 160 parameter/restriction rows and 85 standard
figures are present. Target/weight equality, gaps, loss contributions, fixed
first-child benefits and actual endogenous pensions were checked independently.
The exact configuration loop was smoke-tested before the first new solve.

Later report-only revisions improve normalization wording and labels, add an
incentive table and render checks, and classify unexpected errors honestly.
They do not change any solved case. The final PDF was rebuilt on Torch; its
69 pages were rendered, visually reviewed through contact sheets and full-size
representative pages, and checked for text outside page boundaries.

Source: `code/model/tools/run_e5f_first_child_loading_probe.py`, using the
reporting helpers in `code/model/tools/run_e5f_soft_housing_probe.py` and the
authenticated frozen utility-comparison runtime. The new optional reporting
arguments preserve the original soft-probe defaults. Executed/final driver
hashes are recorded in `final_verification.json`.

Remote root:
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/first_child_loading_probe_20260926`.
Retained run: `run_001`. Runtime/reference contract remains the pinned
`utility_four_arm_preparation_20260925_v2` comparison.

To regenerate, use a Torch allocation with one CPU, 12 GB memory and a 25-minute
limit. Load `anaconda3/2025.06`, set OMP/OPENBLAS/MKL/NUMBA thread counts to 1,
set MPLBACKEND=Agg and PYTHONDONTWRITEBYTECODE=1, and use
`PYTHONPATH=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/utility_four_arm_preparation_20260925_v2/tools:/scratch/td2248/commute_pdf_qa_deps`.
Run the remote `tools/run_e5f_first_child_loading_probe.py --output <fresh-output>`.
To regenerate the PDF and visual QA only, add `--render-only` and point
`--output` to the retained run. That path performs no model solve.
