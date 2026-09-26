# Estate recipient diagnostic

## Receipt-risk finding — completed September 26

**Inheritance uncertainty matters.** The matched fixed-price test completed in
Torch job 18571183 in 3 minutes 35 seconds. It compares no receipts, a certain
age-specific payment, and a zero/positive lottery with exactly the same
age-specific mean as the certain payment. Prices, the receipt-profile scale,
preferences, entry, earnings and the inherited fiscal/target contract are fixed.
All receipt-law changes are experimental and unadopted.

The numerical follow-up completed in job 18571631 in 11 minutes 4 seconds,
using 319 rather than 160 wealth nodes. It includes an independent unpatched
control plus the same three-case comparison. The table below uses this finer
grid; the [original-grid readout](results/receipt_risk/v1/report/README.md) is
retained separately.

| Moment | Target | No receipts | Certain mean | Lottery |
|---|---:|---:|---:|---:|
| Ownership, ages 30–55 (%) | 67.626 | 75.555 | 73.461 | 75.661 |
| Recent-parent ownership gap (percentage points) | 12.761 | 7.984 | 3.036 | 6.369 |
| Childlessness (%) | 19.828 | 17.502 | 14.544 | 16.177 |
| Fertility | 2.100 | 2.126 | 2.217 | 2.166 |
| Wealth / annual earnings | 6.927 | 6.054 | 5.961 | 6.251 |
| Mean rooms | 5.608 | 6.553 | 6.724 | 6.707 |

The complete [comparison and visual packet](results/receipt_risk/finer_grid_v2/report/README.md)
includes [all target fits, weights and loss contributions](results/receipt_risk/finer_grid_v2/report/full_target_fit.csv)
and [all parameter estimates, bounds and restrictions](results/receipt_risk/finer_grid_v2/report/parameters.csv).
The finer-grid fixed-parameter losses are 302.060, 624.166 and 424.947;
the original-grid losses are 281.190, 584.392 and 388.612. They
describe fit under the old experimental target contract, not model selection
or a new calibrated specification.

The original-grid no-receipt control reproduces the retained net-valuation
arrays, pre-choice distribution, price and every target exactly. The finer-grid
control reproduces independently solved unpatched arrays and its calendar
distribution exactly. All 13 current focused tests, the ordered-loop smoke,
and original budget, purchase, fiscal, calendar, value and probability gates
pass. Occupied receipt clipping is zero. The 102 standard diagnostic figures
across both three-case packets are retained unchanged and have been visually
inspected through complete contact sheets.

The matched test establishes that replacing uncertain receipts by their
conditional mean substantially changes household behavior. It does not
attribute every difference from the earlier equal age46–62 payment to risk:
that earlier experiment also had a different age profile, timing and endogenous
prices. The recent-parent statistic is a descriptive ownership gap, not a
causal birth effect. The lottery improves the wealth/earnings comparison but
leaves material housing-size, childlessness and parent-gap discrepancies.

**Lead recommendation:** use stochastic adult receipts as the candidate for
funded equilibrium and recalibration. Do not retain guaranteed equal payments
as the preferred inheritance law, and do not select a receipt law from the
unrecalibrated loss. Keep observed entrant financial positions separate,
with their funding and creditor counterparties explicit. The empirical
age-only pooling, unsupported-age restriction, IID arrival assumption and
start-period timing must remain visible in the final contract.

The fixed-price cases intentionally do not clear housing or rebalance estates.
On the finer grid, housing excess demand is 2.952% under certainty and 2.691%
under risk; the lottery generates more estates than it pays by 6.446e-3 period
model units. Those feedbacks belong in the funded equilibrium stage.

The [nested-grid comparison](results/receipt_risk/finer_grid_v2/report/grid_effect_comparison.csv)
preserves the original entry point masses and their income conditionals
exactly, placing zero entry mass at inserted knots. Its first attempt stopped
at the fixed-entry-grid safeguard; the corrected and tested embedding changes
numerical indices only. Prices and the receipt scale remain common across
grids. The lottery effects relative to each grid's own no-receipt control are
stable on this refinement: the recent-parent gap changes by −1.611 versus
−1.615 percentage points; childlessness falls by 1.325 points on both grids;
fertility rises by 0.040 on both; wealth/earnings rises by 0.198 versus 0.196.
Ownership changes by +0.199 versus +0.105 points, so describe that effect as
small rather than attach importance to its magnitude.

**Separate numerical issue:** baseline levels are not established as converged.
Without receipts, refinement lowers the parent-ownership gap by 0.984 points,
lowers childlessness by 0.825 points and raises fertility by 0.026. This must
be addressed in the final calibration/convergence work. Stable receipt
differences on two grids do not certify converged calibration levels or policy
effects. The main integration task has this finding and owns its inclusion in
the final quantitative acceptance checks.

For integration, the stationary receipt hooks are not yet a dated policy
transition. If the estate scale changes over calendar time, the backward
continuation from date $t$ and the forward advance into date $t+1$ must both
use the **receiving date's** receipt law. Stationary age profiles coincide
across dates in this experiment, so it does not test that dated indexing.
The funded implementation must also derive the donor pool from the accepted
household mortality rule and explicitly state taxes, non-household recipients,
any credit losses at death and the funding of entrant financial positions.

## Resolution work — September 26, 2026

The author requires the open model issues to be scientifically evaluated and
resolved this week. The working integration deadline is Sunday, September 27.
The existing advisor checklist owns the consolidated issue list; this packet
owns estate recipients, their resource account, and compatibility with entrant
wealth. The completed comparison below remains valid experimental evidence.
It does not establish that its uniform, certain payment is the preferred
inheritance specification, nor does a larger fixed-parameter loss select the
appropriate receipt law.

### Evidence that disciplines the next decision

- [Sommer, Sullivan, and Verbrugge (2013), p. 857](https://www.kamilasommer.net/RentPriceRatio.pdf)
  explicitly send estate proceeds to government spending that does not affect
  household utility. Thus no household receipt can be an explicit resource
  closure. Their model does not by itself justify this project's warm-glow
  interpretation or its externally initialized entrant wealth. A government
  recipient cannot be silently relabeled as an inheritance received by children.
- [Feiveson and Sabelhaus (2018), Federal Reserve note](https://www.federalreserve.gov/econres/notes/feds-notes/how-does-intergenerational-wealth-transmission-affect-wealth-concentration-20180601.html)
  supplies inheritance probabilities and conditional amounts by age and usual-
  income group, with [accessible Figure 3 data](https://www.federalreserve.gov/econres/notes/feds-notes/how-does-intergenerational-wealth-transmission-affect-wealth-concentration-accessible-20180601.htm).
  These are pooled SCF 1995–2016 estimates, not a new 2007 estimate. Receipt
  probabilities refer to the previous three years; the income groups are
  within-age usual-income ranks. They are not automatically model labor-income
  ranks. The table covers ages 25–80 and does not supply sampling uncertainty.
- The existing public SCF 2007 codebook, inheritance section lines 27644–27739,
  excludes deceased-spouse receipts, combines life insurance and inherited
  trusts with inheritance in the public type code, and rounds receipt years
  to the nearest five years. A precise three- or four-year receipt hazard
  cannot be reconstructed from those rounded years. Donor-side estate-flow
  estimates in the existing 2007 packet do not identify recipient incidence.

### Finite experimental design and implementation

Use one pinned household specification for a no-receipt control and an initial
fixed-price, fixed-funding-scale pair: a deterministic conditional-mean
transfer and a lottery with the **same age-conditional mean**. This
pair isolates receipt risk. Subsequently solving separately funded equilibria
adds endogenous price and estate-pool feedback and is a different comparison.
The old equal payment to ages46–62 is retained as a separate comparison.
A fixed-price evaluation is an attribution check within this design, not a
substitute for resolving the receipt rule. Do not vary estate valuation,
preferences, earnings, initial
wealth, tax, or targets between the paired new cases. Any later integration
of independently accepted changes needs a new common control and full change
disclosure. The four-utility overnight source remains untouched.

For a published three-year probability $p_3$ and conditional total receipt
$a_3$, the proposed diagnostic period mapping is
\[
p_4=1-(1-p_3)^{4/3},\qquad
\mu_4=\frac43p_3a_3,\qquad
a_4=\mu_4/p_4.
\]
It assumes a constant arrival intensity within the age/income cell and
preserves annual expected receipts. It is not an estimated four-year law.
Representing all positive receipts by $a_4$ also omits variation in size
conditional on receipt. The published accessible amount header has an unclear
unit label; use only its relative amount profile until the physical units are
independently verified. A common source-unit conversion cancels when the
profile is scaled to the model estate pool.

For pre-receipt household mass $g(x)$ and conditional expected relative
receipt $\mu(x)$, the funded scale is
\[
\lambda=\frac{D}{\sum_x g(x)\mu(x)},\qquad
T^{\rm mean}(x)=\lambda\mu(x),\qquad
T^{\rm lottery}(x)=\begin{cases}
\lambda a_4(x)&\text{with probability }p_4(x),\\
0&\text{otherwise.}
\end{cases}
\]
Here $D$ is the estate pool available to the modeled recipients after the
explicitly specified liquidation costs, taxes and external recipients. The
funding identity must be recomputed in equilibrium. For a pure fixed-price,
fixed-transfer risk comparison, hold the same $\lambda$ in both cases and
report the resulting funding residual; do not call that mechanical comparison
a balanced equilibrium. In separate equilibria, changed donor behavior can
change $D$ and hence $\lambda$; this feedback must be reported separately.

Keep positive estates and unpaid liabilities separate. If the signed net
estate is $e=b'+(1-\psi)qh$, define
\[
D^+=\sum_x d(x)g(x)\max(e(x),0),\qquad
L=\sum_x d(x)g(x)\max(-e(x),0),\qquad
D^{\rm signed}=D^+-L.
\]
The positive pool funds household or external recipients and any estate tax.
$L$ instead identifies liabilities left at exit and requires a creditor/recovery
rule; do not automatically deduct $L$ from unrelated positive estates or call
it government expenditure. Selling costs are already inside $e$ and must not
be deducted a second time. If a modeled exit is the death of the entire
household decision unit, transfers to a surviving spouse within that unit are
not an additional cross-household receipt. Reconcile that interpretation with
the mortality mapping being reviewed by its existing owner.

The three-case fixed-price driver now makes these four experimental choices
explicit. It uses the exact published age nodes 26,30,...,78 and restricts
receipts to zero at 18,22,82; the unsupported-age restriction is not an
empirical finding. It pools published usual-income groups with weights
0.5/0.4/0.1, avoiding an unsupported mapping to model labor-income states,
especially during retirement. At each age, the zero/positive approximation
preserves the pooled probability and expected amount. Receipts are liquid
wealth at the start of the receiving period, before interest and choices:
$b\mapsto b+X$. They are not also added to income. This timing differs from
the old uniform income-transfer diagnostic and is held fixed across the new
mean/lottery pair. Receipts are independent across periods; a one-lifetime-
receipt restriction would be a separate economic assumption.

The implemented backward operator is $QV$ and the forward operator is
$Q^\top g$, with identical positive wealth-grid interpolation weights. It adds
no persistent state. The forward pass rejects any occupied wealth clipping;
it records receipt flows and verifies the financial-wealth and mass
identities. The same operator enters the native stationary distribution and
cohort transition used by the calendar and event observers. Exact control
array/target reproduction and the original calendar, budget, fiscal and
probability gates must pass. Housing and estate funding residuals are
reported, because prices and the estate-profile scale are fixed in this
attribution test. They are not treated as cleared markets or a funded GE.

### Decision and verification requirements

1. Reconcile the donor estate definition, surviving-spouse transfers,
   non-household recipients, debt written off at death, and external entry
   endowments. An external resource account is permitted only when explicitly
   stated and used consistently in policy comparisons. Do not automatically
   subtract later inheritance from observed entry wealth.
2. Preserve the exact no-receipt control. Verify zero-probability and
   probability-one limits; mass and wealth conservation at every receipt
   transition; identical backward/forward shock weights; purchase and debt
   feasibility; payroll/income-observer exclusion; and the original value,
   fiscal and market gates. Inspect numerical sensitivity where receipt
   interpolation crosses down-payment thresholds.
3. Report the full target-fit and parameter tables plus the stable diagnostic
   figures. The focal economic diagnostics are ownership and wealth before
   receipt, ownership after receipt, the recent-parent ownership **gap**
   (a descriptive association), childlessness, fertility timing, house prices,
   and wealth dispersion. Do not select on the scalar loss alone.
4. Compare recalibrated specifications using one common target contract and
   identifying restrictions before testing the central policy claim. A
   deadline does not convert a failed gate, underidentification, or an
   unsupported receipt mapping into evidence. Resolution requires a backed
   specification choice or an explicit revision of the quantitative claim,
   rather than another undated diagnostic proposal.

The source extraction and period-mapping builder is implemented in
`code/data/scf/build_inheritance_receipt_profile.py`. Torch job 18569789 passed
all eight tests and retained all 168 published rows. The largest annual-mean
mapping identity residual was 1.776e-15. The HTML source hash and the published
and mapped tables are in `results/recipient_evidence/v1/`.
`code/model/tools/audit_e5f_estate_resource_account.py` audits the retained
signed estates, sale costs, liabilities at exit and entrant financial positions
without a new household solve. Job 18569810 passed the signed liquidation-cost
and estate-observer identities. Negative net estates are exactly zero in all
three retained cases; mean entrant wealth remains 0.186520 annual-earnings
units in each. The full by-age account is in
`results/recipient_evidence/accounts_v1/resource_account.json`. This is a
household-boundary account, not a national resource-constraint certification.

The wealth-jump primitive passed six tests in Torch job 18570774. The new
adapter is `code/model/tools/e5f_estate_receipt_risk_adapter.py`; the bounded
driver is `code/model/tools/run_e5f_estate_receipt_risk.py`, submitted with
`code/cluster/submit_e5f_estate_receipt_risk.sh`. It permits three fixed-price
household solves, has a 25-minute execution cap, writes a heartbeat every
55 seconds and complete-case/best-fit receipts, and produces the unchanged
17-plot packet plus full fit and parameter tables for each case. The exact
case-loop fixture and adapter tests run before household execution. Reference
and source identities are pinned in its contract. No new recipient law,
entry law or production target has been adopted by this experimental work.

## Status — September 26, 2026

Author requested the three-way test discussed in the bequest-accounting task.
Torch access is restored. Native source/checkpoint preflight passed in job
18563842 (`results/v2/preflight`), and all nine focused fixture tests passed in
18564037 (`results/v2/smoke`). Comparison job 18564157 completed all three cases
in 63 minutes, using 36 equilibria (34 recipient iterations). The control
matched all saved solution arrays and targets exactly. Every case passed the
unchanged numerical gates with zero budget-excess mass; estate receipts balance
net deaths to 5.839e-11 model units. Both economic changes remain unadopted.
The first preflight (18563806) stopped before
model import because Slurm relocated its script; the launcher now preserves the
submission project root explicitly. Model work, imports, tests, source hashing
and rendering stay on Torch. Existing reference results are untouched.

## Completed comparison

The [64-page report](results/v2/final_report/estate_receiver_comparison.pdf)
contains all 13 target rows with weights, gaps and loss contributions, all 26
parameter/restriction rows per case, and the stable 17 figures for each case.
All 315 numeric target/parameter cells were independently extracted from the
PDF and matched to the CSVs. All 64 pages were reviewed as contact sheets, with
the tables and flagged figure pages also checked at full size; exact coverage
and the remaining readable inset-legend placement are recorded in the sidecar. The first
render exposed inherited crowded income legends; job18567191 corrected only
their formatting from saved solutions, with no equilibrium solves. All 24
modified legends fit inside their axes and preserve plotted data. Original
figures and the preliminary report are retained separately.

| Moment | Target | Control | Net valuation | Net valuation + receipts |
|---|---:|---:|---:|---:|
| Inherited weighted loss | — | 280.411 | 281.190 | 506.158 |
| Fertility measure | 2.100 | 2.100 | 2.100 | 2.155 |
| Childlessness (%) | 19.828 | 18.317 | 18.327 | 16.420 |
| Ownership ages 30–55 (%) | 67.626 | 75.676 | 75.694 | 72.216 |
| Wealth / earnings | 6.927 | 6.045 | 6.077 | 5.923 |
| Mean rooms | 5.608 | 6.544 | 6.547 | 6.630 |
| Mean first-birth age | 25.976 | 26.342 | 26.345 | 26.413 |

Complete machine-readable tables:

| Case | All targets, gaps, weights and loss contributions | All estimates, bounds and restrictions |
|---|---|---|
| Control | [Target fit](results/v2/run/control/target_fit.csv) | [Parameters](results/v2/run/control/parameters.csv) |
| Net valuation | [Target fit](results/v2/run/net_valuation/target_fit.csv) | [Parameters](results/v2/run/net_valuation/parameters.csv) |
| Net valuation + receipts | [Target fit](results/v2/run/net_valuation_transfer/target_fit.csv) | [Parameters](results/v2/run/net_valuation_transfer/parameters.csv) |

Net valuation alone has small effects. Relative to that case, funded receipts
raise the fertility measure by 0.056, raise house prices 2.143%, and lower
ownership ages 30–55 by 3.478 percentage points. The balanced transfer is 0.385
per eligible household per four-year period, or 0.096 per year, in units of
mean annual gross working earnings. The native birth-entry replacement gap is
0.026; these fixed-preference treatments are not closed demographic equilibria.

Effects precede actual receipt: at the age-42 node, mean net financial wealth
falls from 0.410 to -0.030 between valuation alone and receipts, while eligibility
starts at 46. Those lifecycle CSVs are retained alongside the target tables.
This documents changed choices before inheritance; the test does not separate
anticipated-income effects from housing-price feedback. Recalibration and policy
transitions are needed before claiming a policy effect or adopting the change.

For the original young-wealth concern, dividing the unchanged **control** estate
pool by entrant flow gives a mechanical entry grant of 2.003 mean annual earnings
using gross estates, or 1.770 net of selling costs. This is a one-time wealth
stock; do not annualize it by dividing by four. It is not a solved age-18
redistribution counterfactual: an endogenous estate pool and prices would change.

The larger inherited loss is not a model-selection decision. Existing measurement
approximations include the first-birth observer versus the empirical event-study
estimator, the recent-parent residence proxy, and all-positive estates versus
child-directed empirical estates. The unchanged experimental tax, mortality and
entry law also condition this comparison.

## Reference and economic changes

Reference: September 25 paired PAYGO experiment, `oasi_087510`, actual selected
case `worker_05/point_04`. Remote parent is
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/nightpair_20260925_v1`.
The control retains that experiment's tax rate 0.08751017424959717, all selected
structural parameters, fertility preference, B15 earnings, entrant joint wealth
distribution, wealth grid, housing supply, target rows/weights and numerical
gates. It does not implement the newly adopted pension ratio/tax rule and is
not an adopted calibration.

| Case | Bequest utility estate | Receipt by living households |
|---|---|---|
| control | `max(b' + q*h, 0)` | none |
| net_valuation | `max(b' + (1-psi)*q*h, 0)` | none |
| net_valuation_transfer | same net estate | equal anticipated period transfer to age cells starting in [45,65] |

Both economic changes are **experimental, not adopted**. The existing four-year
grid pays age-start cells 46, 50, 54, 58 and 62; this is not a precisely
age-integrated recipient profile. No family matching or random inheritance draw
is added. The existing current-income purchase rule includes the receipt in
`Y/R`, with no second credit to transaction wealth or next-period assets.

The receiving case solves `T*M = D`, where `M` is recipient household mass and
`D` is the period death-weighted positive estate flow net of selling costs.
Estates use post-saving wealth and the native post-transaction current measure;
terminal death has probability one. Entry wealth is not funded or changed by
this pool. There is no transfer at age 18.

The empirical bequest target remains the original **gross** observer in all
three cases. Net flows and receipts are supplemental accounting. Receipts do
not enter measured labor earnings. Thus a target-definition change cannot
masquerade as an economic response. The existing all-positive-estates versus
child-directed empirical-target mismatch remains disclosed and unchanged.

The fixed-preference comparison keeps normalized entrant/age composition and
reports the birth-based adult-entry replacement residual. Treatment cases are
not claimed to be closed-population stationary equilibria. It does not impose
fertility 2.1 by changing preferences. Policy transitions and any subsequent
recalibration remain a separate stage after this diagnostic.

### Joint entry/estate interpretation

The coordinated entry review is
`../native_financing_diagnostic_20260919/specification_followup/target_review_v1/entry_wealth_rationale_review/report.md`.
Observed entry wealth is a stock containing prior saving and any prior help;
later inheritances are subsequent flows. They are not automatically the same
transfer. Double counting occurs if a simulated pre-entry grant is added to a
stock already containing that grant, or if the same estate pool is spent twice.

The present test balances only $D_t=M_t T_t$. Entrant endowments remain a
separate external injection $N_t\mathbb{E}[b_0]$, where $N_t$ counts entering
households. A closed family-transfer account would need to identify the
pre-entry transfer component and its donors/funding. The common-scale candidate
and its empirical income dependence are unadopted and are not used here.

A transition must recompute dated deaths and recipient mass, and household
decisions must anticipate the corresponding transfers. A fixed per-entrant
wealth law still produces changing aggregate entrant resources when cohort
entry changes. The 16/20 calendar queue approximates mean 18-year maturation;
the household age grid in this experiment still begins at 18.

## Implementation and verification

- `code/model/tools/e5f_estate_receiver_adapter.py`: isolated in-memory income
  and Bellman hooks, supplemental death accounting and bounded transfer loop.
- `code/model/tools/run_e5f_estate_receiver_probe.py`: pins the parent lock,
  source, targets and selected checkpoint; reuses its verified purchase runtime,
  observers, fiscal/market/household/purchase/mass/value gates and 17 figures.
- `code/cluster/submit_e5f_estate_receiver_probe.sh`: separate Torch stages.
- Two `test_*estate_receiver*.py` files test period/annual accounting, death
  timing, income exclusion, funded receipts, control nesting, loop/cap failures.
- `code/model/tools/build_e5f_estate_receiver_report.py` assembles a completed
  three-case packet into a PDF with all tables and 51 retained figures; it
  checks numeric columns and renders every page on Torch.
- `code/model/tools/format_e5f_estate_receiver_figures.py` re-exports the same
  figures from saved checkpoints with compact income legends, verifies that
  curve data do not change during formatting, and checks legend containment.

Source review corrected the worker draft's ancestor-versus-selected preference
comparison, remote paths, solve counter, unconditional transfer balance gate,
and confusion of mass conservation with demographic replacement. Source
preflight and all nine focused tests passed on Torch. The native control
reproduced the selected arrays and every target model value exactly. No tolerance
was widened, and the comparison remained within its contracted budget.

Budget: 20-minute preflight allocation, 20-minute smoke allocation (tests capped
at 10 minutes), then a three-hour run deadline with 42 native solves maximum:
one control, one valuation case, at most 40 receiver iterations. The reference
chosen equilibrium took 99.811 seconds, implying about 70 minutes at that speed
for 42 solves; nonlinear convergence, audits and graph export add uncertainty.
The earlier different-model receiver test took 54 minutes. No search or repeated
full-case battery is planned. Every iteration writes a receipt; a 55-second
heartbeat, latest completed case and descriptive best-so-far remain readable.
Nonconvergence, baseline mismatch or a numerical gate stops the run.

Completed cases retain checkpoints, full target/parameter tables, estate
accounts and all figures on Torch. Compact receipts, tables and the final PDF
are also retained locally. Large model checkpoints remain on Torch.

## Reproduce the report on Torch

The numerical stage is already complete; report reproduction needs no new
equilibrium solve. Inside an allocation in the remote task tree, with the
Anaconda 2025.06 module and `PYTHONPATH` including `report_dependencies` and
`/scratch/td2248/commute_pdf_qa_deps`, use a **new** output filename:

```bash
python3 code/model/tools/build_e5f_estate_receiver_report.py \
  --run-root results/v2/run --figure-root results/v2/figure_formatting \
  --output results/v2/NEW_REPORT/estate_receiver_comparison.pdf
```

Report dependencies installed only in this task tree are ReportLab 5.0.1 and
pypdf 6.19.0; the existing Torch PDF runtime provides PyMuPDF and Pillow. The
recorded report job is `results/v2/final_report_job.sh` (18567191). Numerical
results and their fingerprints are unchanged by the display-only export.

## Resume on Torch

Place only the four new Python files and
launcher in a **separate** task tree, preserving their `code/model/tools/` and
`code/cluster/` paths. Suggested remote task root:
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/estate_receiver_probe_20260925_v1`.
Never copy them into or modify the frozen `nightpair_20260925_v1/source` tree.

The launched comparison uses `results/v2/`; those directories must not be
reused. For a separately justified future run, choose a fresh version below.
From the separate task root, the stage sequence is:

```bash
export E5F_ESTATE_PROBE_STAGE=preflight
export E5F_ESTATE_PROBE_OUTPUT="$PWD/results/NEW_VERSION/preflight"
bash code/cluster/submit_e5f_estate_receiver_probe.sh --submit
```

Inspect completion, then:

```bash
export E5F_ESTATE_PROBE_PREFLIGHT_RECEIPT="$PWD/results/NEW_VERSION/preflight/preflight.json"
export E5F_ESTATE_PROBE_STAGE=smoke
export E5F_ESTATE_PROBE_OUTPUT="$PWD/results/NEW_VERSION/smoke"
bash code/cluster/submit_e5f_estate_receiver_probe.sh --submit
```

The smoke runs the exact driver loop and actual adapter fixed point with
deterministic fixtures; no full native equilibrium is solved in that stage.
After its success, the first production case is the exact native baseline
replay, and treatments may proceed only after it passes:

```bash
export E5F_ESTATE_PROBE_SMOKE_RECEIPT="$PWD/results/NEW_VERSION/smoke/complete.json"
export E5F_ESTATE_PROBE_STAGE=run
export E5F_ESTATE_PROBE_OUTPUT="$PWD/results/NEW_VERSION/run"
bash code/cluster/submit_e5f_estate_receiver_probe.sh --submit
```

Receipts pin the complete parent source/target/checkpoint and the new diagnostic
files. A source change invalidates smoke/preflight. Do not retry or enlarge the
budget without a changed, stated diagnosis. The launcher does not poll or create
a monitor; active job status is recorded above.
