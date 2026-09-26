# Estate recipient diagnostic

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
