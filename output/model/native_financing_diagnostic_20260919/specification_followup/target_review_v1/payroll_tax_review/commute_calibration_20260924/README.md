# Experimental commute calibration — September 24

The detached Torch chain completed as smoke `18471839`, eight-worker array
`18471840`, and selected-export job `18471841`. The remote run
root is `/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/commute_calibration_20260924_v1/results/run_001`.
The shared deadline, measured from first smoke start and covering search and
export, was **2026-09-24 20:56:08 EDT**. Smoke passed its full loop and
17 standard figures. Production allowed at most five proposals per worker, so
the chain planned at most 41 objectives including smoke. It scored 23 cases;
one further worker case timed out inside a native solve, and 17 were unrun.
Seven worker jobs completed; worker 6 exited 124 at the shared cutoff. The
selected export completed in 26 seconds, before the deadline. No numerical
retry, target edit, or rescore occurred.

The objective contract is `objective.json` (SHA-256
`c28e5d620d463dc3a592c038ca2771b47bf37ed1a9cf5ba31f7dc793b58b61db`):
eight free structural coordinates, twelve positive-weight moments with every
previous working weight retained, and a separate completed-fertility
normalization at $2.1$. The selected B-floor checkpoint SHA-256 is
`83a28e46b36e2fbe30338d366611f3ec209f0c5a68309ee4ee9fa8523b66adee`.
The frozen native source inventory SHA-256 is
`76406fcc10206d9e30bcc29d4e18accdf7fffdd9219e50e6e11c1c3360f01336`.
The manifest pins the new driver, Slurm scripts, objective, selected plan,
receipts, and reviewed PAYGO comparison driver. `preflight_receipt.json` records
the Torch zero-solve checks. Source and numerical gates are otherwise unchanged.

Relative to the selected September 23 configuration, the frozen experiment applies
the author-accepted national housing targets, PSID wealth/earnings ratio
$6.92658379107299$, SCF bequest-flow/wealth ratio $0.007291023472616158$,
first-birth rooms response $1.465$ provisionally, externally fixed $\theta_1=
0.008193084126995582$, annual depreciation $0.01416143718381309$, and annual
property tax $0.010598360773872594$. The PAYGO payroll tax
$0.08751017424959717$ is **experimental, not adopted**. The floor utility,
B15 income process, inherited heterogeneous entry wealth, supply elasticity
$0.63$, other unchanged targets and weights, and stationary entry law remain.
The adopted 16/20 birth-entry queue has not been implemented, so this is not a
new demographic steady state. The model's stationary housing observer is not
the PSID panel estimator; its bequest observer counts positive estates across
all child states while the accepted SCF target is child-directed. These two
measurement mismatches remain explicit approximations; no model observer or
recipient rule was changed. Do not compare this objective's loss with old
target systems.

## Completed result and review

The **tax question** has two distinct pieces of evidence. The earlier
[fixed-parameter tax comparison](../tax_comparison_20260924/README.md) solves
the same selected parameter vector and utility scale at payroll rates 17.9%
and proposed 8.751%, balancing the pension and re-solving equilibrium in each
case. Under that controlled change, annual pension / working gross earnings
falls from 51.159% to 25.011%, wealth / earnings rises from 5.241 to 7.263,
ownership at ages 30–55 rises from 39.493% to 41.881%, and childlessness at
ages 40–44 falls from 16.909% to 12.538%. Its
[complete 13-row target comparison](../tax_comparison_20260924/full_target_comparison.csv) uses the older target
contract and does not re-normalize fertility. Tonight's lower-tax calibration
also changes empirical targets, externally fixes $\theta_1$, updates
depreciation and property tax, and re-normalizes fertility. It demonstrates a
feasible recalibration of that **bundle** under the proposed tax; no
high-tax refit under the same updated contract was run, so differences in fit
cannot be attributed to tax alone.

The selected case is remote `worker_03/point_02`, with weighted loss
**381.0485741441105**. The 13-row [full target fit](final_selected/target_fit.csv)
and 24-row [parameter and restriction table](final_selected/parameters.csv)
are the exact machine outputs; [case receipt](final_selected/receipt.json) and
[export receipt](final_selected/final_receipt.json) pin the objective, native
source, selected checkpoint, and 17 selected-case figures. The
[reviewed PDF](../../../../../../pdf/commute_calibration_review_20260924.pdf) contains
four corrected report/table pages and the unchanged 17 selected diagnostic
pages. Its [report-only builder](final_selected/make_reviewed_report.py) and
[review receipt](final_selected/report_review_receipt.json) verify all 13 target
and 24 parameter rows and pixel identity of the original 17 figure pages.
The original 21-page export remains untouched as
`final_selected/experimental_commute_calibration.pdf`.

The first-birth target choice is **pending author clarification**. The run's
$1.465$ target was our interpretation of “the new one”; the author has since
questioned that reading and recalled the earlier $1.025$ household estimate.
The run remains a frozen experimental diagnostic, not an adopted target
contract. Its first-birth model value is $0.9786043315988904$ and its
ownership-at-ages-30–55 value is $0.40575669234977646$ against target
$0.6762604168538028$. The largest loss contributions are ownership
$171.17644338031667$ and mean rooms $95.98362055094853$.

**Preflight/review error:** the actual frozen annual discount-factor search
upper bound was **0.9995**, inherited from the legacy contract, although the
requested cap was **0.99**. The objective and original PDF were not altered;
the latter rounds the bound to 1. All 23 scored candidate values were at most
$0.9762497404726131$, and the selected value is the same, so no observed case
crossed $0.99$. This does **not** make the launched search contract compliant.
The reviewed PDF and exact parameter CSV state the discrepancy.

The remote ledger records **144 stationary solve starts, 143 completions, and
one incomplete solve** in timed-out worker 6 point 03. Across all 23 scored
receipts, household budget-excess mass, end-mortgage-floor violation mass,
purchase-threshold violation mass, transaction-outside-grid mass, and
occupied negative value-step count are zero. Fiscal gates passed; maximum
scaled pension-budget residual is $3.382765664248042\times10^{-12}$ and
maximum housing-market residual is $1.9347972762613373\times10^{-5}$.
All scored receipts match the pinned source, checkpoint, and objective hashes,
and all eight files in the remote launch manifest still match their hashes.
The source `owner_rungs` figure concentrates owned demand at the 10-room
maximum; this is an observation for later diagnosis, not an established cause
of the ownership miss. Income-state figure legends remain crowded.
