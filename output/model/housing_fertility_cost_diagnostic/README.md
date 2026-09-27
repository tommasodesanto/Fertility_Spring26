# Housing-cost sensitivity of fertility

Completed experimental comparison, September 26, 2026. No specification is adopted by this diagnostic. The subsequently chosen estate-funded entrant assets and residual sink are **not implemented here**: estate rules and entry remain frozen at the authenticated reference.

The [PDF](../../pdf/housing_fertility_cost_diagnostic.pdf) has 167 pages: ten concise comparison/method pages, supplemental mechanism plots, complete target and parameter tables, and the unchanged 17 diagnostic panels for each of 12 reported solutions. The complete current compact packet is **`combined_001/`**. Files directly in this parent directory are the retained earlier partial packet; do not combine them with current summaries. Full native checkpoints and grids remain on Torch.

## Main comparison

Both house asset prices and rents rise 10%, with the user-cost mapping unchanged. Earnings, pensions, payroll taxes, credit primitives, estate rules, entry, shock scales and all other preferences stay fixed. The share utility reference-rent normalization remains fixed. These are permanent fixed-price partial equilibria, not market-clearing GE or transition paths. Current owners also receive asset revaluation.

| Benefit treatment | Housing specification | Baseline fertility | Fertility after +10% prices | Fertility change | Birth-age change |
|---|---|---:|---:|---:|---:|
| Fixed | Floor | 2.041 | 1.940 | -4.976% | +0.334 years |
| Fixed | Share loading 0.100 | 2.618 | 2.580 | -1.464% | +0.174 years |
| Fixed | Share loading 0.200 | 2.590 | 2.528 | -2.368% | +0.278 years |
| Matched | Floor | 2.101 | 2.000 | -4.808% | +0.339 years |
| Matched | Share loading 0.100 | 2.101 | 2.038 | -2.955% | +0.193 years |
| Matched | Share loading 0.200 | 2.100 | 2.006 | -4.498% | +0.297 years |

The larger share loading produces meaningful fertility price sensitivity, close to the floor in this experiment. This weakens a claim that a floor is necessary for the mechanism; it does not select a calibrated winner. The floor is not uniformly stronger across comparisons: in the fixed-benefit common renter states the 0.200 loading produces a larger probability decline than the floor. The unchanged first-birth utility cost (0.506777), first-birth shock scale (0.209276), and later-birth scale (0.482069) remain important conditioning assumptions.

Child benefits are \(b m^{0.86}\), with \(m\) children at home. The fixed-benefit trials use \(b=0.1339982278467087\). Matching changes only positive \(b\), within [0.001, 0.400], to fertility 2.1 within 0.002, and holds that value fixed under the price shock. Matched values are 0.143611328125 (floor), 0.0574990234375 (0.100 loading), and 0.0645126953125 (0.200 loading), all interior. The CRRA coefficient is \(\psi=0.86b\), distinct from the code field `psi_child`, which stores \(b\).

This is one-dimensional diagnostic normalization, not full SMM. The inherited target system remains frozen, including the historical ACS quantity target and first-birth room target 1.465, measured against an unmatched stationary proxy. Matched baseline room responses are 0.785, 0.989 and 1.527. Complete descriptive fit tables, including every target, model value, gap, weight and contribution, are in each `combined_001/<case>/target_fit.csv` and the PDF. Every parameter, inherited restriction, matching bound and near-bound flag is in the corresponding `parameters.csv`.

## Household mechanisms

The following responses use identical control weights for childless households aged 18–34 and optimize fertility and tenure. Probability changes are percentage points. Housing and nonhousing are model units.

| Matched specification | Initial tenure | Birth probability change | Owner probability change | Housing change | Nonhousing change |
|---|---|---:|---:|---:|---:|
| Floor | Renter | -2.729 | -2.647 | -0.380 | -0.007 |
| Floor | Owner | -0.307 | +0.861 | -0.148 | +0.027 |
| Share 0.100 | Renter | -1.720 | -2.359 | -0.374 | -0.002 |
| Share 0.100 | Owner | -0.380 | +0.881 | -0.159 | +0.034 |
| Share 0.200 | Renter | -2.639 | -2.814 | -0.405 | +0.003 |
| Share 0.200 | Owner | -0.822 | +0.867 | -0.141 | +0.028 |

Households adjust housing and births jointly. Nearly unchanged mean nonhousing spending does not mean births are unaffected. These endpoint comparisons do not decompose lifetime lost births from transition postponement.

The attempt-now-minus-wait value gap excludes the current fertility taste shock and retains future shock-inclusive continuation values. Under the verified simple binary logit it is recovered as \(\log p_{try}-\log p_{wait}\), in shock-scale units, only when both probabilities are positive. Endpoints and unavailable choices are separately classified. A negative gap does not imply negative direct utility from children or preference for lifetime childlessness.

At matched baselines, positive-gap first/later birth shares are respectively 55.82%/27.58% (floor), 13.16%/5.32% (0.100), and 12.04%/2.78% (0.200). These use **case-specific stationary weights**. First-birth risk median normalized gaps are -0.980, -0.716 and -0.739; lower positive-gap birth shares therefore do not establish that one specification uniformly makes children less desirable. Benefits and stationary distributions differ.

Lead with realized common-state birth responses, not average normalized-gap magnitudes: tiny probabilities can dominate log gaps, and finite-in-both summaries exclude endpoint mass. For example, the floor's common first-birth gap comparison omits about 8.06% of common risk mass. The illustrative liquid-wealth-at-most-one subset is not an empirical bottom quantile. Supplemental grid panels use ages 22, 26 and 30 and fixed control income deciles. Conditional renter housing for “child 1” is the current childless (n,m)=(0,0) family-state policy, not post-birth housing or an event response; “child 2” uses (1,1). Realized allocations after tenure choice are separately reported.

## Verification and provenance

All 12 reported points pass unchanged probability, occupied-value monotonicity, household budget, feasibility, transition accounting and distribution checks. Stressed common-state budgets, unavailable-choice mass, and feasibility are separately audited. All 32 parameter rows are checked against saved actual parameters/derived values or the two explicitly named inherited external restrictions (pension ratio and adult-entry birth conversion). All 13 target rows have checked gap/loss arithmetic. Every nonprivate parameter remains identical within each baseline/shock pair. Actual fiscal residuals and PE market/replacement residuals are retained; market clearing is not imposed. Intermediate bracket points receive native solve/reconstruction checks, not the full final-case audit.

- Current evidence: `combined_001/final_verification.json`, `finalization_checks.json`, `common_state_gate_checks.json`, `price_sensitivity.csv`, `common_state_changes.csv`, `common_grid_gap_changes.csv`, `complete.json`, `solve_history.csv`, and `pdf_receipt.json`.
- Original Torch directory: `/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/housing_fertility_cost_diagnostic_20260926/run_001`, preserved unchanged. Original production job 18605876 stopped at the numerical time limit after 32:40; scientific gates did not fail. Saved-result audit/report jobs were 18609345 and 18609458.
- Authorized continuation: sibling `continuation_001`, combined packet `combined_001`. Smoke 18611752 passed; first production setup attempt 18611768 failed in 16 seconds because smoke and run reused an existing setup directory, before any solve. Stage/job-specific setup paths fixed that issue.
- Continuation job **18611784 completed in 5:41**, exit 0, under a 20-minute cap. Four new native evaluations: three midpoint refinements and one price shock. The original 12-evaluation matching cap sufficed; the authorized maximum of 16 was unused. Every new packet was saved immediately; no completed solve was repeated.
- Source checkpoint pins, continuation smoke/manifest, exact launch script and logs are in `continuation_001/`. Original native numerical functions are AST-identical to the retained executed driver; report changes are separate. All model work, hashing, plots and PDF generation ran on Torch with one numerical thread.
- PDF SHA256: `460f3b3f8b545b3938f34d709a8d4136995afa1d9ad7bbc8fd74be5075e44644`. All 167 pages have no out-of-page text. All eleven contact sheets were visually inspected; first comparison pages also inspected at full size. No automatic preview was opened.

## Reproduction

Owned sources are `code/model/tools/run_e5f_housing_fertility_cost_diagnostic.py`, `e5f_housing_fertility_cost_audit.py`, and `e5f_housing_fertility_cost_resume.py`. Native science is loaded from the authenticated comparison snapshot, never current solver defaults.

Use Torch module `anaconda3/2025.06`, one numerical thread, `PYTHONDONTWRITEBYTECODE=1`, `MPLBACKEND=Agg`, and the exact `PYTHONPATH` in `continuation_001/launch.sh`. The resume helper authenticates all pins, bracket history, no-duplicate points and isolated outputs before solving. Do not rerun a completed continuation or overwrite retained output folders.

The final PDF alone regenerates without any new solves using:

```sh
python tools/run_e5f_housing_fertility_cost_diagnostic.py --stage render --output combined_001
```

Run this from the remote diagnostic root inside a Torch Slurm allocation with the recorded environment. Full saved-result auditing uses `tools/e5f_housing_fertility_cost_audit.py --output <fresh isolated copy>`; its runtime setup directory must be fresh. Reproduction preserves the original case checkpoints and graph sets, using copied compact tables in the combined directory and links to immutable heavy artifacts.
