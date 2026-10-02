# First-birth rooms income sensitivity (v1)

## Facts

- Reference estimator: byte-pinned canonical `sa_rooms_first_birth_v2.do`, arm `A2h` (SHA-256 `49ebdb6c4c780eaa6eb00be0c66c0f8523af85ed0f2c049c48d286f2281d8b8e`); the main target and reference are unchanged. Prepared v2 data are pinned to `9fd181c48e1951a7a0052836926260a89d056d666aad3e49e2ae9cca05f1ba00`.
- Design: current adults ages 18+, IW weights, biological first birth, confirmed-childless controls, pre-treatment reference-person/spouse treated status, ID/year fixed effects, age/education indicators, and two-year windows. Headline window is (+3/+4) relative to (-3/-2).
- Income: `INCFAMR` is FU total family income, PCE-adjusted to 2022 USD, for tax year `interview_year-1`. Common fits retain only observed, strictly positive income with `INCFAMR != 9999999`; the third fit adds `log_family_income = ln(INCFAMR)`. The top-code marker is excluded. Income is not imputed, winsorized, or adjusted for family size.
- Comparison: `baseline_full` is byte-identical A2h; `common_no_income` applies the income-validity restriction immediately before estimation; `common_income` uses the identical restriction and adds log income. Missing, nonpositive, and top-coded exclusions are separately recorded.

## Verified results

The headline is the change in rooms at event years $+3/+4$ relative to $-3/-2$.

| Specification | Rooms | Standard error | Fitted observations | Individuals |
|---|---:|---:|---:|---:|
| Original sample, no income | 1.465293 | 0.050270 | 117,853 | 9,310 |
| Income-observed sample, no income | 1.461421 | 0.050325 | 117,450 | 9,302 |
| Same income-observed sample, log real family income | 1.299762 | 0.047462 | 117,450 | 9,302 |

The income control lowers the common-sample coefficient by 0.161660 rooms
(11.062%). The sample restriction alone lowers it by 0.003872 rooms. These
comparisons are descriptive changes between specifications, not a causal
mediation decomposition. The broad post-birth housing increase remains in
both specifications. The long pre-birth bin is positive in both curves; this
sensitivity does not establish parallel trends or a causal direct effect.

The original candidate sample has 119,334 rows; 407 fail income validity:
0 missing, 395 nonpositive and 12 top-code markers. After the estimator's own
sample restrictions, fitted observations decline by 403. All 3,302 treated
individuals remain; controls decline from 6,008 to 6,000. All remaining treated
cohorts retain reference-window support.

The unchanged baseline matches every saved coefficient, standard error and
non-runtime fit-receipt field within $10^{-7}$, with exactly the original
fitted keys. Both common fits share digest
`4ef0659f23c52cabdf14f618773e6358cd98f0650907c9695b7fb2b8ed11cf78`.
Finite, symmetric positive-semidefinite covariance, coefficient/SE consistency,
headline consistency, native completion and zero unsupported cohorts all pass.
See [verification.json](verification.json), [comparison.csv](comparison.csv),
[full_window.csv](full_window.csv) and [the graph](full_window.png).

The full three-fit run completed in 21.6 minutes. The current driver differs
from the executed version only in post-run receipt validation, graph
labeling and generated CSV line endings; `executed_driver.py`, `executed_launcher.sh`, and `run_config.json`
preserve native estimation provenance. Source microdata and individual keys
are excluded from this results folder. No main target or paper specification
was replaced.

## Readiness

All three fits and their verification passed. The driver stages literal anchored copies, streams hashes, extracts private income with a 300-second bound and source-stat guard, and requires Stata markers, preserved cohort reference support, equal common-fit final-sample hashes, finite symmetric PSD covariance, and canonical baseline equality. It deletes private keys only after all three fits pass; `collect-verify` writes only aggregate coefficients, covariance, receipts, key hashes, support files, CSVs, PNG, and `verification.json` here.

## Reproduction

Use `extract-income`, then `stage` with explicit reference/sample/ado paths as needed, then the cluster wrapper with `full`; finally run `collect-verify` against the supplied canonical A2h output directory. The full-window graph sets the omitted −3/−2 baseline to zero and shows original, income-observed, and income-observed-plus-income curves with confidence bands. Do not add this sensitivity to the main target without an explicit research decision.

## Interpretation

Time-varying family income may respond to childbirth. Including it changes the
housing comparison and may condition on a consequence of the event. Treat this
as sensitivity evidence rather than an estimate of a causal direct effect. See
[Caetano et al., discussion of post-treatment controls](https://hsantanna.org/badcontrols/articles/bad-controls-conceptual.html)
and [Kleven, Landais and Søgaard on earnings after childbirth](https://www.nber.org/papers/w24219).
A fixed pre-birth income scalar is absorbed by person fixed effects; this does
not eliminate differences in housing trends associated with pre-birth income.

## Execution record

The native Torch synthetic loop passed in job `19002586` before full job
`19002614`. A launcher-path failure (`19002444`) occurred before estimation;
the first native synthetic attempt (`19002574`) exposed missing Mata-library
initialization. Both are preserved under the private task root and recorded
in `plan.json`. Each issue was repaired explicitly before a fresh smoke test;
there is no automatic retry loop.

### Exact pinned inputs

The full task root is
`/scratch/td2248/projects/Fertility_Spring26_rooms_income_sensitivity_v1_20261001`.
It contains the frozen reference estimator, prepared sample, income extract,
ado directory and staged estimator copies. On a fresh task root with these
inputs, stage with:

```bash
python3 code/data/psid_followup_mar2026/rooms_income_sensitivity.py stage "$task_root"   --reference-estimator "$task_root/reference_estimator.do"   --analysis-sample "$task_root/analysis_sample.dta"   --income-dta "$task_root/income.dta" --ado-dir "$task_root/ado"
```

Run the native toy phase on a separate synthetic task root first. The full
launcher takes `full "$task_root"`; set `ROOMS_PYTHON_BIN` to a Python with
NumPy. The collector additionally needs Matplotlib, the saved canonical
`output/sa_rooms_first_birth_v2/A2h` directory and the expected original fitted
key digest recorded in `plan.json`. Never collect private sample-key CSVs.
