# Python Model Codebase

This folder contains the active Python implementation of the current
discrete-time center-periphery fertility model. The former MATLAB code is now
archived for historical reference and parity checks.

Use the local virtualenv created for the port:

```bash
cd code/model
.venv/bin/python -m dt_cp_model.cli smoke
```

The virtualenv uses conda's NumPy and has `numba==0.58.1` installed for the
compiled golden-section savings kernels. If the environment has to be rebuilt:

```bash
/Users/tommasodesanto/miniconda3/bin/python -m venv --system-site-packages .venv
.venv/bin/python -m pip install 'numba==0.58.1'
```

Useful commands:

```bash
.venv/bin/python -m dt_cp_model.cli smoke --quiet
.venv/bin/python -m dt_cp_model.cli solve-theta --setup fast --theta x0 --max-iter-eq 120 --quiet
.venv/bin/python -m dt_cp_model.cli benchmark --setup fast --max-iter-eq 1 --force-full --quiet
.venv/bin/python -m dt_cp_model.cli objective --setup fast --max-iter-eq 1 --inv-iters 1 --quiet
.venv/bin/python tools/compare_benchmarks.py benchmarks/fast_ge_iter1.json benchmarks/matlab_fast_ge_iter1.json
.venv/bin/python tools/compare_benchmarks.py benchmarks/solve_theta_x0_fast_compiled_forward.json benchmarks/matlab_solve_theta_x0_fast.json
.venv/bin/python tools/check_population_closure.py
.venv/bin/python tools/export_room_pattern_diagnostics.py --case hR8=../../output/model/reduced_target_overnight_20260527/records/hR8_default_best.json,8.0 --case hR6=../../output/model/reduced_target_overnight_20260527/records/hR6_micro_best.json,6.0
.venv/bin/python tools/run_policy_counterfactuals_from_record.py --outdir ../../output/model/policy_counterfactuals_live_20260528_hR8_default
.venv/bin/python -m dt_cp_model.cli accounting-scale --setup fast --max-iter-eq 120 --quiet
.venv/bin/python -m dt_cp_model.cli scaled-equilibrium --setup fast --baseline-max-iter-eq 120 --max-iter-eq 80 --quiet
.venv/bin/python -m dt_cp_model.cli scaled-equilibrium --setup fast --baseline-max-iter-eq 1 --max-iter-eq 80 --outside-value -35.98609192692343 --force-full --quiet
```

## Fast Intergen One-Run Review

For the one-market intergenerational strand under
`intergen_housing_fertility/`, the default interactive inspection path should
stay small: one candidate, one solve or trusted solution cache, one quick plot
packet. Do not add `--run-policy-cases` unless the goal is explicitly a GE
counterfactual exercise.

From this directory:

```bash
PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/build_intergen_mechanics_packet.py \
  --source ../../output/model/intergen_current_review/source_candidate.json \
  --outdir ../../output/model/intergen_current_review/quick \
  --target-set candidate_replacement_roomgap_14moment_tfr192_v1 \
  --J 17 --Nb 60 --income-states 5 --n-house 5 --hR-max 6.0 \
  --max-iter-eq 10 --interp-method linear \
  --clean-outdir --no-csv \
  --skip-standard-diagnostics --skip-contact-sheet --quick-first-look-only
```

If a trusted `solution_cache.pkl` already exists, rebuild the quick plots
without re-solving:

```bash
PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/build_intergen_mechanics_packet.py \
  --source ../../output/model/intergen_current_review/source_candidate.json \
  --outdir ../../output/model/intergen_current_review/quick \
  --target-set candidate_replacement_roomgap_14moment_tfr192_v1 \
  --J 17 --Nb 60 --income-states 5 --n-house 5 --hR-max 6.0 \
  --max-iter-eq 10 --interp-method linear \
  --solution-cache ../../output/model/intergen_current_review/solution_cache.pkl \
  --refresh-plots-from-cache \
  --clean-outdir --no-csv \
  --skip-standard-diagnostics --skip-contact-sheet --quick-first-look-only
```

With `--no-csv`, the quick packet writes only the target table, moments,
solution summary, liquid- and total-wealth first-look policy figures, and
liquid- and total-wealth density figures. Use the full packet only for a slower
audit.

To audit the entry-wealth atom and owner-rung total-wealth jaggedness from the
current point-entry and five-node-entry caches:

```bash
PYTHONPATH=$PWD .venv/bin/python tools/audit_intergen_entry_atom.py \
  --outdir ../../output/model/intergen_current_review/atom_audit
```

To audit wealth/income denominator conventions and the price/rent scale before
changing wealth targets:

```bash
PYTHONPATH=$PWD .venv/bin/python tools/audit_intergen_wealth_units.py \
  --outdir ../../output/model/intergen_current_review/wealth_units_audit
```

The intended model boundary is:

```python
from dt_cp_model import solve_theta

sol, P, p_eq = solve_theta(theta, setup_mode="fast")
```

## Parent-Gated Bequest Calibration Exercise

The bounded bequest experiment profiles the normalized parent-gated luxury
warm glow, with and without an externally pinned post-retirement owner-LTV
taper. Production defaults remain unchanged. A local exact-loop smoke is:

```bash
cd code/model
PYTHONPATH=$PWD .venv/bin/python tools/run_intergen_bequest_calibration_exercise.py \
  --smoke \
  --outdir ../../output/model/intergen_bequest_calibration_exercise_smoke
```

The production-sized diagnostic uses the verified clean-frontier seed, the
active 15-moment objective, a 2-by-2 profile in `(theta0, theta1)` under each
LTV arm, and an optional bounded 11-coordinate polish. It reports the
late-life decumulation ratio as diagnostic-only because no empirical target
has yet been approved. Use
`code/cluster/submit_intergen_bequest_calibration_exercise.sh` for the real run.

The follow-up proper joint calibration uses
`tools/run_intergen_bequest_exit_chain.py`. Unlike the bounded diagnostic, it
re-estimates all 11 clean-frontier coordinates in every arm and adds a free
soft-zero `theta0` in the 12-parameter headline arms. The main array, primary
winner selector, nuisance-reoptimized profile, and fresh A3 Jacobian use:

```text
code/cluster/submit_intergen_bequest_exit_battery.sh
code/cluster/submit_intergen_bequest_exit_selector.sh
code/cluster/submit_intergen_bequest_exit_profile.sh
code/cluster/submit_intergen_bequest_exit_jacobian.sh
```

See `output/model/intergen_bequest_exit_battery_20260714/README.md` for the
complete identification contract, external variants, acceptance rules, job
IDs, and budgets.

## Clean Mortality Plus Bequest Test

The follow-up clean specification uses SSA post-retirement survival, no owner-
LTV taper, and a normalized child-blind warm glow. Arm `M2` in
`tools/run_intergen_bequest_exit_chain.py` re-estimates the 11 clean-frontier
parameters plus `theta0`; `theta_n=0`, `tenure_choice_kappa=0`, and
`theta1=0.25` are external restrictions. The completed mortality-only `M1`
winner is injected as the exact nested `theta0=0` seed.

Cluster and collection entry points:

```text
code/cluster/submit_intergen_mortality_bequest_recalibration.sh
code/model/tools/collect_intergen_mortality_bequest_recalibration.py
```

The recoverable run contract is
`output/model/intergen_mortality_bequest_recalibration_20260715/README.md`.
For the strict M1 identification and leave-one-moment-out audit, use
`tools/audit_intergen_bequest_exit_jacobian.py --arm M1 --winner-arm M1`; the
completed contract is
`output/model/intergen_mortality_identification_20260715/README.md`.

## Standard Internally Calibrated Bequest Block (M4)

Arm `M4` retains the clean SSA-survival specification and normalized
child-blind De Nardi warm glow but estimates both remaining bequest parameters:
the 11 clean-frontier coordinates plus `theta0` and `theta1` are disciplined by
14 moments. The child multiplier `theta_n=0` and deterministic tenure
`tenure_choice_kappa=0` remain external restrictions. The two late-life wealth
levels are the age-76--84 median total estate-to-income ratio and the age-65--75
reference-person median nonhousing-wealth-to-income ratio; estate dispersion
and the family-size estate gap are diagnostic-only.

Production entry points are:

```text
code/cluster/submit_intergen_standard_bequest_nested_reference.sh
code/cluster/submit_intergen_standard_bequest_recalibration.sh
code/cluster/submit_intergen_standard_bequest_recalibration_collector.sh
code/cluster/submit_intergen_bequest_exit_jacobian.sh
code/cluster/submit_intergen_standard_bequest_theta1_profile.sh
code/cluster/submit_intergen_standard_bequest_theta1_profile_collector.sh
```

The six-chain collector requires two bit-identical strict tight solves and
dominance over the exact `theta0=0` M1 seed under the M4 objective. The fresh
13-column Jacobian and five-cell conditional `theta1` profile are the
post-estimation identification gate; `theta1=0.25` is only one dispersed start.
The recoverable contract and results live under
`output/model/intergen_standard_bequest_recalibration_20260716/`.

Regenerate the exact M4 winner's full visual diagnostic packet with:

```bash
PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/build_intergen_mechanics_packet.py \
  --source ../../output/model/intergen_standard_bequest_recalibration_20260716/report/results.json \
  --outdir ../../output/model/intergen_standard_bequest_recalibration_20260716/diagnostic_packet \
  --target-set candidate_replacement_bequest_median_composition_v1 \
  --J 17 --Nb 120 --income-states 5 --n-house 5 \
  --max-iter-eq 40 --tol-eq 2.5e-5 \
  --m4-standard-bequest --clean-outdir
```

That exact solve writes the paper's two established quantitative figures under
`diagnostic_packet/classic_draft/`. Redraw only those two figures from the
trusted solved cache (about five seconds in the verified local smoke) with:

```bash
PYTHONPATH=$PWD .venv/bin/python tools/plot_intergen_draft_figures_from_cache.py \
  --solution-cache ../../output/model/intergen_standard_bequest_recalibration_20260716/diagnostic_packet/solution_cache.pkl \
  --lifecycle-out ../../output/model/intergen_standard_bequest_recalibration_20260716/diagnostic_packet/classic_draft/quant_lifecycle_equilibrium_repaired_nb120.png \
  --decision-rules-out ../../output/model/intergen_standard_bequest_recalibration_20260716/diagnostic_packet/classic_draft/quant_decision_rules_repaired_nb120.png
```

Redraw the full packet, including those figures, without re-solving with:

```bash
PYTHONPATH=$PWD .venv/bin/python tools/build_intergen_mechanics_packet.py \
  --outdir ../../output/model/intergen_standard_bequest_recalibration_20260716/diagnostic_packet \
  --target-set candidate_replacement_bequest_median_composition_v1 \
  --m4-standard-bequest --refresh-plots-from-cache
```

Regenerate the matched PSID/model late-life portfolio decomposition with:

```bash
PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/diagnose_intergen_bequest_distribution.py \
  --winner-json ../../output/model/intergen_standard_bequest_recalibration_20260716/report/results.json \
  --winner-arm M4 --arm M4 \
  --outdir ../../output/model/intergen_standard_bequest_recalibration_20260716/distribution_diagnostic \
  --quiet
```

## Intergen Entrant-Feasibility Diagnostic

The intergenerational model permits legacy unsecured debt to roll forward,
tapers that capacity between ages 42 and 62, and rejects any parameter vector
that places positive population mass on a Bellman-infeasible state. The
collateralized owner floor remains separate, and the cash down-payment test
continues to use `(1 - phi) * p * H`.

Regenerate the fixed-theta acceptance packet with:

```bash
NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  PYTHONPATH=$PWD .venv/bin/python \
  tools/run_intergen_feasibility_fix_diagnostic.py --overwrite
```

The driver writes both prescribed `Nb=120` records, the complete 15-moment
comparison, and an acceptance summary under
`output/model/feasibility_fix_diagnostic_<date>/`. A failed acceptance gate is
written to the packet and returned as a nonzero exit status.

The SMM objective should remain a wrapper around this same parameter-vector
application path, not a separate interpretation of `theta`.

## Population Closure

The live calibration uses the benchmark-normalized outside-option closure:
`P.population_closure = "outside_option_benchmark_normalized"`.

In the benchmark, the model solves for a normalized stationary distribution,
calibrates the outside value to a target city-entry probability \(q^E\), and
sets the outside-born potential entrant mass mechanically so that benchmark
scale is \(S=1\):
\[
S E_0(p)=q^E(p)\left[M+S B_0(p)\right],
\qquad
M=\frac{E_0(p)}{q^E(p)}-B_0(p).
\]

Counterfactuals should hold the benchmark \(M\), \(\bar W^E\), and entrant
taste scale fixed, then let the implied stationary scale move with prices and
fertility. The older `normalized`, `renewal_valve`, and
`accounting_scale_prices` closures remain in code for diagnostics and
comparison, but they are not the current benchmark calibration closure.

The policy counterfactual runner follows that boundary: it first re-solves the
record under `outside_option_benchmark_normalized`, recovers the benchmark
outside objects, and then solves policy cases with
`population_closure="accounting_scale_prices"` while holding those outside
objects fixed. By default it runs the two current housing-supply cases:
all-location `H0 +10%` and center `H0 +10%`. Add repeated `--case ...`
arguments to run credit, tax, or transfer cases from the same driver.

Run the regression guard after changes touching entry, fertility, housing
demand, or equilibrium iteration:

```bash
.venv/bin/python tools/check_population_closure.py
```

See `PLAN.md` for the full implementation and optimization plan.
