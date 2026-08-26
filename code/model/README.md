# Python Model Codebase

This folder contains the project's active Python model implementations. The
former MATLAB code is archived for historical reference and parity checks.

## Current paper architecture (August 18)

### Integrated demographic-transition repair (August 26)

The certified baseline still uses the historical `births / 2.1` age-20
household-renewal queue and has not been overwritten. An isolated replacement
person-cohort layer under `demographic_transition/` is now integrated with an
experimental perfect-foresight household and housing-market solver. The new law separately advances births,
age-specific survival, age-sex net migration, nonhead-to-head household
formation, head dissolution, and household-head mass.

The verified path packet is
`../../output/model/e5f_coherent_person_cohort_path_20260826a/`. It starts from
the exact Vintage 2025 age-sex population, reproduces the cohort-anchored Census
main path cell by cell, and writes the proposed four-year model interface. The
terminal diagnostic is
`../../output/model/e5f_coherent_terminal_demographic_diagnostic_20260826a/`.
It shows that the old queue-based terminal birth rates imply a corrected-law
renewal ratio of about `1.0625`; hence the imposed terminal preference
`psi_child=0.27460458049447606` cannot remain the terminal condition after the
queue is removed. The household policy, price, transfer, demographic scale,
and terminal preference must be solved jointly under the new law.

The coupled one-period regression and terminal-mapping tests pass under
`tools/run_e5f_perfect_foresight_person_demography.py`. Conditional terminal
price/rebate/demographic roots are solved by
`tools/solve_e5f_person_demography_terminal_root.py`; the verified 1%/2%
comparison is
`../../output/model/e5f_person_demography_terminal_comparison_20260826a/`.
The roots imply `184.278m` and `185.360m` terminal resident persons and pass
the declared market, fiscal, transfer, finite-renewal, and person/head one-step
gates. `tools/run_e5f_perfect_foresight_person_demography_policy.py` is the
isolated full-transition driver. Its H32 convergence stage is in progress; no
corrected transition has yet been promoted.

`tools/build_e5f_person_demography_terminal_spectrum.py` computes the exact
four-year linear operator implied by the frozen terminal demographic closure.
The two policy cases have dominant eigenvalues `0.960754` and `0.961022`, so
their demographic half-lives are roughly `69--70` years and one-percent modal
attenuation requires `116` model periods (`464` years). H32 therefore cannot
be a terminal certification horizon; it is a path-solver stage, H64 is a seed
bridge, and H128 is the first planned terminal test. The fingerprinted packet
is `../../output/model/e5f_person_demography_terminal_spectrum_20260826a/`.

The production model is the one-market sequential-fertility model shown in the
August presentation. It lives in `intergen_eqscale_seq_optimized` and retains
literal parity, independent child maturation, the child-room floor, tenure and
housing-size choice, saving, moving costs, income risk, and warm-glow bequests.
The circulated Stone--Geary one-shot model remains available below as a
fallback; it has not been overwritten.

The paper calibrates this model at a point along a dated transition. For each
candidate, `tools/run_e5f_transition_calibration.py`:

1. derives the old fertility intercept so old-steady-state completed fertility
   is the declared replacement normalization `2.1`;
2. reweights only the age marginal to the 2007 national ACS household-head
   distribution;
3. phases in a linear fertility-preference change through 2023;
4. imposes Census HH-3 household totals and ACS age shares through 2023 as an
   explicit household-formation/migration bridge; and
5. evaluates the full twelve-moment target system along the simulated
   transition, principally on its 2023 state.

The bridge matches household totals and ages by construction; these are inputs,
not fitted paths. Birth cohorts enter a four-slot vintage queue and generate
adult households twenty years later. The conversion is exactly adjusted births
divided by `2.1`; it is an external sex/survival/household-formation
normalization, not a hidden child-death state. Housing clears date by date with
the maintained static supply elasticity. This is a sequence of temporary
equilibria, not a perfect-foresight welfare transition.

The calibrated household block crosses the five-state earnings process with a
three-node discretization of measured permanent-income dispersion, adds one
estimated utility cost paid only when a first birth succeeds, and allows a
separate housing-service floor when the first child is at home. This last
mechanism is nested: setting `hbar_first_child_jump=0` reproduces the previous
child-room floor exactly. The selected
refinement and two exact execution-identity repeats are under:

```
output/model/e5f_transition_ridge_refinement_jump11_polish_r2_20260818/
```

The strict loss is `36.0992231622`. Both evaluations reproduce all eleven
parameters, twelve moments, and the loss to machine precision. The active
target set is `e5_fullhistory_roomsfix_h1_20260817`, with fingerprint
`3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497`.
The two dominant misses are the provisional four-year first-birth rooms
contrast (`0.380` versus `0.720`) and mean occupied rooms (`6.710` versus
`5.780`); together they contribute `26.304` to the loss. The complete row-level target receipt is
`intergen_eqscale_seq_optimized/e5_target_provenance.csv`: six rows are
measured, five remain provisional, and one is an external normalization.

The historical validation packet is:

```
output/model/e5f_transition_historical_validation_jump11_polish_r2_20260818/
```

It labels the Census/ACS population bridge as imposed and keeps fertility,
timing, ownership, rooms, and price/rent comparisons as untargeted holdouts.
The paired post-2023 exercise contains no policy. It compares a national
closed finite-horizon benchmark with an open sensitivity. The value `0.169`
normalizes `M` and `rho` in the old state; those two flows are then fixed while
the realized outside-origin share can change. The closed stationary audit finds no usable
root on the declared price grid: the adjusted-birth renewal ratio lies between
`0.5570040098` and `0.7084847727`. The audited 2183 closed/open population
indices are `0.159890`/`0.378074` relative to 2023, and the corresponding
asset-price ratios are `0.417498`/`0.642780`. The open stationary endpoint has
old-steady-state population scale `0.380834`, price ratio `0.647571`, and
realized outside-origin entry share `0.443763`; these are sensitivity results,
not forecasts.
In a separate conditional profile, both fertility-logit noise scales vary from
`0.25` to `2` times their pre-polish values and the preference decline is
re-solved at every factor to keep 2023 completed fertility on target. No factor
has a closed root; the largest renewal ratio is `0.725974`, and the best fit is
at the unscaled noise. This rules out a simple common-slope repair, not every
possible fertility mechanism.
The open endpoint is conditional on the external entrant flow. It is not a
national forecast.

Rebuild the single canonical no-policy report from the repository root with:

```bash
PYTHONPATH=code/model:code/model/tools code/model/.venv/bin/python \
  code/model/tools/build_e5f_no_policy_transition_report.py \
  --selected-report output/model/e5f_transition_ridge_refinement_jump11_polish_r2_20260818/report/summary.json \
  --repeat-task output/model/e5f_transition_ridge_refinement_jump11_polish_r2_20260818/repeat_001 \
  --repeat-task output/model/e5f_transition_ridge_refinement_jump11_polish_r2_20260818/repeat_002 \
  --refinement-plan output/model/e5f_transition_ridge_refinement_jump11_polish_r2_20260818_plan/candidate_plan.json \
  --target-provenance-csv code/model/intergen_eqscale_seq_optimized/e5_target_provenance.csv \
  --historical-packet output/model/e5f_transition_historical_validation_jump11_polish_r2_20260818 \
  --continuation-packet output/model/e5f_post2023_no_policy_continuations_jump11_polish_r2_20260818_production \
  --slope-frontier-packet output/model/e5f_target_preserving_closed_root_frontier_jump11_production_20260818 \
  --output-dir output/model/e5f_no_policy_transition_report_jump11_polish_r2_20260818
```

The builder does not solve or rescale the model. It fails unless the selected
calibration, both repeats, target provenance, historical classification,
40-date paired continuation, numerical residuals, endpoint accounting, and the
target-preserving slope diagnostic all match their recorded hashes. The current
paper reports no funded policy counterfactual. The older funded driver remains
fail-closed until its renewal and eligibility architecture is migrated and
re-estimated.

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

## Preserved one-shot fallback and one-run review

The parity-verified one-shot implementation is under
`intergen_housing_fertility_optimized/`. It does not replace or redirect the
production sequential package. Its `README.md` and `REFACTOR_REPORT.md` contain the exact
parity commands, benchmark evidence, correctness changes, and promotion gates.

The fixed-M5-parameter owner-ladder density robustness uses that parity-verified
solver to compare the canonical five rungs with one-room and half-room ladders:

```bash
PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/run_intergen_owner_ladder_robustness.py
```

Use `--custom-points N` to run one evenly spaced (N)-point ladder over the
same `[2,10]` support.

The fixed-M5 funded property-tax test solves the equal lump-sum transfer from
the stationary government budget in a rebated 1% baseline, a rebated 2% tax
case, and a rebated 2% tax plus targeted-grant case:

```bash
PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/run_intergen_funded_property_tax_test.py
```

Use `--smoke` for the 40-node testing grid. This driver is a fixed-parameter
mechanism test, not a funded-baseline recalibration. To run the same fiscal
cases at the selected current-M estimate, use:

```bash
PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/run_intergen_funded_property_tax_test.py \
  --profile new-moment \
  --source ../../output/model/intergen_new_moment_unrestricted_overnight_20260723_w2/report/results.json \
  --outdir ../../output/model/intergen_current_m_figures_policy_20260724/funded_policy
```

Those normalized-population results are decomposition rows only. The headline
funded exercise must let entry and stationary population scale adjust. Its
driver is the balanced-budget extension of the established Phase-9b protocol:

```bash
PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/run_intergen_funded_policy_with_entry.py
```

It recovers the outside-entry objects at the rebated 1% baseline, holds them
fixed across policies, and jointly solves the house price and lump-sum rebate.
Use `--smoke` only with a matching 40-node fixed-population packet.

To isolate the property-tax increase without any rebate or grant, starting
from the actual current-M fiscal convention, use:

```bash
PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/run_intergen_tax_no_rebate_with_entry.py
```

This compares a 1% and 2% annual property tax, discards revenue in both cases,
and resolves the Phase-9b entry/scale price closure.

Reproduce the two established draft figures from that exact current-M result
with:

```bash
PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/build_intergen_current_m_draft_figures.py
```

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

For the repaired E5/E5F comparison, `tools/run_e5_repaired_policy_with_entry.py`
accepts either a certified `E5_MATURATION_REPAIR` or `E5F` collector report and
routes the selected mechanism through the funded entry closure. The source must
carry an exact strict repeat. Its production input is the approved
outside-origin entrant share (currently `0.169`), not a model entry probability.
For each funded baseline it computes

\[
q^{E,*}=\frac{1-s^{E,\mathrm{out}}}{B_0/E_0},
\]

and fails if that candidate-specific probability is infeasible. The companion
`intergen_eqscale_seq_optimized/build_e5f_psinneg_iteration_packet.py` builds
the common-target and funded-policy comparison tables and figures, and refuses
policy packets that do not reproduce the empirical share and baseline scale.

The same repaired-E5 driver now also accepts the explicit experimental flag
`--closure-mode quota`; `logit` remains the default so prior rows remain
reproducible. Quota mode replaces the outside value and taste scale with the
baseline identities

\[
\bar R=(1-s^{E,\mathrm{out}})E_0/B_0,
\qquad
\bar M=s^{E,\mathrm{out}}E_0,
\]

then holds \((\bar R,\bar M)\) fixed and solves
\(S E_0=\bar R S B_0(\mathrm{policy})+\bar M\) jointly with the funded house
price and transfer. It fails if \(\bar R>1\) and never evaluates the entry
logit. For the sequential `3+` architecture, quota mode defines (B_0) in the
same top-bin child units as measured completed fertility and writes the raw
three-child flow separately. The quota closure is a fixed-flow policy-reporting
arm, not an equilibrium model of migration.

The live calibration uses the benchmark-normalized outside-option closure:
`P.population_closure = "outside_option_benchmark_normalized"`.

In the benchmark, the model solves for a normalized stationary distribution,
maps the empirical outside-origin share into \(q^E\), calibrates the outside
value to that probability, and sets the outside-born potential entrant mass
mechanically so that benchmark scale is \(S=1\):
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

The August 11 closed reproductive-closure audit evaluates the maintained E5b
entrant requirement, mature-child flow, reproductive residual, and housing
scale map over a fixed-price grid. It is a stationary diagnostic, not a
calibration or transition solver:

```bash
.venv/bin/python tools/audit_closed_reproductive_closure.py --smoke \
  --outdir ../../output/model/closed_reproductive_closure_audit_20260811/smoke
.venv/bin/python tools/audit_closed_reproductive_closure.py
```

The report and visual packet are under
`../../output/model/closed_reproductive_closure_audit_20260811/`; the readable
memo is `../../output/pdf/fertility_population_housing_closure_audit.pdf`.

The August 14 paper update wires the diagnostic `renewal_valve` scale into the
active Markov-income price loop and reruns the retained shared-clock benchmark
and legacy policy cases under both normalized and endogenous stationary scale:

```bash
.venv/bin/python tools/build_population_closure_update.py --solve-policies
```

Outputs are under `../../output/model/population_closure_update/`. At the
retained theta, fixed outside flow `M=E0-B0` nests the normalized benchmark at
scale one. This establishes an open stationary comparison, not a calendar-time
transition; it also assumes outside and locally born entrants share the fixed
entrant-state distribution. This is a historical diagnostic; it is not the
circulated one-shot exercise documented next.

A separate August 15 driver applies the same fixed-inflow accounting directly
to the saved circulated one-shot model, without changing its household states
or Bellman problem:

```bash
.venv/bin/python tools/build_current_one_shot_stationary_closure.py
```

It solves open and closed stationary roots and reruns the circulated funded
purchase-grant policies under the new closure. Outputs are under
`../../output/model/current_one_shot_stationary_closure/`. This is a
fixed-parameter mechanism packet, not a promoted recalibration: the saved
estimate's original first-child housing-response target was later withdrawn.
The packet reports adult-household mass, four-year birth flows, the complete
fiscal ledger, and full fixed-versus-renewal baseline nesting checks.

The generic calendar-time illustration used in the comprehensive theory report
is intentionally separate from the maintained solver. It propagates adult age
cohorts and a four-stage child pipeline, compares the baseline with a
same-initial-state no-fertility-decline path, and verifies the fixed-price
housing-demand decomposition:

```bash
.venv/bin/python tools/build_dynamic_population_housing_paper_illustrations.py
```

Its packet is under
`../../output/model/dynamic_population_housing_paper/`. The reported `price`
is a contemporaneous housing-cost index, not a forward-looking asset price.

See `PLAN.md` for the full implementation and optimization plan.
