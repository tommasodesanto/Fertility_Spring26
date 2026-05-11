# Population Closure Experiment Copy

This folder contains the initial experiment notes for the paper-facing
population closure. The first pass is now wired into the live Python package
behind the guarded mode
`P.population_closure = "outside_option_local_births"`. The default Python
solver path remains `P.population_closure = "normalized"`.

## Proposed Solver Change

The current Python port mirrors the MATLAB solver:

- `P.E_total = 1 / P.J`
- `P.entry_by_loc = P.E_total * P.entry_shares`
- `forward_distribution` rescales the full cross-section to `P.N_target`
- housing demand is divided by `P.N_target`

For the outside-option closure, use the guarded mode
`P.population_closure = "outside_option_local_births"`, keeping the current
behavior as the default parity mode.

The new mode should:

1. Set `P.normalize_population_mass = False`.
2. In `forward_distribution`, skip the `P.N_target / tm` rescaling.
3. Report housing demand in aggregate room-equivalent levels, not per unit of
   normalized population.
4. After each Bellman/distribution pass, compute entry values
   `V_1(b_entry, renter, location, childless, no_child_state)`.
5. Use `stats.entrants_mature_total` as the city-born mature-child flow.
6. Update `P.entry_by_loc` with
   `outside_entry_flow + local_birth_entry_weight * entrants_mature_total`
   times the city-location logit probabilities.
7. Update `P.E_total = sum(P.entry_by_loc)` and
   `P.entry_shares = P.entry_by_loc / P.E_total`.

The existing `stats.entrants_mature_total` object is already close to the needed
model object. The substantive change is to stop normalizing it away and to use it
to update entry masses, not only entry shares.

## 2026-05-05 Smoke Results

Artifacts:

- `benchmarks/fast_ge_iter1_2026_05_05_postclosure_default.json`
- `benchmarks/fast_ge_iter1_2026_05_05_postclosure_default_logupdate.json`
- `benchmarks/fast_ge_iter1_no_mass_rescale_2026_05_05.json`
- `benchmarks/accounting_scale_fast_iter1_2026_05_05.json`
- `benchmarks/accounting_scale_fast_live_stop_2026_05_05.json`
- `benchmarks/accounting_scale_fast_live_stop_fixed_vo_2026_05_05.json`
- `benchmarks/accounting_scale_psi_sensitivity_iter1_2026_05_05.json`
- `benchmarks/scaled_equilibrium_fast_iter80_fixed_vo_2026_05_05.json`
- `benchmarks/scaled_equilibrium_psi_sensitivity_fixed_vo_2026_05_05.json`
- `benchmarks/scaled_equilibrium_outside_flow_sensitivity_fixed_vo_2026_05_05.json`
- `benchmarks/outside_closure_fast_iter8_2026_05_05.json`
- `benchmarks/outside_closure_fast_damped_iter30_2026_05_05.json`
- `benchmarks/outside_closure_fast_calibrated_iter50_2026_05_05.json`

Checks:

- Default mode is unchanged: post-edit versus pre-edit one-iteration fast GE has
  max absolute moment difference `0`.
- Default mode still matches the saved MATLAB one-iteration fast reference with
  max absolute moment difference `0.0026325232`, the old documented gap.
- Turning off population rescaling while keeping fixed `E_total=1/J` changes the
  one-iteration fast GE moments only at floating-point mass error
  (`max_abs_diff = 1.43e-12`). This confirms that the explicit normalization is
  not doing substantive work when entry mass is fixed at the stationary cohort
  size.
- The fast accounting path `dt_cp_model.cli accounting-scale` calibrates
  `V^O=-42.56021268867121` from the baseline one-iteration fast solve. It
  recovers implied population `1.000000000000002`, mature city-born children per
  entrant `0.8361542407`, and aggregate housing demand
  `[6.9161589955, 0.1372367846]`.
- With the fast setup under the live stopping rule (`max_iter_eq=120`), the same
  accounting path calibrates `V^O=-35.98609192692343`, recovers implied
  population `0.9999999999999994`, mature city-born children per entrant
  `0.8404686131`, and aggregate housing demand
  `[4.3923897537, 2.6921642495]`.
- Holding that same `V^O` fixed, a one-iteration `psi_child` perturbation moves
  implied scale in the expected direction:
  `psi_child=0.03` gives TFR `1.643` and scale `0.889`;
  baseline `psi_child=0.07` gives TFR `1.702` and scale `1.000`;
  `psi_child=0.09` gives TFR `1.730` and scale `1.063`.
- Important bug fix: the first accounting-scale-price implementation
  accidentally recalibrated `V^O` whenever the helper was called with
  `outside_value=None`. That mechanically forced the implied scale back to one
  inside counterfactual runs. The helper now uses `P.outside_value` by default;
  recalibration happens only when explicitly requested by the CLI baseline
  calibration path.
- The preferred endogenous-housing-demand closure is now
  `P.population_closure = "accounting_scale_prices"`. It keeps the distribution
  normalized for speed, computes the implied stationary population scale from
  the outside-option accounting equation, and clears housing markets with
  aggregate demand `N(p) * d_i(p)`.
- With fixed `V^O=-35.98609192692343`, the corrected scaled-price baseline
  reaches the fast solver's soft-accept region with `best_eq_error=0.00519`,
  TFR `1.709`, prices `[0.52215, 0.67307]`, and implied population `1.00197`.
- Corrected fertility-preference scaled-price sensitivity with the same fixed
  `V^O`:
  `psi_child=0.03` gives TFR `1.656`, scale `0.9896`;
  baseline `psi_child=0.07` gives TFR `1.709`, scale `1.0020`;
  `psi_child=0.09` gives TFR `1.734`, scale `1.0021`;
  `psi_child=0.12` gives TFR `1.772`, scale `1.0155`.
- Corrected outside-flow scaled-price sensitivity with the same fixed `V^O`:
  `M^O=0.75/J` gives scale `0.9883` but had not fully settled by 80 iterations
  (`best_eq_error=0.0303`);
  `M^O=1/J` gives scale `1.0020`;
  `M^O=1.25/J` gives scale `1.0200`.
- The fully joint outside-option/local-birth GE mode now uses log entry-mass
  updates, but it should still be treated as experimental. With a calibrated
  outside value and `50` fast-setup iterations, it is finite and improving but
  not accepted (`best_eq_error = 0.3305`). The robust production path for now is
  the fast accounting scale from a solved normalized equilibrium.
- The accounting-scale object now reports explicit stationarity residuals:
  `stationary_entry_residual`, `stationary_scale_residual`, and
  `stationary_entry_relative_residual`.
- `accounting_scale_prices` is guarded against accidental use of the placeholder
  `outside_value=0`. A direct solve must set `P.outside_value` and
  `P.outside_value_is_calibrated = True`, or deliberately set
  `P.allow_uncalibrated_outside_value = True`.
- `tools/check_population_closure.py` now runs the regression checks for this
  block. On 2026-05-05 it verifies that the default normalized one-iteration
  path is unchanged against `fast_ge_iter1_2026_05_05_preclosure.json`, the
  calibrated scale is exactly one up to numerical tolerance, uncalibrated
  scaled-housing mode is rejected, fixed-outside-value fertility sensitivity
  moves city scale from `0.889394` to `1.16818`, and a three-iteration
  scaled-housing smoke solve remains finite.
- The low-outside-flow case is not a hidden accounting failure. A longer
  `160`-iteration run with fixed `V^O=-35.98609192692343` keeps exact
  stationarity residuals but remains in a price-map cycle with
  `best_eq_error=0.03025`. A force-full `80`-iteration run also cycles
  (`best_eq_error=0.03936`), so the issue is the non-smooth discrete
  tenure/housing price map rather than stale Howard evaluation. The solver now
  reports `accepted`, `strict_converged`, `convergence_reason`,
  `best_eq_iter`, and `final_eq_error` so this case is no longer silently
  treated as converged.

## Calibration Implication

This should not change the household problem. The parameters
`outside_entry_flow`, `outside_value`, and `local_birth_entry_weight` only close
the level of the stationary city population. They can be calibrated or normalized
to match the baseline population/housing-stock scale, then held fixed in
counterfactuals.
