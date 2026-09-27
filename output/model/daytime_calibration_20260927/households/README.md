# Saved household audit and three illustrative lives

Source: authenticated overnight main continuation `de_0093`. No equilibrium or
household problem was solved. Standard 17 diagnostic figures are unchanged;
the two PNGs here are supplemental. Source code:
`code/model/tools/inspect_e5f_saved_households.py`.

## Findings

- Current stationary mass is one to 7.2e-14; no negative probability mass.
- Bottom liquid-grid mass is 5.72e-7, top mass zero. No occupied savings policy
  outside the grid, nonpositive consumption, nonfinite consumption or negative
  adjacent-wealth value step in occupied pre-choice states.
- Liquid wealth median is zero; 99th percentile is 18.795, far below the upper
  grid value 3000. Debt is held by 42.919% of households. Current net worth adds
  housing value to financial assets: median 2.309; 99th percentile 27.059.
- Largest owner product contains 23.113% of households (35.236% of owners),
  the largest single owner category, not a majority of owners.
- Ownership reaches 95.352% at age 82. This needs economic/empirical review;
  passing numerical screens does not establish realistic late-life tenure.
- Capped completed children average 1.869. The top-bin weight 3.602359422009
  converts this to 2.100105, matching normalization. These are distinct measured
  objects, not a contradiction. The age CSV reports both explicitly.

## Files and interpretation

`distribution_audit.json`: full-distribution screens and discrete inverse-CDF
quantiles (zero atoms are retained). `distribution_by_age.csv`: age profiles.
`supplemental_distribution.png`: distribution and selected savings-policy zooms.
`three_household_lives.csv` / `.png`: exactly three lives. Seed 20260927; entrant
wealth grid nodes at cumulative entry ranks 10%, 50%, 90%, then remaining entry
states sampled from their conditional distribution. They are illustrative,
not representative averages or handpicked outcomes. All three happened to
survive to the terminal age in this fixed draw. The model's four-year age grid
is retained in the source CSV. The updated display interpolates financial
series linearly between those dates, with dots marking the actual simulated
observations. Housing, children and tenure remain step functions. Age markers
and a separate table at ages 20, 30, 40, 50, 60, 70 and 80 make comparisons easier. These are
display interpolations, not additional annual model simulations. Income is annual gross earnings while working and pension in
retirement; consumption is per model period. Financial quantities are model
units. Children are capped at three (last state represents three or more).

## Exact transition verification

Simulation calls the authenticated runtime's sequential fertility operator,
then its current location/tenure realization. Given the sampled current state,
the original forward operator advances saving, income and child aging with
identity current-choice maps. Thus wealth interpolation lotteries are the same
as those in the model's distribution law, not a new continuous-wealth policy.
Survival is drawn with `P.survival_probs[j]`, matching the calendar operator.
At three occupied age states, composing current realization with conditional
forward evolution exactly reproduces the original forward kernel (L1 = 0).
`simulation_verification.json` records source paths and tests. These checks do
not validate the model's economic assumptions or finer-grid accuracy.

Reproduce from repository root (one process, no solves):

```sh
EXPECTED_UTILITY_OVERNIGHT_SHA256=3b770d8c8c22d2b0449b34a575d6353b063bc015d74ce11016dad7e22ed7ca5e \
E5F_LOCAL_EXECUTION_AUTHORIZATION=tommaso_authorized_20260927_local_primary_continuation_v1 \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MPLBACKEND=Agg \
code/model/.venv/bin/python code/model/tools/inspect_e5f_saved_households.py \
--contract tmp/e5f_overnight_local_20260927/portable/night_launch_v4/primary_continuation/production_contract.json \
--case tmp/e5f_overnight_local_20260927/portable/night_launch_v4/primary_continuation/search/de_0093/case \
--output output/model/daytime_calibration_20260927/households
```

Display-only regeneration (no checkpoint loading or simulation):

```sh
MPLBACKEND=Agg code/model/.venv/bin/python code/model/tools/inspect_e5f_saved_households.py --plot-only --output output/model/daytime_calibration_20260927/households
```

`three_household_lives_annual_display.csv` and
`three_household_lives_selected_ages.csv` contain the interpolated display
values. `display_interpolation.json` records their interpretation and the
unchanged source CSV hash. All 51 original observations match exactly.
