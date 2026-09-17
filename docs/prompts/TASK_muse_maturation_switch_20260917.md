# Task for Muse Spark (OpenCode): parent-age maturation switch, default off

Repository: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`,
branch `main`. Read `CLAUDE.md` first. Python: `code/model/.venv/bin/python`.
Do not change any default behavior: the switch must be off by default and the
off value must reproduce the current arrays bit for bit. Do not touch targets,
calibration profiles, drivers under `code/model/tools/`, or anything under
`output/`. Do not use git.

## What to build

In `code/model/intergen_eqscale_seq_optimized/parameters.py`, the child
maturation law is built by `make_independent_child_count_transition_matrix`
(around line 888): each child at home leaves independently with probability
`mu = 1 / (A_m / period_years)` per period (`A_m = 18` years, period 4 years,
so `mu = 4/18`), and the transition is `Binomial(m, 1 - mu)` on the count of
children at home `m` (child state index `cs`), the same at every age. It is
applied in `solver.py` by `apply_child_aging` (around lines 6975–7009), which
indexes the transition by the number of children ever born only, and in the
population step of `code/model/tools/run_e5f_open_population_transition.py`
(around lines 847–895). A newborn enters the first maturation draw.

Add a default-off switch `child_maturation_mode` with values:

- `"constant"` (default): current behavior, unchanged code path.
- `"parent_age"`: the exit probability depends on the parent's model age `a`
  (the age index `j`): `mu(a) = mu_young` for `a < a_rise`, rising linearly
  to `1` at `a_full`, and `1` for `a >= a_full`. New parameters with defaults
  `mu_young = 0.05`, `a_rise = 34`, `a_full = 62` (tuned so that, with a newborn exemption, the implied years at home are 19.7 and 76 percent of child-years fall while the parent is under 46, matching an 18-year benchmark; see docs/model/f3_maturation_tradeoff_note_20260917.md) (calendar ages; convert to
  period indices the way the package does elsewhere). In addition, a child born
  in the current period is exempt from that period's maturation draw. Implement
  the exemption with a one-bit flag `born_this_period` carried alongside the
  child state, or, if the state layout makes that expensive, by applying the
  draw to `m - d` where `d` is the birth indicator available at the point where
  births are realized (`solver.py` around lines 2719–2768 computes the
  post-birth values; the birth destination is `birth_destination_child_state`,
  line ~236). State clearly which of the two you implemented and why.

The switch must be honored consistently in: the Bellman continuation
(`apply_child_aging` and wherever `Pa` is used in `solver.py`, including the
KFE / distribution operators around lines 6001–6100 and 7056–7079), and the
population transition in `run_e5f_open_population_transition.py`. Grep for
every use of the child transition matrix (`Pa`, `child_transition`,
`make_independent_child_count_transition_matrix`) and cover each one.

## Tests (put them in `code/model/intergen_eqscale_seq_optimized/tests/`)

1. With `child_maturation_mode="constant"`, every array the switch touches is
   bitwise identical to the current code on a tiny configuration
   (`Nb=20, J=17, n_parity=4`); assert with `np.array_equal`, not `allclose`.
2. With `"parent_age"`, the transition matrix at each age has rows summing to
   one, exit probability equal to `mu(a)` per child, and equals the identity
   for the newborn component when the exemption applies.
3. The population operator conserves mass under both modes.
4. A smoke that solves the tiny configuration in both modes without error.

## Deliverable

A short report (under 40 lines): files changed with line ranges, which
exemption implementation was used, test results, and the exact override keys
to set in a spec (`child_maturation_mode: parent_age`, `mu_young`, `a_rise`,
`a_full`) so the sandbox at `code/model/sandbox/` can run it with
`make sandbox SPEC=...`. Do not run the sandbox yourself.
