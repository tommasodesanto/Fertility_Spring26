# Children ever born by household age: fixed-price credit comparison

The [count-share figure](child_count_shares_by_age.png) and [CDF figure](child_count_cdf_by_age.png)
compare the saved \(\phi=0.8\) and \(\phi=1.0\) household distributions. The CDF
shows \(P(N\leq k\mid\text{age cell})\) for \(k=0,1,2\), plus the \(\phi=1.0\)
minus \(\phi=0.8\) difference in percentage points. The complete observations
are in [children_by_age.csv](children_by_age.csv) and
[differences_by_age.csv](differences_by_age.csv).

\(N\) is **children ever born**, the sixth axis of the saved seven-dimensional
`g` array. This differs from the seventh-axis count of children currently at
home. All households, including the childless, enter each age denominator.
Each arm is normalized by its **own** realized household mass at that age. The
saved `g` is a post-birth, current-period household distribution: births are
assigned within the age period, then current location and tenure are realized.
Age means the modeled reproductive household member's age: 17 four-year cells
starting at 18, 22, ..., 82. For example, the point marked 42 covers ages
42–45; the point marked 22 covers ages 22–25, not exact age 25 or an
interview-age projection. The top child-count state is labeled **3 or more**.

At ages 42–45, the \(\phi=0.8\) shares of households with 0, 1, 2, and 3+
children ever born are 18.888%, 14.068%, 28.448%, and 38.596%. Under
\(\phi=1.0\), they are 19.611%, 13.958%, 28.133%, and 38.298%. Thus
childlessness is 0.723 percentage points higher in the relaxed-credit arm at
that age. The distribution is constant after the modeled childbearing window.
These whole-model age shares should not be read as a recalculation of an
empirical target: the target observers can use different age and sample
mappings. This is a **fixed-price, partial-equilibrium** comparison. The financed
share \(\phi=1\) removes the nominal down-payment requirement, while the
model's collateral floor and other screens remain. It is not a claim of
unrestricted borrowing or a recalibrated equilibrium.

The script uses `g`, the realized current cross-section constructed after
`realize_current_cross_section` in the active production
[`distribution.py`](../../../../../../../../code/model/production/engine/distribution.py).
The source's age formula is `age_start + j*da`, and executed `P` has
`age_start=18`, `da=4`, `J=17`, `n_parity=4`, and
`use_postdecision_current_distribution=true`. Source array hashes match the
two arms in the parent [`completed.json`](../completed.json). The
[provenance and checks](checks_and_provenance.json) also record total mass,
minimum cell and age mass, age-share row sums, monotone CDFs, and equality of
the current and beginning distributions' age/count margins. No solve was run.

Regenerate from the saved arrays:

```sh
output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python \
  output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/credit_mechanism/credit_relaxation/children_by_age/plot_children_by_age.py
```

For another matched \(\phi\) pair, pass `--phi-080 PATH`, `--phi-100 PATH`,
`--output DIR`, and optionally `--title TEXT`. Each `PATH` must point to a
`solution_arrays.npz` beside its executed `executed_P.json`.
