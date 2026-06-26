# Intergen Housing Smoothing Audit, June 26 2026

This note records a narrow diagnostic pass on the one-market intergenerational
housing block. The goal was to test whether the jagged aggregate consumption
function is mainly a tenure-logit smoothness problem, a discrete tenure value
crossing problem, or something deeper in grid/policy accounting. No solver or
calibration logic was changed.

## Inputs

- Baseline cache:
  `output/model/intergen_fixedstats_overnight_review_20260626/best_de_g044_i022_packet/solution_cache.pkl`
- Dense-grid GE cache:
  `output/model/intergen_fixedstats_overnight_review_20260626/housing_block_grid_audit_20260626/ge_dense_core_nb120/solution_cache.pkl`
- Target set:
  `candidate_replacement_young_old_roomgap_v1`
- Dense-grid audit configuration:
  `Nb=120`, liquid grid core roughly `[-8,11]`, `H_own=[2,4,6,8,10]`,
  `hR_max=6.0`, `interp_method=linear`, `max_iter_eq=10`

## Value-Gap Audit

New tool:
`code/model/tools/audit_intergen_active_tenure_values.py`.

The tool reconstructs conditional tenure-branch values from the saved policies
and next-period value function, then applies the same transaction and
feasibility accounting as the tenure-choice kernel. It is a read-only check of
the saved solution object.

Commands run:

```bash
code/model/.venv/bin/python code/model/tools/audit_intergen_active_tenure_values.py \
  --cache output/model/intergen_fixedstats_overnight_review_20260626/best_de_g044_i022_packet/solution_cache.pkl \
  --outdir output/model/intergen_fixedstats_overnight_review_20260626/best_de_g044_i022_packet/active_tenure_value_audit \
  --wealths 0,0.146514,4.2487 --top-per-wealth 3 --extra-top-mass 3

code/model/.venv/bin/python code/model/tools/audit_intergen_active_tenure_values.py \
  --cache output/model/intergen_fixedstats_overnight_review_20260626/housing_block_grid_audit_20260626/ge_dense_core_nb120/solution_cache.pkl \
  --outdir output/model/intergen_fixedstats_overnight_review_20260626/housing_block_grid_audit_20260626/ge_dense_core_nb120/active_tenure_value_audit \
  --wealths 0,0.146514,4.2487 --top-per-wealth 3 --extra-top-mass 3
```

Key validation:

- Baseline stored-vs-reconstructed tenure probability max difference:
  `1.91e-08`.
- Dense-grid stored-vs-reconstructed tenure probability max difference:
  `2.65e-08`.

Interpretation: the diagnostic matches the stored solution closely enough to
trust the value-gap readout.

Dense-grid value-gap findings:

- At the high-mass near-zero states (`b=0` and `b=0.146514`), renting is the
  only relevant option in practice. Owner branches are infeasible or far below
  the renter value. This means the big high-mass kink near the entry liquid
  wealth is not mainly a smooth owner/renter indifference crossing.
- At `b=4.214286`, age `38`, childless renter, the renter branch beats `own_H4`
  by only `0.0222` value units and the stored owner probability is `0.1134`.
  This is a genuine close tenure margin.
- At `b=4.214286`, age `38`, one-child renter, `own_H8` beats `own_H6` by only
  `0.0135` value units, with stored probabilities about `0.794` on `H8` and
  `0.206` on `H6`. This is a genuine close owner-rung margin.

So some mid-wealth jaggedness is discrete tenure/rung economics, but the
high-mass near-zero consumption drop is not explained away by near-indifference.

## Kappa Sweep

New tool:
`code/model/tools/run_intergen_kappa_smoothing_audit.py`.

Command run:

```bash
code/model/.venv/bin/python code/model/tools/run_intergen_kappa_smoothing_audit.py \
  --source-cache output/model/intergen_fixedstats_overnight_review_20260626/housing_block_grid_audit_20260626/ge_dense_core_nb120/solution_cache.pkl \
  --outdir output/model/intergen_fixedstats_overnight_review_20260626/housing_block_kappa_smoothing_audit_20260626 \
  --kappas 0,0.005,0.01,0.02,0.05
```

The script varies only `tenure_choice_kappa`; theta, targets, owner ladder,
renter cap, and dense grid are inherited from the source cache. Each case writes
a `solution_cache.pkl` and standard `housing_block_audit/` packet. Full target
moment output is in
`output/model/intergen_fixedstats_overnight_review_20260626/housing_block_kappa_smoothing_audit_20260626/target_moments_by_kappa.csv`.

Summary:

| kappa | loss | own | young own | old own | TFR | entry c drop | entry mass |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 14.6247 | 0.6037 | 0.1784 | 0.9752 | 1.8875 | -0.3436 | 0.1029 |
| 0.005 | 14.5271 | 0.6032 | 0.1774 | 0.9740 | 1.8865 | -0.3411 | 0.1027 |
| 0.010 | 15.0154 | 0.6030 | 0.1724 | 0.9743 | 1.8860 | -0.3426 | 0.1025 |
| 0.020 | 15.9149 | 0.5986 | 0.1579 | 0.9719 | 1.8869 | -0.3412 | 0.1021 |
| 0.050 | 15.5045 | 0.5918 | 0.1614 | 0.9002 | 1.8932 | -0.3076 | 0.0994 |

Interpretation:

- Increasing tenure-choice smoothness does not remove the high-mass entry-wealth
  consumption drop. The drop at `b_entry_fixed≈0.146514` stays around `-0.34`
  for `kappa` from `0` through `0.02`, and remains `-0.308` even at `0.05`.
- `kappa=0.05` changes economics materially: old-age ownership falls from about
  `0.974` to `0.900`, renter large-room share rises to `0.111`, and owner
  large-room share falls to `0.636`. That is a substantive discrete-choice
  specification change, not a harmless numerical smoothing fix.
- `kappa=0.005` and deterministic tenure have slightly lower fixed-theta loss
  than the current `kappa=0.01`, but neither fixes the near-zero kink. This is
  not enough to justify changing the model before a recalibration.

## Implication

The current evidence points away from a simple solver failure or a simple
tenure-logit smoothness failure. The next housing-block audit should focus on:

1. The liquid grid and entry/liquid wealth atom around `b_entry_fixed`.
2. The budget/accounting reason why average consumption drops sharply exactly at
   the high-mass entry node even when owners are infeasible.
3. Owner-rung discreteness around the mid-wealth purchase region, especially
   `H4/H6/H8` thresholds for parents.
4. Whether the owner ladder should be refined locally only after retuning, since
   removing or shifting the low rung changed ownership moments materially in the
   earlier grid audit.

