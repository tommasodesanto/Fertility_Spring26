# Reproducible saved-array review

`code/model/tools/analyze_e5f_rental_wedge_saved_arrays.py` reads the six retained `policy_arrays.npz` files sequentially with NumPy only. It verifies each NPZ SHA against its outer case receipt before analysis, requires the renderer axis order `(wealth, tenure, location, age, income, children_ever_born, child_state)`, checks bitwise common `g_pre`, and writes `lead_saved_array_review.json` outside the immutable experiment results tree.

The remote read completed with six cases and a common finite-value mask of 2,717,403 states. The renter statistic uses `g_current[:, 0, ...]` and `hR_pol[:, 0, ...]`; this agrees with every outer receipt's renter mass. The prior hand-entered supplement is preserved as `results/saved_array_supplement_pre_axis_fix.json`, with the axis correction recorded in `results/saved_array_supplement.json`.

Reproduction command:

```text
python code/model/tools/analyze_e5f_rental_wedge_saved_arrays.py \
  --results /scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_exhaustive_v1 \
  --output /scratch/td2248/projects/Fertility_Spring26_specification_20260920/lead_saved_array_review/output
```

This is a saved-array read only; it performs no household solve, replay, plot regeneration, or gate change.
