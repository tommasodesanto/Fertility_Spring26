# Reproducible housing profiles and model overlay

Command:

```bash
code/model/.venv/bin/python code/empirical/housing/plot_initial_housing_profiles.py --input output/model/native_financing_diagnostic_20260919/specification_followup/housing_profiles_v1/full --model-dir output/model/native_financing_diagnostic_20260919/specification_followup/housing_profiles_v1/model_collection_shared_v3 --output output/model/native_financing_diagnostic_20260919/specification_followup/housing_profiles_v1/reproducible_overlay_v1
```

The command emits the empirical-only plots and, with `--model-dir`, physical model overlays. Model rooms use renter `hR_pol` and owner `H_own`, capped at 9 before `g_current` weighting, in age cells `18+4j`. The model has no DUE classification; ownership is all-household/all-structure. Checkpoint family prices and preferences may differ, so this is descriptive and not a causal effect.
