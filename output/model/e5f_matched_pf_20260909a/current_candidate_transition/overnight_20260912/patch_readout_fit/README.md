# Refreshed patch figure packet

This packet replays the selected Torch forecast `forecast_6_0` from job
`17499630`, using the approved initial checkpoint and the stationary 2011,
2015, and 2019 source fits retained from the original patch packet.

The replay source receipt is `source/root_receipt.json` (SHA-256
`8b69b1e874e0f69030035d18a9681f6ec9a1b63565df628201646a9d66d47136`). The
exact replay check is recorded in `source/verification.json`: status `PASS`,
maximum absolute aggregate discrepancy `0.0`, and runtime 191.58 seconds.
The selected 2019 window has psi `0.10239514522037683`, model fertility
`1.6416808867927979`, and target `1.64575`. The finite forecast's terminal
distance check remains failed and `horizon_verified=false`; this packet does
not promote the forecast to a completed history or long-run result.

The exact Torch collector command was:

```text
python -B collect_e5f_patch_readout_refresh.py --batch /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/stationary_2019_pf_test_20260912 --plan /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/stationary_2019_pf_test_20260912/plan.json --forecast-stage /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/finite_sequences_20260912/patch_retry/results/arm_3/trial_00_2019/forecast_6_0 --out /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/finite_sequences_20260912/patch_readout_fit
```

The local figure command was:

```text
python3 -B code/model/tools/build_e5f_patch_readout.py --base output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/patch_readout_fit --pdf output/pdf/e5f_patch_review_fit.pdf
```

The resulting review has seven rendered pages in
`output/pdf/e5f_patch_review_fit.pdf`; individual PDF and PNG figures are
under `figures/`.
