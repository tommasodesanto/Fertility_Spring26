# Income-feasibility frontier (2026-07-18, Torch job 14170408)

The (sigma, floor, c_bar_0) feasibility frontier at the FIXED M5 winner — the computed version of the paper footnote's claim that the calibration sits on the feasibility frontier. 28 cells, annual persistence held at the calibrated rho=0.9602, tight M5 evaluator, Torch platform (canonical). Variant A@0.0645 reproduces the canonical M5 loss exactly: 9.044422069071352.

Variants: A = no floor, c_bar_0 = M5 (0.315 annual); B = no floor, c_bar_0 = 0.10 annual; C = measured floor (G = 0.13 + 0.10/child annual, debt-blind means test), c_bar_0 = M5; D = measured floor, c_bar_0 = 0.10.

| sigma_annual | A | B | C | D |
| ---: | :--- | :--- | :--- | :--- |
| 0.0645 | ok loss=9.04 | ok loss=1333.08 | ok loss=9.04 | ok loss=1333.08 |
| 0.09 | DEAD@forward_age_22 | ok loss=925.25 | DEAD@forward_age_22 | ok loss=925.24 |
| 0.12 | DEAD@forward_age_22 | ok loss=639.81 | DEAD@forward_age_22 | ok loss=639.86 |
| 0.15 | DEAD@forward_age_22 | ok loss=509.43 | DEAD@forward_age_22 | ok loss=513.10 |
| 0.18 | DEAD@forward_age_22 | ok loss=459.97 | DEAD@forward_age_22 | ok loss=466.22 |
| 0.20 | DEAD@forward_age_22 | DEAD@forward_age_22 | DEAD@forward_age_22 | ok loss=452.52 |
| 0.22 | DEAD@forward_age_22 | DEAD@forward_age_22 | DEAD@forward_age_22 | DEAD@forward_age_22 |


## Readings

1. **The current calibration sits against the frontier even harder than documented**: without a floor, at the calibrated persistence, sigma = 0.09 already kills it (A row; dead mass 2.0e-4 at age 22 — the bottom state of the 0.09 grid pays 0.295 of mean income annually, below the 0.33 bundle plus entry-debt service).
2. **The measured floor buys nothing at the current c_bar_0**: C = A cell for cell (identical losses, identical dead cells, identical boundary census). A guarantee of 0.13 cannot clear a 0.33 bundle — the decision memo's units point, now as a table row.
3. **The bundle is the binding object over most of the range**: lowering c_bar_0 to 0.10 alone (B) extends the frontier from <0.09 to 0.18. The floor's entire feasibility contribution is the last step, 0.18 -> 0.20 (D); inside the frontier B and D are near-identical (receipt <= 0.25%, outlays <= 3bp).
4. **The tail closes smoothly in sigma along D**: estate p90/50 = 1.81, 2.00, 2.18, 2.40, 2.75, 3.07 for sigma = 0.0645...0.20 (target 3.45; M5 baseline 1.75). No heterogeneity anywhere.
5. **The yliq trade-off is now a measured curve**: young liquid wealth crosses its 0.179 target near sigma ~ 0.09 (where the tail is only ~2.0) and overshoots to 1.33 at sigma = 0.20 (where the tail is 3.07). Sigma alone cannot deliver both; that is the refit's problem statement.
6. Fixed-theta loss falls monotonically in sigma along the feasible low-bundle rows (1333 at 0.0645 -> 452 at 0.20): risk substitutes for the lost bundle discipline (rebuilds young saving, suppresses fertility: TFR 2.87 -> 2.16, childless 0.05 -> 0.28). The honest-sigma end is the best fixed-theta point in the family.
7. D's own frontier at 0.22 is a knife-edge dynamic failure (dead mass < 5e-7, census slack +0.013 with the transfer paid): the debt-amortization schedule, not the flow bundle, binds — the M6 forbearance-margin note again.

## Key moments along the feasible rows

| cell | p90/50 | yliq_2535 | tfr | childless | own | receipt | outlays/inc |
|---|---:|---:|---:|---:|---:|---:|---:|
| A_s0.0645 | 1.751 | +0.328 | 2.034 | 0.189 | 0.658 | 0.0000 | 0.00000 |
| B_s0.0645 | 1.808 | +0.153 | 2.873 | 0.050 | 0.742 | 0.0000 | 0.00000 |
| B_s0.09 | 1.996 | +0.324 | 2.790 | 0.059 | 0.743 | 0.0000 | 0.00000 |
| B_s0.12 | 2.177 | +0.676 | 2.569 | 0.092 | 0.750 | 0.0000 | 0.00000 |
| B_s0.15 | 2.398 | +1.172 | 2.307 | 0.162 | 0.717 | 0.0000 | 0.00000 |
| B_s0.18 | 2.753 | +1.646 | 2.087 | 0.255 | 0.660 | 0.0000 | 0.00000 |
| C_s0.0645 | 1.751 | +0.328 | 2.034 | 0.189 | 0.658 | 0.0000 | 0.00000 |
| D_s0.0645 | 1.808 | +0.153 | 2.873 | 0.050 | 0.742 | 0.0000 | 0.00000 |
| D_s0.09 | 1.996 | +0.324 | 2.790 | 0.059 | 0.743 | 0.0000 | 0.00000 |
| D_s0.12 | 2.177 | +0.676 | 2.569 | 0.092 | 0.750 | 0.0001 | 0.00000 |
| D_s0.15 | 2.401 | +0.928 | 2.367 | 0.158 | 0.742 | 0.0012 | 0.00011 |
| D_s0.18 | 2.753 | +1.172 | 2.244 | 0.238 | 0.712 | 0.0025 | 0.00027 |
| D_s0.20 | 3.065 | +1.328 | 2.157 | 0.277 | 0.683 | 0.0091 | 0.00054 |

Full 15-row target-fit tables and censuses per cell: results.jsonl / results.json. Cross-checks: A@0.0645 = canonical M5 record; D@0.20 = probe cell 4 (platform band ~1e-12 vs local). Reproduce:

```
cd code/cluster && sbatch submit_intergen_income_frontier.sh   # or locally: PYTHONPATH=code/model python3 code/model/tools/run_intergen_income_feasibility_frontier.py
```

