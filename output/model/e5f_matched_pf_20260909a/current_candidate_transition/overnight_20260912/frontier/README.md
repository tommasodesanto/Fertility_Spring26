# Smaller-shock native forecasts

Both runs clear their finite six-date housing and PAYGO roots and retain all household gates. Neither is a fitted historical sequence or a horizon-verified solution. The original historical initial calibration and all empirical targets/weights remain fixed.

| Preference decline | First-window model fertility | Data | Gap | Maximum housing residual | Maximum PAYGO residual | Standard graphs |
|---:|---:|---:|---:|---:|---:|---:|
| -0.0225 | 1.90752953 | 1.97487500 | -0.06734547 | 1.98e-05 | 2.194e-07 | 17 |
| -0.0275 | 1.86726183 | 1.97487500 | -0.10761317 | 2.38e-06 | 1.967e-08 | 17 |

The smaller decline gets closer to the first empirical window. Failures of still smaller shocks at the old starting price/pension pair do not prove that no terminal equilibrium exists. A separate pinned warm-start experiment tests -0.01414 and -0.0175 from the nearby converged terminal, without modifying the household model, fiscal closure or numerical tolerances.

Array 17488021 runs those terminal probes (30 minutes each). Conditional array 17488254 reproduces a passed terminal and automatically runs its short forecast (90 minutes each); a failed terminal skips its dependent forecast. The original long transition remains active.
