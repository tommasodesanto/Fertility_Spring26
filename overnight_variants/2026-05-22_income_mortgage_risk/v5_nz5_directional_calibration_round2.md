# V5 Nz=5 Directional Calibration Audit

This is a small same-grid calibration probe for the isolated V5
benchmark-normalized outside-option closure. It is not an optimizer. The goal
is to test whether the bad `Nz=5` moments are locally movable in the right
direction before running a larger search.

Baseline reference: `income_mortgage_risk_v5_hank_z_outside_closure_nz5.log`
with loss `311.3`.

## Summary

| Case | OK | iters | GE err | seconds | loss | loss minus baseline | improved moments |
|---|---:|---:|---:|---:|---:|---:|---:|
| fertility_logit_mid | True | 12 | 0.000314 | 102.21 | 286.19 | -25.11 | 11/19 |
| fertility_mid_finance_high | True | 11 | 0.000218 | 97.29 | 252.09 | -59.21 | 12/19 |
| beta_mid_finance_high | True | 15 | 0.000191 | 126.17 | 275.52 | -35.78 | 11/19 |
| fertility_mid_finance_beta_mid | True | 12 | 0.00025 | 121.18 | 253.22 | -58.08 | 13/19 |
| center_mild | True | 20 | 0.000473 | 166.21 | 312.46 | 1.16 | 10/19 |
| alpha_cons_high | True | 11 | 0.00027 | 114.08 | 206.95 | -104.35 | 14/19 |
| alpha_high_fertility_mid_finance | True | 16 | 0.00041 | 107.73 | 166.72 | -144.58 | 16/19 |

Lowest-loss probe is `alpha_high_fertility_mid_finance` with loss `166.722` versus baseline `311.3`.

## Key Moment Directions

Positive direction entries mean the probe reduced the absolute target gap
relative to the `Nz=5` baseline.

### fertility_logit_mid
| Moment | Target | Baseline | Probe | Direction |
|---|---:|---:|---:|---:|
| `tfr` | 1.700 | 1.558 | 1.734 | +0.108 |
| `childless_rate` | 0.150 | 0.351 | 0.283 | +0.068 |
| `mean_age_first_birth` | 26.000 | 34.793 | 34.468 | +0.325 |
| `tfr_gradient` | 0.133 | -0.239 | -0.228 | +0.011 |
| `own_rate` | 0.627 | 0.575 | 0.562 | -0.013 |
| `own_gradient` | 0.170 | -0.082 | -0.032 | +0.050 |
| `young_liquid_wealth_to_income` | 0.600 | 1.215 | 1.367 | -0.152 |
| `center_share_nonparents` | 0.494 | 0.316 | 0.311 | -0.005 |
| `old_age_own_rate` | 0.863 | 0.709 | 0.711 | +0.002 |
| `old_age_parent_childless_gap` | 0.070 | 0.291 | 0.291 | -0.000 |

### fertility_mid_finance_high
| Moment | Target | Baseline | Probe | Direction |
|---|---:|---:|---:|---:|
| `tfr` | 1.700 | 1.558 | 1.737 | +0.106 |
| `childless_rate` | 0.150 | 0.351 | 0.282 | +0.069 |
| `mean_age_first_birth` | 26.000 | 34.793 | 34.460 | +0.333 |
| `tfr_gradient` | 0.133 | -0.239 | -0.236 | +0.003 |
| `own_rate` | 0.627 | 0.575 | 0.571 | -0.003 |
| `own_gradient` | 0.170 | -0.082 | -0.018 | +0.064 |
| `young_liquid_wealth_to_income` | 0.600 | 1.215 | 1.366 | -0.151 |
| `center_share_nonparents` | 0.494 | 0.316 | 0.311 | -0.005 |
| `old_age_own_rate` | 0.863 | 0.709 | 0.730 | +0.021 |
| `old_age_parent_childless_gap` | 0.070 | 0.291 | 0.258 | +0.033 |

### beta_mid_finance_high
| Moment | Target | Baseline | Probe | Direction |
|---|---:|---:|---:|---:|
| `tfr` | 1.700 | 1.558 | 1.551 | -0.007 |
| `childless_rate` | 0.150 | 0.351 | 0.335 | +0.016 |
| `mean_age_first_birth` | 26.000 | 34.793 | 34.257 | +0.536 |
| `tfr_gradient` | 0.133 | -0.239 | -0.219 | +0.020 |
| `own_rate` | 0.627 | 0.575 | 0.521 | -0.054 |
| `own_gradient` | 0.170 | -0.082 | -0.084 | -0.002 |
| `young_liquid_wealth_to_income` | 0.600 | 1.215 | 1.097 | +0.118 |
| `center_share_nonparents` | 0.494 | 0.316 | 0.326 | +0.010 |
| `old_age_own_rate` | 0.863 | 0.709 | 0.666 | -0.043 |
| `old_age_parent_childless_gap` | 0.070 | 0.291 | 0.288 | +0.003 |

### fertility_mid_finance_beta_mid
| Moment | Target | Baseline | Probe | Direction |
|---|---:|---:|---:|---:|
| `tfr` | 1.700 | 1.558 | 1.722 | +0.120 |
| `childless_rate` | 0.150 | 0.351 | 0.267 | +0.084 |
| `mean_age_first_birth` | 26.000 | 34.793 | 33.937 | +0.855 |
| `tfr_gradient` | 0.133 | -0.239 | -0.204 | +0.036 |
| `own_rate` | 0.627 | 0.575 | 0.506 | -0.069 |
| `own_gradient` | 0.170 | -0.082 | -0.061 | +0.021 |
| `young_liquid_wealth_to_income` | 0.600 | 1.215 | 1.209 | +0.005 |
| `center_share_nonparents` | 0.494 | 0.316 | 0.321 | +0.005 |
| `old_age_own_rate` | 0.863 | 0.709 | 0.666 | -0.043 |
| `old_age_parent_childless_gap` | 0.070 | 0.291 | 0.289 | +0.002 |

### center_mild
| Moment | Target | Baseline | Probe | Direction |
|---|---:|---:|---:|---:|
| `tfr` | 1.700 | 1.558 | 1.614 | +0.056 |
| `childless_rate` | 0.150 | 0.351 | 0.332 | +0.018 |
| `mean_age_first_birth` | 26.000 | 34.793 | 34.620 | +0.173 |
| `tfr_gradient` | 0.133 | -0.239 | -0.334 | -0.095 |
| `own_rate` | 0.627 | 0.575 | 0.628 | +0.051 |
| `own_gradient` | 0.170 | -0.082 | -0.168 | -0.086 |
| `young_liquid_wealth_to_income` | 0.600 | 1.215 | 1.091 | +0.123 |
| `center_share_nonparents` | 0.494 | 0.316 | 0.360 | +0.044 |
| `old_age_own_rate` | 0.863 | 0.709 | 0.772 | +0.063 |
| `old_age_parent_childless_gap` | 0.070 | 0.291 | 0.256 | +0.035 |

### alpha_cons_high
| Moment | Target | Baseline | Probe | Direction |
|---|---:|---:|---:|---:|
| `tfr` | 1.700 | 1.558 | 1.476 | -0.082 |
| `childless_rate` | 0.150 | 0.351 | 0.385 | -0.034 |
| `mean_age_first_birth` | 26.000 | 34.793 | 34.808 | -0.015 |
| `tfr_gradient` | 0.133 | -0.239 | -0.224 | +0.015 |
| `own_rate` | 0.627 | 0.575 | 0.642 | +0.037 |
| `own_gradient` | 0.170 | -0.082 | -0.062 | +0.020 |
| `young_liquid_wealth_to_income` | 0.600 | 1.215 | 1.258 | -0.043 |
| `center_share_nonparents` | 0.494 | 0.316 | 0.346 | +0.030 |
| `old_age_own_rate` | 0.863 | 0.709 | 0.718 | +0.009 |
| `old_age_parent_childless_gap` | 0.070 | 0.291 | 0.211 | +0.079 |

### alpha_high_fertility_mid_finance
| Moment | Target | Baseline | Probe | Direction |
|---|---:|---:|---:|---:|
| `tfr` | 1.700 | 1.558 | 1.673 | +0.116 |
| `childless_rate` | 0.150 | 0.351 | 0.310 | +0.041 |
| `mean_age_first_birth` | 26.000 | 34.793 | 34.445 | +0.348 |
| `tfr_gradient` | 0.133 | -0.239 | -0.251 | -0.011 |
| `own_rate` | 0.627 | 0.575 | 0.644 | +0.035 |
| `own_gradient` | 0.170 | -0.082 | 0.046 | +0.128 |
| `young_liquid_wealth_to_income` | 0.600 | 1.215 | 1.322 | -0.108 |
| `center_share_nonparents` | 0.494 | 0.316 | 0.342 | +0.026 |
| `old_age_own_rate` | 0.863 | 0.709 | 0.727 | +0.018 |
| `old_age_parent_childless_gap` | 0.070 | 0.291 | 0.199 | +0.092 |

