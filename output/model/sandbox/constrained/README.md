# Who faces a binding down payment, and does scarcer space change it?

Two steady states at the same retained parameters and the same fixed child benefit level are compared: a baseline with down payment share 1 - financed share = 0.20, and a no-down-payment benchmark (financed share 1.0). A second pair repeats the comparison after lowering the supply scale H0 until mean occupied rooms falls to the data target. Rows are model age (18 to 62) by children at home (m = 0, 1, 2, 3 or more).

## Definitions

- Renter share: realized renter mass / cell mass, from the stationary distribution.
- Renter at cap: share of renters whose planned rooms sit at the rental cap (6 rooms, within 0.0001); renter rooms differ across previous-tenure states, so this uses pre-choice mass times the rent probability at each state (up to the small smoothing share).
- Owner-only sizes: share of owners in rungs at or above 8 rooms (owner grid [2.0, 4.0, 6.0, 8.0, 10.0]).
- Binds share: cell mass whose modal (tenure, rooms) choice differs across the two solutions at the same pre-choice state, with a 0.05 room tolerance; pre-choice weights use the baseline distribution (tenure axis as proxy, exact up to the small smoothing share).
- Median renter wealth: median liquid wealth over realized renters in the cell, also scaled by the down payment on a 8-room unit at the solved price (baseline down payment 1.265; no-down-payment benchmark has none).
- Attempt rate: first-birth attempt rate over not-yet-parent mass (one-shot setup: all births flow from n = 0, so non-zero-m rows show n/a; outside the fertile decision ages 18-42 the menu is unset, so n/a there too).
- First births from cap: age-level share of first-birth children coming from renters at the cap (birth-weighted by the model's own expected-children menu).

## Part A: baseline (down payment kept)

| Age | m | Cell mass | Renters | Renters at cap | Owners 8+ rooms | Binds | Median renter wealth | Median / down payment | Attempt (n=0) | First births from cap |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 18 | 0 | 0.0484 | 0.978 | 0.090 | 0.856 | 0.065 | 0.000 | 0.00 | 0.186 | 0.216 |
| 18 | 1 | 0.0133 | 0.901 | 0.667 | 0.950 | 0.162 | 0.163 | 0.13 | n/a | 0.216 |
| 18 | 2 | 0.0000 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | 0.216 |
| 18 | 3+ | 0.0000 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | 0.216 |
| 22 | 0 | 0.0391 | 0.950 | 0.038 | 0.852 | 0.066 | 0.163 | 0.13 | 0.200 | 0.078 |
| 22 | 1 | 0.0192 | 0.771 | 0.590 | 0.964 | 0.078 | 0.302 | 0.24 | n/a | 0.078 |
| 22 | 2 | 0.0034 | 0.690 | 0.666 | 0.971 | 0.128 | 0.163 | 0.13 | n/a | 0.078 |
| 22 | 3+ | 0.0000 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | 0.078 |
| 26 | 0 | 0.0324 | 0.938 | 0.030 | 0.809 | 0.044 | 0.000 | 0.00 | 0.207 | 0.050 |
| 26 | 1 | 0.0216 | 0.736 | 0.714 | 0.950 | 0.068 | 0.163 | 0.13 | n/a | 0.050 |
| 26 | 2 | 0.0073 | 0.642 | 0.788 | 0.965 | 0.100 | 0.000 | 0.00 | n/a | 0.050 |
| 26 | 3+ | 0.0004 | 0.585 | 0.836 | 0.975 | 0.111 | 0.000 | 0.00 | n/a | 0.050 |
| 30 | 0 | 0.0281 | 0.911 | 0.024 | 0.766 | 0.037 | 0.163 | 0.13 | 0.186 | 0.030 |
| 30 | 1 | 0.0218 | 0.644 | 0.699 | 0.938 | 0.042 | 0.302 | 0.24 | n/a | 0.030 |
| 30 | 2 | 0.0106 | 0.544 | 0.761 | 0.958 | 0.049 | 0.302 | 0.24 | n/a | 0.030 |
| 30 | 3+ | 0.0012 | 0.473 | 0.799 | 0.971 | 0.052 | 0.163 | 0.13 | n/a | 0.030 |
| 34 | 0 | 0.0259 | 0.871 | 0.024 | 0.721 | 0.059 | 0.163 | 0.13 | 0.157 | 0.034 |
| 34 | 1 | 0.0212 | 0.565 | 0.681 | 0.922 | 0.031 | 0.442 | 0.35 | n/a | 0.034 |
| 34 | 2 | 0.0126 | 0.479 | 0.733 | 0.947 | 0.040 | 0.302 | 0.24 | n/a | 0.034 |
| 34 | 3+ | 0.0021 | 0.423 | 0.763 | 0.962 | 0.058 | 0.302 | 0.24 | n/a | 0.034 |
| 38 | 0 | 0.0252 | 0.807 | 0.023 | 0.649 | 0.063 | 0.442 | 0.35 | 0.133 | 0.037 |
| 38 | 1 | 0.0208 | 0.485 | 0.659 | 0.904 | 0.014 | 0.581 | 0.46 | n/a | 0.037 |
| 38 | 2 | 0.0130 | 0.407 | 0.694 | 0.935 | 0.016 | 0.581 | 0.46 | n/a | 0.037 |
| 38 | 3+ | 0.0027 | 0.368 | 0.724 | 0.951 | 0.039 | 0.442 | 0.35 | n/a | 0.037 |
| 42 | 0 | 0.0262 | 0.696 | 0.022 | 0.534 | 0.104 | 0.581 | 0.46 | 0.142 | 0.024 |
| 42 | 1 | 0.0212 | 0.426 | 0.636 | 0.888 | 0.010 | 0.860 | 0.68 | n/a | 0.024 |
| 42 | 2 | 0.0120 | 0.358 | 0.674 | 0.924 | 0.008 | 0.721 | 0.57 | n/a | 0.024 |
| 42 | 3+ | 0.0024 | 0.341 | 0.681 | 0.938 | 0.067 | 0.581 | 0.46 | n/a | 0.024 |
| 46 | 0 | 0.0315 | 0.583 | 0.024 | 0.510 | 0.108 | 0.581 | 0.46 | n/a | n/a |
| 46 | 1 | 0.0209 | 0.365 | 0.626 | 0.888 | 0.006 | 0.860 | 0.68 | n/a | n/a |
| 46 | 2 | 0.0082 | 0.310 | 0.634 | 0.922 | 0.009 | 0.721 | 0.57 | n/a | n/a |
| 46 | 3+ | 0.0011 | 0.296 | 0.637 | 0.936 | 0.109 | 0.721 | 0.57 | n/a | n/a |
| 50 | 0 | 0.0366 | 0.492 | 0.025 | 0.493 | 0.124 | 0.721 | 0.57 | n/a | n/a |
| 50 | 1 | 0.0192 | 0.322 | 0.602 | 0.886 | 0.005 | 0.860 | 0.68 | n/a | n/a |
| 50 | 2 | 0.0054 | 0.275 | 0.599 | 0.921 | 0.054 | 0.721 | 0.57 | n/a | n/a |
| 50 | 3+ | 0.0005 | 0.263 | 0.595 | 0.935 | 0.131 | 0.721 | 0.57 | n/a | n/a |
| 54 | 0 | 0.0411 | 0.413 | 0.025 | 0.472 | 0.114 | 0.721 | 0.57 | n/a | n/a |
| 54 | 1 | 0.0169 | 0.290 | 0.584 | 0.882 | 0.019 | 0.860 | 0.68 | n/a | n/a |
| 54 | 2 | 0.0035 | 0.251 | 0.577 | 0.920 | 0.080 | 0.860 | 0.68 | n/a | n/a |
| 54 | 3+ | 0.0003 | 0.241 | 0.569 | 0.934 | 0.155 | 0.721 | 0.57 | n/a | n/a |
| 58 | 0 | 0.0450 | 0.341 | 0.022 | 0.451 | 0.120 | 0.581 | 0.46 | n/a | n/a |
| 58 | 1 | 0.0144 | 0.268 | 0.533 | 0.875 | 0.048 | 1.000 | 0.79 | n/a | n/a |
| 58 | 2 | 0.0022 | 0.237 | 0.525 | 0.915 | 0.119 | 0.860 | 0.68 | n/a | n/a |
| 58 | 3+ | 0.0001 | 0.230 | 0.560 | 0.931 | 0.181 | 0.860 | 0.68 | n/a | n/a |
| 62 | 0 | 0.0483 | 0.283 | 0.019 | 0.433 | 0.099 | 0.302 | 0.24 | n/a | n/a |
| 62 | 1 | 0.0120 | 0.260 | 0.496 | 0.870 | 0.053 | 1.000 | 0.79 | n/a | n/a |
| 62 | 2 | 0.0014 | 0.236 | 0.486 | 0.912 | 0.173 | 1.000 | 0.79 | n/a | n/a |
| 62 | 3+ | 0.0001 | 0.232 | 0.509 | 0.930 | 0.195 | 1.000 | 0.79 | n/a | n/a |

## Part A: no-down-payment benchmark (same rows where they change)

| Age | m | Renters | Renters at cap | Owners 8+ rooms | Median renter wealth | Attempt (n=0) |
|---:|---|---:|---:|---:|---:|---:|
| 18 | 0 | 0.917 | 0.053 | 0.590 | 0.000 | 0.185 |
| 18 | 1 | 0.732 | 0.590 | 0.923 | 0.000 | n/a |
| 18 | 2 | n/a | n/a | n/a | n/a | n/a |
| 18 | 3+ | n/a | n/a | n/a | n/a | n/a |
| 22 | 0 | 0.882 | 0.022 | 0.484 | 0.163 | 0.198 |
| 22 | 1 | 0.673 | 0.499 | 0.879 | 0.163 | n/a |
| 22 | 2 | 0.535 | 0.561 | 0.880 | 0.163 | n/a |
| 22 | 3+ | n/a | n/a | n/a | n/a | n/a |
| 26 | 0 | 0.860 | 0.016 | 0.456 | 0.000 | 0.205 |
| 26 | 1 | 0.633 | 0.564 | 0.846 | 0.163 | n/a |
| 26 | 2 | 0.491 | 0.704 | 0.856 | 0.000 | n/a |
| 26 | 3+ | 0.393 | 0.757 | 0.871 | 0.000 | n/a |
| 30 | 0 | 0.826 | 0.012 | 0.449 | 0.163 | 0.184 |
| 30 | 1 | 0.561 | 0.630 | 0.819 | 0.302 | n/a |
| 30 | 2 | 0.429 | 0.653 | 0.840 | 0.302 | n/a |
| 30 | 3+ | 0.335 | 0.714 | 0.863 | 0.163 | n/a |
| 34 | 0 | 0.782 | 0.015 | 0.450 | 0.163 | 0.155 |
| 34 | 1 | 0.485 | 0.614 | 0.795 | 0.302 | n/a |
| 34 | 2 | 0.371 | 0.625 | 0.818 | 0.163 | n/a |
| 34 | 3+ | 0.295 | 0.648 | 0.841 | 0.163 | n/a |
| 38 | 0 | 0.710 | 0.018 | 0.432 | 0.302 | 0.131 |
| 38 | 1 | 0.413 | 0.596 | 0.774 | 0.581 | n/a |
| 38 | 2 | 0.312 | 0.574 | 0.795 | 0.442 | n/a |
| 38 | 3+ | 0.256 | 0.562 | 0.809 | 0.302 | n/a |
| 42 | 0 | 0.574 | 0.019 | 0.374 | 0.442 | 0.140 |
| 42 | 1 | 0.354 | 0.564 | 0.757 | 0.721 | n/a |
| 42 | 2 | 0.262 | 0.517 | 0.774 | 0.442 | n/a |
| 42 | 3+ | 0.222 | 0.488 | 0.771 | 0.302 | n/a |
| 46 | 0 | 0.459 | 0.023 | 0.383 | 0.442 | n/a |
| 46 | 1 | 0.289 | 0.522 | 0.751 | 0.721 | n/a |
| 46 | 2 | 0.206 | 0.433 | 0.758 | 0.442 | n/a |
| 46 | 3+ | 0.167 | 0.364 | 0.752 | 0.302 | n/a |
| 50 | 0 | 0.372 | 0.026 | 0.388 | 0.581 | n/a |
| 50 | 1 | 0.239 | 0.479 | 0.744 | 0.581 | n/a |
| 50 | 2 | 0.162 | 0.346 | 0.745 | 0.442 | n/a |
| 50 | 3+ | 0.126 | 0.241 | 0.735 | 0.163 | n/a |
| 54 | 0 | 0.299 | 0.027 | 0.389 | 0.581 | n/a |
| 54 | 1 | 0.200 | 0.435 | 0.737 | 0.581 | n/a |
| 54 | 2 | 0.127 | 0.259 | 0.733 | 0.442 | n/a |
| 54 | 3+ | 0.092 | 0.127 | 0.720 | 0.163 | n/a |
| 58 | 0 | 0.235 | 0.026 | 0.387 | 0.442 | n/a |
| 58 | 1 | 0.170 | 0.363 | 0.730 | 0.721 | n/a |
| 58 | 2 | 0.100 | 0.144 | 0.720 | 0.442 | n/a |
| 58 | 3+ | 0.065 | 0.040 | 0.704 | 0.163 | n/a |
| 62 | 0 | 0.181 | 0.022 | 0.385 | 0.302 | n/a |
| 62 | 1 | 0.150 | 0.267 | 0.720 | 0.581 | n/a |
| 62 | 2 | 0.080 | 0.042 | 0.704 | 0.302 | n/a |
| 62 | 3+ | 0.047 | 0.011 | 0.687 | 0.163 | n/a |

## Part B: scarcer space

| Solution | Completed fertility | Childless share | Mean first-birth age | First births 30+ | Ownership 30-55 | Mean rooms | First-birth rooms response | Rooms gap 3+ vs 1-2 | Price | Binds 26-38, m 0-1 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline_dp | 1.8717 | 0.2386 | 26.961 | 0.2847 | 0.4585 | 5.7180 | 1.0719 | 0.3080 | 0.7906 | 0.0454 |
| baseline_nodp | 1.8658 | 0.2410 | 26.973 | 0.2857 | 0.5569 | 5.8233 | 1.1493 | 0.3486 | 0.7989 | 0.0454 |
| scarce_dp | 1.8401 | 0.2467 | 27.097 | 0.2929 | 0.4369 | 5.5625 | 1.0741 | 0.3371 | 0.8234 | 0.0279 |
| scarce_nodp | 1.8353 | 0.2486 | 27.107 | 0.2938 | 0.5317 | 5.6550 | 1.1648 | 0.3534 | 0.8312 | 0.0279 |

## Reading (ten lines)

1. At family-forming ages (26-38, m 0-1) the binds share is 0.045: removing the down payment changes the modal housing choice for about that fraction of mass.
2. Renting at those ages averages 0.766, and 0.313 of renters sit at the 6-room cap, so the cap rather than the down payment is the visible wall.
3. Fully 0.823 of owners at those ages hold 8 or more rooms: owners who buy, buy big; the down payment screens entry rather than size.
4. Median renter wealth is only 0.21 times the 8-room down payment (1.265): the median renter cannot cover it, yet most still would not buy without it -- wealth is short but renting is preferred.
5. The n=0 attempt rate at 26-38 averages 0.173, and only 0.038 of first births come from capped renters: constrained renters contribute few births.
6. Ownership 30-55 rises from 0.459 to 0.557 without the down payment, while completed fertility moves only from 1.872 to 1.866: tenure responds, births barely do.
7. Scarcity (H0 8.11 to 7.35) cuts mean rooms from 5.718 to 5.563, inside 0.05 of the 5.56 target, with the price up from 0.791 to 0.823.
8. Yet the binds share at 26-38 falls from 0.045 to 0.028 under scarcity: higher prices push more young households into renting, but the marginal own-vs-rent choice moves less, not more.
9. Fertility edges down under scarcity (completed fertility 1.872 to 1.840, first-birth age 26.96 to 27.10) while the no-down-payment gaps stay put: space costs delay births slightly without making the down payment pivotal.
10. Bottom line: the down payment binds for a small minority of family-forming households at the baseline and does not start to bind once space is scarce; it moves tenure, not births.

