# Income-disciplined recalibration (M5)

The reported winner estimates the 11 clean-frontier parameters plus `theta0`, `theta1`, and `tenure_choice_kappa` against 15 moments; `theta_n=0` is the only external restriction. Income persistence is the literature anchor (annual rho 0.90, innovation s.d. 0.20) and entrant wealth uses the PSID 18-24 distribution. The winner beats the strict M1 `theta0=0`/`tenure_choice_kappa=0` nested seed under the M5 objective. All reported moments use a strict, exactly repeated winner solve.

| Loss | Residual | theta0 | theta1 | theta_n | Eligible chains |
|---:|---:|---:|---:|---:|---:|
| 9.044422 | 1.48e-05 | 0.311765 | 0.397282 | 0.000000 | 7 |

Complete target-fit, parameter-bound, lifecycle, and plot artifacts are adjacent.

Established 12-moment loss: 3.166207 (ceiling 4.50).
All predeclared acceptance criteria pass: False.
Estimated tenure_choice_kappa: 0.010017.
