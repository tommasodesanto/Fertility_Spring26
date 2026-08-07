# Population-closure arithmetic

This folder contains deterministic demographic arithmetic only; it does not call a model solver.  The quota anchor is provisional: $s=0.169$ is the current across-CBSA outside-origin entrant share.  All percentages below use the supplied mature entrant-flow changes, not the births-per-household changes.

## Inputs and accounting

The supplied baseline identity holds exactly: $E_0=\bar R B_0+\bar M$, with $\rho=\bar R B_0/E_0=0.831$.  The arithmetic multiplier is $(1-s)/s$; at $s=0.169$, its raw products are $1.7385602\%$ (tax) and $3.5052957\%$ (tax+grant).  The packet's requested four-decimal displays are **1.7385%** and **3.5053%**; because the first differs from conventional rounding of the raw input product (which is 1.7386%), `*_packet_4dp` records that supplied display separately from the raw calculation.

## Artifact 1: convergence after a permanent birth shift

Method A is an exact-$J=17$ age-structured transition.  A household cohort enters at model age 18, survives precisely 17 four-year periods, and the child-at-home stock matures at hazard $2/9$.  Baseline household population is normalized to one, so every supplied flow is rescaled by $[(1/J)/E_0]$ before the simulation; this makes the baseline stationary by construction (maximum one-period state drift: $6.939\times10^{-18}$).  This is a **bound**, since it ignores early mortality even though $E_0>1/17$ in the supplied model.

| Method A bound | $d\ln B$ | New population gap | Years to 25% | Years to 50% | Years to 75% |
|---|---:|---:|---:|---:|---:|
| tax | 0.35357% | 1.769320886% | 108 | 232 | 444 |
| tax+grant | 0.71287% | 3.632630107% | 112 | 236 | 452 |

Method B is the generational closed form.  With generation length $g$, the remaining population gap is multiplied by $\rho$ each generation, so $t_p=g\ln[1/(1-p)]/\ln(1/\rho)$.

| $g$ (years) | Closure | Years | Half-life (years) |
|---:|---:|---:|---:|
| 26 | 25% | 40.403588512 | 97.349248158 |
| 26 | 50% | 97.349248158 | 97.349248158 |
| 26 | 75% | 194.698496315 | 97.349248158 |
| 30 | 25% | 46.619525206 | 112.326055567 |
| 30 | 50% | 112.326055567 | 112.326055567 |
| 30 | 75% | 224.652111133 | 112.326055567 |

**Prominent consistency flag.** Method A's 50% closure time is 232--236 years, while Method B gives 97.35 years ($g=26$) or 112.33 years ($g=30$).  Thus the exact-$J$ bound is about 2.1--2.4 times the generational half-life, slightly exceeding the requested factor-of-two consistency screen; this is not smoothed over and should be resolved before treating either as a forecast.  Both methods nevertheless put the adjustment on the order of decades to a century-plus, not an immediate stationary response.

## Artifact 2: anchor sensitivity

`anchor_sensitivity.csv` evaluates the requested grid, $s=0.10,0.105,\ldots,0.30$.  The shaded region in `anchor_sensitivity.png` is $s<0.1432$, infeasible for both certified arms because $\bar R>1$ (individual thresholds: floor 0.1426; tilt 0.1432).

| $s$ | Tax arithmetic $d\ln S$ | Tax+grant arithmetic $d\ln S$ |
|---:|---:|---:|
| 0.1000 | 3.1821% | 6.4158% |
| 0.1200 | 2.5928% | 5.2277% |
| 0.1432 | 2.1155% | 4.2653% |
| 0.1690 | 1.7385% | 3.5053% |
| 0.2000 | 1.4143% | 2.8515% |
| 0.2500 | 1.0607% | 2.1386% |
| 0.3000 | 0.8250% | 1.6634% |

The multiplier makes the closure effect sharply sensitive to the outside-origin anchor, especially at low $s$.  At the provisional anchor, the requested packet displays are reproduced as 1.7385% and 3.5053%.  The speed comparison is a diagnostic bound versus a generational approximation, not a population forecast or a model solve.

## Files and reproduction

- `build_population_closure_arithmetic.py` — self-contained generator.
- `method_a_closure_times.csv` and `method_a_transition_paths.csv` — exact-$J$ bound results.
- `method_b_closure_times.csv` — generational closed form.
- `anchor_sensitivity.csv` and `anchor_sensitivity_selected.csv` — full requested grid and the table above.
- `anchor_sensitivity.png` — 240-dpi figure.

Run with the project environment:

```sh
MPLCONFIGDIR=/private/tmp/population_closure_mpl ../../../code/model/.venv/bin/python build_population_closure_arithmetic.py
```
