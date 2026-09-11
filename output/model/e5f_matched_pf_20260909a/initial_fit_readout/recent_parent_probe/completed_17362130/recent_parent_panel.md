# Recent-parent observation panel: verified collection

All four read-only batches passed: panel jobs 17363528, 17363529 and 17363530, plus joint-candidate job 17363531. Each passed 13 compiled startup tests with zero Bellman or equilibrium solves. All 641 source hashes, all staged input hashes, every selected checkpoint hash and every result/receipt hash were independently checked. The 18 new panel observations and the joint observation each have one pass: repeat equality is **untested**, not failed. The original baseline has two verified identical passes.

The lead accepts the named synchronized post-fertility, uniform-age and residence-proxy observation as the calibration model analogue under the approved plan. This collection certifies implementation and arithmetic only. Original observer packets retain their diagnostic flags; exact annual ACS equivalence is not claimed. No target, weight, objective or architecture was changed, and no further job was launched by this worker.

The baseline ownership gap is **−5.324359 percentage points**; all 19 local panel values lie between **−5.501821 and −5.131456 points**. The unchanged empirical reference is **+16.289551 points**. The joint candidate gives **−6.342477 points**, 1.018119 points below baseline. These are measured values under the accepted named approximation; there is no weighted loss or objective selection here.

| Observation | Selected-birth ownership | Current-empty ownership | Difference |
|---|---:|---:|---:|
| Baseline | 49.898224% | 55.222583% | −5.324359 pp |
| Joint candidate from original smoke 17360699 | 55.203706% | 61.546183% | −6.342477 pp |

The complete 19-row panel and actual nine parameter values are in `recent_parent_panel.csv` and `recent_parent_panel.json`. The joint candidate is separate in `recent_parent_joint.csv/json`. Each row points to its original observation file and checkpoint/parameter-source hashes. Original source identity remains `7e872053`; observer identity remains `70abd4a8`.

## Raw local sensitivities

Each derivative uses the actual original candidate values: D = (f_plus − f_minus)/(u_plus − u_minus), where u = log(parameter), except u_beta = log(−log(beta_annual)). The beta “plus/minus” case labels refer to beta levels; the discount-rate coordinate reverses that ordering. Each panel point separately renormalized fertility to 2.1, so these derivatives follow that normalized path and do not hold psi fixed.

One-sided slopes use the baseline and the actual transformed steps. Relative asymmetry is the absolute difference between those slopes divided by the absolute central slope; no empirical precision or calibration weight enters this calculation. The ranking below is solely by the absolute raw central derivative.

| Parameter | Central derivative | Minus-side slope | Plus-side slope | Relative asymmetry |
|---|---:|---:|---:|---:|
| H0 | 0.092579117 | 0.087841128 | 0.097412829 | 0.103389 |
| chi | 0.040030467 | 0.007918572 | 0.072791131 | 1.620580 |
| kappa_fert | -0.039314781 | -0.040170856 | -0.038441411 | 0.043990 |
| h_P | 0.030961997 | 0.029048309 | 0.032914347 | 0.124864 |
| kappa_fert_continuation | -0.029134482 | -0.030595659 | -0.027643784 | 0.101319 |
| first_birth_fixed_cost | 0.018806270 | 0.018673704 | 0.018941513 | 0.014240 |
| beta_annual | -0.005156584 | -0.003853970 | -0.006433403 | 0.500221 |
| theta0 | -0.002485273 | -0.003987497 | -0.000952700 | 1.221112 |
| theta1 | 0.000639704 | 0.000635369 | 0.000644126 | 0.013688 |

Housing supply scale H0 has the largest raw response. Chi and the first-birth choice-shock scale kappa_fert follow; chi is substantially asymmetric across the two directions. The first-birth fixed cost is a separate coordinate. This table does not establish parameter identification or recommend an objective or parameter change.

`recent_parent_derivatives.csv` records all exact coordinate values, steps, raw outcomes, central and one-sided slopes, signed/absolute/relative asymmetry and source hashes. Each batch's `completed_JOBID/verification.json` preserves its separate scheduler, numerical, test and hash evidence. Observed batch driver times were about 37–38 seconds for six panel observations and 24 seconds for the joint candidate.
