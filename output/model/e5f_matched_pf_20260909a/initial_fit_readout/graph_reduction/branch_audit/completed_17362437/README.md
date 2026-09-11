# Completed two-node branch audit: no local saving error found

Job **17362437 completed, exit 0:0**, with 17 scalar branch optimizations in 6.496 seconds (13 seconds of job wall time; one CPU; reported peak RSS 694812K). The collected audit, driver and input contract match remote SHA256 hashes; the driver and contract also match the reviewed local files. The completed driver verified all 625 source pins before loading the pinned checkpoint. See `collection_verification.json`, `sacct.txt`, `stdout.txt`, and the unchanged `audit.json`.

This is the **initial stationary `new_balanced`, fixed-preference diagnostic** from source c6dd3508. Age 30 is a household age, not a 2023 calendar observation. No full Bellman or equilibrium solve was performed. The inspector used the exact saved pension 2.0463613896121218, child-preference intercept 0.2900515293650047, beta 0.981239759614208, alpha 0.733 and owner service premium 1.0537270178404272. The inspected no-child state has zero current consumption/housing floors and equivalence scale one. Its current survival probability is one; its saved continuation retains subsequent survival and bequest incentives.

## Conclusion for this conditional slice

**Saved saving choices are optimal under the inspected saved continuation to floating-point precision. The ownership-probability decrease survives direct off-grid branch optimization. Current-state owner-value interpolation slightly dampens the decrease, rather than causing it.** This supports a relative-value explanation for the two inspected nodes. It does not prove that every policy is regular or that the saved continuation/grid is accurate everywhere.

The state is renter/location0/parity0/child-state0, age 30, combined income index 8. The two wealth nodes are 0.8604651163 and 1. At both nodes, the affordable menu contains renting and 2-, 4- and 6-room ownership. The six-room and eight-room down-payment thresholds are 0.7599943224 and 1.0133257632, so the pair does not cross either threshold. Conditional rental housing is six rooms at both nodes.

## Actual formula versus the off-grid diagnostic

The production formula interpolates conditional owner values after subtracting the house purchase price from current wealth. The direct calculation instead reoptimizes each feasible owner branch at that exact post-purchase wealth, using the **same** saved continuation. It is a diagnostic, not a replacement for the production formula.

| Quantity | Wealth 0.8604651163 | Wealth 1 | Change |
|---|---:|---:|---:|
| Renter value | 2.342719654225 | 2.368434043784 | +0.025714389559 |
| Six-room value, saved interpolation | 2.337941349223 | 2.363319611301 | +0.025378262078 |
| Six-room value, direct off-grid optimum | 2.337960595705 | 2.363338285949 | +0.025377690244 |
| Six-room minus renter, saved formula | -0.004778305002 | -0.005114432483 | -0.000336127481 |
| Six-room minus renter, direct diagnostic | -0.004759058520 | -0.005095757835 | -0.000336699315 |

The renter value rises faster: its secant slope is 0.184286459 per wealth unit, versus 0.181877545 for the interpolated six-room branch and 0.181873447 for the direct branch. The six-room disadvantage therefore widens in both comparisons.

Interpolation lowers six-room value by 1.924648205e-05 and 1.867464779e-05 at the two nodes. These are 0.4028% and 0.3651% of the respective absolute six-room value disadvantage. Its change across the pair is only -5.71834267e-07: removing current-state interpolation makes the decline in relative owner value about 0.1698% larger. It does not reverse the sign.

Using the actual saved probabilities, total ownership falls **27.966025% to 26.642793%**, a decline of **1.323232 percentage points**. If one merely feeds the direct off-grid values of all three feasible owner products into the same softmax, the diagnostic probabilities are **28.043807% to 26.716054%**, a decline of **1.327752 pp**. This arithmetic illustration changes neither stored probabilities nor population/market outcomes. Current-state interpolation changes this two-node probability decrease by only about 0.004520 pp.

## Optimization, allocation and probability checks

- All **11 saved on-grid branch actions** (two renters and three bracket nodes for each of the 2-, 4- and 6-room owner products) match independently reconstructed scalar optima. Maximum best-minus-saved objective gain: **8.882e-16**. Maximum saving difference: **8.882e-16**.
- All **11 consumption gaps are exactly zero** at stored precision; both renter housing reconstruction gaps are exactly zero.
- All **20 feasible crossed saving actions** perform worse in the target state's own objective. The least-negative gain is **-0.0013458585**; the most-negative is **-0.0083426167**. No feasible crossed action reveals an improvement.
- The six-room/renter log-odds identity reproduces to **1.715e-10 utility units**. Maximum full-softmax probability discrepancy is **2.627e-08**, consistent with saved float32 probabilities.
- The six direct off-grid optimizations are reported separately. They do not overwrite the production value interpolation or imply a re-solved equilibrium.

## Population relevance and limits

Post-fertility/pre-tenure household masses are **0.000100590479779** and **7.44851632209e-05** at the lower and upper nodes. With global household mass approximately one, these represent **0.01005905%** and **0.00744852%** of global households, or **0.01750756% together**. The lower node is 0.732380% of the selected age-30 childless-renter slice. These are exposure masses, not an aggregate ownership or welfare change.

The evidence rules out a suboptimal saving action and a probability-reconstruction error **for this saved conditional slice**. The fall remains in independently optimized branches using the saved continuation; it is consistent with renters' relative value rising faster between owner-product entry thresholds. The audit does not identify every deeper component of that relative-value slope, prove a global economic monotonicity theorem, or certify continuation-grid accuracy. Other occupied reversals and late-life income-gradient crossings remain separate questions. No further computation is needed to answer this specific two-node optimizer-versus-current-state-interpolation question.

Full arithmetic and source-output hash are in `interpretation.json`. No code, source, parameter, target, numerical gate, job or original graph changed during this collection and readout.
