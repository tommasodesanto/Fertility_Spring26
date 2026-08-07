# Experimental repaired-E5 policy packet: quota population closure

## Bottom line

This packet is an **experimental E5-repair arm**, not a promoted model. It reports long-run stationary comparisons, not forecasts. The 0.169 outside-origin entrant-share anchor is provisional and still needs a national ACS re-anchor.

The change is deliberately small: the funded baseline now identifies a fixed retention rate and a fixed outside inflow. Policy population moves only because the policy changes mature city-born household formation. The unidentified logit entry response is bypassed and retained only as the final labeled sensitivity row.

That is why the population effects fall from the old roughly 9–10 percent range to roughly 1–3 percent: the former was driven mostly by the external taste scale; the latter is the transparent renewal multiplier (1−0.169)/0.169 = 4.917, modestly damped by housing feedback.

The maturation-repair chain-6 winner is shown in `candidate_screening.csv` but has no corrected policy rows: the same 0.169 anchor implies Rbar = 1.0550, which violates the explicit Rbar <= 1 feasibility guard. It is excluded rather than assigned a fallback.

## Main results

| Arm | Policy | ΔTFR (%) | Births/HH (%) | Household population (%) | Required immigration at fixed population (identical to percent of fertility gap closed) (%) | Price (%) | Rent (%) | Damping decomposition: arithmetic → realized (damping, pp) | Unidentified entry-response sensitivity (lambda=2) (%) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| floor | Rebated 2% tax | 0.378 | 0.334 | 1.695 | 1.742 | -17.534 | -6.060 | 1.739 → 1.681 (0.058) | 9.980 |
| floor | Rebated 2% tax + purchase grant | 0.750 | 0.675 | 3.464 | 3.518 | -17.290 | -5.782 | 3.505 → 3.405 (0.100) | 10.104 |
| tilt | Rebated 2% tax | 0.202 | 0.178 | 0.897 | 0.918 | -17.710 | -6.260 | 0.917 → 0.893 (0.024) | 8.988 |
| tilt | Rebated 2% tax + purchase grant | 0.489 | 0.441 | 2.231 | 2.263 | -17.210 | -5.692 | 2.258 → 2.207 (0.051) | 9.062 |

The merged required-immigration figure is reported as the reduction in required immigration, so it is numerically identical to the percent of the fertility gap closed. Aggregate stationary births = population + births/HH (accounting identity); see `policy_summary.csv`. The difference between the arithmetic renewal response and the realized population response is the housing-feedback damping reported in the detailed tables.

## Closure and interpretation

At baseline, \(\bar R=(1-s^{out})E_0/B_0\) and \(\bar M=s^{out}E_0\). In each counterfactual, \((\bar R,\bar M)\) are fixed and the driver solves \(S E_0=\bar R S B_0(\text{policy})+\bar M\) jointly with the house price and balanced-budget transfer. Households are the population unit. Renter-paid rent is \((q+\delta+\tau_H)p\), so it differs from the price change when the property-tax rate changes.

The baseline retention rate is a derived accounting identity, not a directly estimated behavioral elasticity. The experiment assumes quota-rationed outside inflow and policy-invariant retention. A balanced-growth-path closure remains an explicit design question to reconsider after this meeting; it is not implemented here.

## Acceptance tests and residuals

All 14 encoded checks passed. Quota baselines reproduce S=1 and the 0.169 share to 1e-12; feasible-arm fixed-population CSVs are byte-identical across modes and to the prior packets, while the two chain-6 screening CSVs are byte-identical to each other; quota objects contain no outside value or taste-scale parameter; and all policy market/fiscal residuals are below 2.6e-05.

- `policy_summary.csv`: one wide row per arm and policy.
- `policy_summary_long.csv` and `tables/`: requested advisor-facing metric rows; the unidentified logit sensitivity is last.
- `residuals.csv`: market and fiscal residuals for every baseline and policy under each closure, plus fixed-population decompositions.
- `calibration/`: complete free/fixed parameter and target-fit tables for all three unchanged certified theta vectors.
- `diagnostics/`: the unchanged standard solution plots for each certified arm plus the policy comparison plot.
- `raw_runs/`: deterministic driver outputs and checkpoints from Torch.

## Run record

- Local loose-tolerance smoke: completed before Torch submission (local diagnostic).
- Torch production array(s): 15447628.
- No calibration, target, policy, solver-numerics, ACS, half-life, CE-unit, promoted-model, or paper-LaTeX object was changed.
