# V5 Nz=5 Directional Calibration Audit

Verdict: **promising but not calibrated**

This audit uses the isolated V5 benchmark-normalized outside-option closure at
`Nb=30`, `Nz=5`, \(\rho_z=0.95\), unconditional \(\sigma_z=0.35\), and
\(\kappa_E=10^6\). It is not an optimizer. The goal is to test whether the
bad HANK-\(z\) moments are locally movable before spending time on a larger
search.

Baseline reference: `income_mortgage_risk_v5_hank_z_outside_closure_nz5.log`
with loss `311.299772`.

## Main Read

The economy is movable. The best probe,
`alpha_high_fertility_mid_finance`, sets
\(\alpha_c=0.80\), \(\kappa_f=4.0\), and \(\phi=0.90\). It cuts the loss from
`311.30` to `166.72`, improves 16 of 19 target gaps, moves the ownership
gradient from negative to positive, improves the room moments, improves the
family ownership gap, and gets TFR close to target.

This is not a solved calibration. The first-birth age remains far too late,
childlessness is still too high, young liquid wealth is worse, and the
fertility gradient remains wrong-signed. The useful lesson is that the
prototype is not stuck mechanically; it needs a structured calibration block,
and the fixed \(\alpha_c=0.70\) housing-share assumption is probably too tight
for this HANK-\(z\) economy.

## First-Round Probes

| Case | Seconds | Loss | Loss \(\Delta\) | Improved Moments | Read |
|---|---:|---:|---:|---:|---|
| `fertility_utility_max` | 103.20 | 322.89 | +11.59 | 8/19 | Raises TFR modestly, worsens wealth and family gap. |
| `child_cost_low` | 101.25 | 303.29 | -8.01 | 7/19 | Big fertility response, but overshoots TFR and worsens ownership/wealth. |
| `fertility_logit_high` | 99.58 | 265.47 | -45.82 | 11/19 | Strong fertility improvement, better gradients, but wealth rises sharply. |
| `finance_high` | 104.49 | 280.48 | -30.82 | 15/19 | Small but broad ownership improvements; good sign for \(\phi\). |
| `center_pull` | 132.80 | 484.98 | +173.68 | 8/19 | Geography is movable but this step badly overshoots rent ratio and gradients. |
| `beta_low` | 82.34 | 260.14 | -51.16 | 10/19 | Lowers wealth and age at first birth, but destroys ownership. |
| `bequest_child_low` | 98.62 | 299.31 | -11.99 | 7/19 | Mildly improves old parent gap, not enough elsewhere. |
| `combo_soft` | 138.69 | 329.35 | +18.05 | 8/19 | Invalid directionally: overshoots fertility, damages ownership, and gives negative \(M\). |

## Second-Round Probes

| Case | Seconds | Loss | Loss \(\Delta\) | Improved Moments | Read |
|---|---:|---:|---:|---:|---|
| `fertility_logit_mid` | 102.21 | 286.19 | -25.11 | 11/19 | \(\kappa_f=4.0\) gets TFR close but worsens wealth. |
| `fertility_mid_finance_high` | 97.29 | 252.09 | -59.21 | 12/19 | Better than fertility alone; \(\phi=0.90\) helps ownership gaps. |
| `beta_mid_finance_high` | 126.17 | 275.52 | -35.78 | 11/19 | Wealth improves, but ownership falls too much. |
| `fertility_mid_finance_beta_mid` | 121.18 | 253.22 | -58.08 | 13/19 | Better timing/childlessness, but ownership remains damaged. |
| `center_mild` | 166.21 | 312.46 | +1.16 | 10/19 | Hits ownership level and center share, but rent ratio and gradients deteriorate. |
| `alpha_cons_high` | 114.08 | 206.95 | -104.35 | 14/19 | Housing-share lever is powerful; fixes ownership level and room moments but hurts fertility. |
| `alpha_high_fertility_mid_finance` | 107.73 | 166.72 | -144.58 | 16/19 | Best probe: moves ownership gradient positive and improves most non-wealth moments. |

## Best Probe Moment Table

| Moment | Target | Baseline | Best Probe |
|---|---:|---:|---:|
| `tfr` | 1.700 | 1.558 | 1.673 |
| `childless_rate` | 0.150 | 0.351 | 0.310 |
| `mean_age_first_birth` | 26.000 | 34.793 | 34.445 |
| `tfr_gradient` | 0.133 | -0.239 | -0.251 |
| `own_rate` | 0.627 | 0.575 | 0.644 |
| `own_gradient` | 0.170 | -0.082 | 0.046 |
| `own_family_gap` | 0.110 | 0.383 | 0.184 |
| `prime_childless_renter_median_rooms` | 4.000 | 5.767 | 5.092 |
| `prime_childless_owner_median_rooms` | 6.000 | 6.800 | 5.400 |
| `housing_increment_0to1` | 0.664 | 0.591 | 0.693 |
| `housing_increment_1to2` | 0.566 | 1.492 | 1.375 |
| `young_liquid_wealth_to_income` | 0.600 | 1.215 | 1.322 |
| `center_share_nonparents` | 0.494 | 0.316 | 0.342 |
| `center_share_newparents` | 0.416 | 0.393 | 0.418 |
| `migration_rate` | 0.032 | 0.034 | 0.034 |
| `old_age_own_rate` | 0.863 | 0.709 | 0.727 |
| `old_age_parent_childless_gap` | 0.070 | 0.291 | 0.199 |
| `inv_pop_share_C` | 0.450 | 0.409 | 0.431 |
| `inv_rent_ratio_C_over_P` | 1.140 | 1.196 | 1.182 |

Best-probe closure diagnostics remain normalized:
\(S=1\), \(q^E=0.9\), outside probability `0.1`, residual outside-born flow
\(M=0.00485336\), outside value `-1504146.387742`, and final GE error
`4.10132e-4`.

## Calibration Implications

- Use `Nz=5` for the exploratory search and keep `Nz=7` as validation.
- Add \(\alpha_c\) or an equivalent housing-share/room-demand parameter to the
  active search block. The fixed `alpha_cons=0.70` assumption is too restrictive
  once idiosyncratic \(z\) risk is active.
- Keep \(\kappa_f\) and \(\phi\) in the near-term search. They jointly improve
  fertility and ownership gradients.
- Do not lower \(\beta\) aggressively without an offsetting tenure channel:
  it helps wealth and timing but damages ownership and old-age ownership.
- Treat geography separately. Center amenity/rent shifter moves are powerful,
  but even the mild probe worsened the rent ratio and ownership gradient.
- The remaining hard failures are first-birth timing, childlessness, fertility
  gradient, and liquid wealth. Those likely need age-profile fertility costs or
  income/asset-return discipline, not just scalar fertility taste changes.

## Outputs

- First round: `v5_nz5_directional_calibration.csv`,
  `v5_nz5_directional_calibration_moments.csv`, and
  `v5_nz5_directional_calibration.json`.
- Second round: `v5_nz5_directional_calibration_round2.csv`,
  `v5_nz5_directional_calibration_round2_moments.csv`, and
  `v5_nz5_directional_calibration_round2.json`.
