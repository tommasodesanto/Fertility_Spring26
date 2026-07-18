# Transfer-floor probe (2026-07-18)

Question: can realistic income risk survive once a *measured* means-tested floor replaces the oversized Stone-Geary intercept as the bottom-of-distribution backstop — with no new heterogeneity? Fixed-theta probe at the M5 winner; nothing recalibrated.

Design: `T = min(max(0, G(n,s) - x), G(n,s))` with **debt-blind means test** `x = y + R*max(b,0)` (v2; v1 tested on cash-on-hand `R*b+y` and made every mortgaged owner eligible — 79% receipt, floor became a mortgage subsidy; see `tmp/transfer_floor_probe_spec_20260718.md` + v2 addendum, archived in this folder). G(n,s) = G0 + Gn*n with children present. Implemented default-off in the markov Bellman kernels + feasibility census; solver kernels reject nothing new when off: **cell 0 reproduces gate0 bitwise** (loss 8.99896940330704, residual 9.324104321982198e-06; local platform levels per the same-run-delta rule — Torch remains canonical for reported levels).

Units: annual figures are fractions of mean annual household income; period = annual x 4. M5 c_bar_0 = 0.3149 annual; probe cells lower it by hand (no refit).

| cell | sigma_ann | c_bar_0 ann | G0 ann (+Gn/kid) | status | loss | key result |
|---|---|---|---|---|---:|---|
| 0 base_repro | 0.0645 | 0.3149 | none | ok | 9.00 | = gate0 bitwise |
| 1 floor_only | 0.0645 | 0.3149 | 0.13 (+0.10) | ok | 9.00 | **bitwise = baseline; zero receipt** — measured floor never touched at M5 |
| 2 lowc_nofloor | 0.20 | 0.10 | none | infeasible | - | dead at 22 (entry-debt renters, bottom z): low c_bar_0 alone is NOT enough |
| 3 s12_lowc | 0.12 | 0.10 | 0.13 (+0.10) | ok | 639.89 | feasible; thin receipt |
| 4 s20_lowc | 0.20 | 0.10 | 0.13 (+0.10) | ok | 452.52 | **flagship**: feasible, receipt 0.9%, outlays 5bp, p90/50 1.75->3.07 |
| 5 s12_midc | 0.12 | 0.14 | 0.18 (+0.10) | ok | 507.88 | feasible; thin receipt |
| 6 s20_midc | 0.20 | 0.14 | 0.18 (+0.10) | infeasible | - | dead at 26 (deep-debt bottom-z, dynamic amortization infeasibility): corridor needs low c_bar_0 |
| 7 gate_check_smallG | 0.20 | 0.3149 | 0.13 (+0.10) | infeasible | - | canary: too-small floor cannot mask infeasibility (dead mass 0.006) |

## Findings

1. **Backstop property (cell 1).** Adding the measured US-scale floor (13% of mean income + 10%/child, debt-blind) to the *current* calibration changes nothing: zero receipt, zero outlays, bitwise-identical equilibrium. The current calibration never touches the safety net.
2. **Necessity and the feasible corridor (cells 2, 4, 6).** At sigma=0.20 the floor is necessary (cell 2 dies at age 22 even with c_bar_0=0.10) and sufficient only with the low bundle: c_bar_0=0.10 annual works (cell 4); c_bar_0=0.14 still dies (cell 6) through deep-entry-debt renters who cannot meet the debt-amortization schedule out of floored income — the dynamic version of the feasibility conflict (points to the M6 forbearance margin note).
3. **The wealth tail (cells 3-5).** Honest persistent risk builds the late-life tail with no added heterogeneity: estate p90/50 1.751 -> 2.18 (sigma=0.12) / **3.065 (sigma=0.20; target 3.45)** at fixed parameters, through the Rouwenhorst grid's own upper states — consistent with the July-17 type-probe direction, but from the process itself.
4. **The costs a refit must absorb (cell 4).** Young liquid wealth 0.323 -> 1.328 (target 0.179): honest persistent risk quadruples young precautionary saving, and the thin floor rightly does not undo it — this is the binding identification tension (note beta_annual >= 0.94 external floor). Childlessness 0.189 -> 0.277 (risk postpones/suppresses fertility); TFR 2.035 -> 2.157; own_rate_2534 overshoots from below (0.276 -> 0.434 vs 0.341). The 9435-weight 65-75 nonhousing>=1x-own-income share falls 0.610 -> 0.415 partly mechanically (higher income dispersion raises denominators for high-z retirees) — moment-definition sensitivity worth flagging before any refit.
5. **Fiscal size (cell 4).** Receipt 0.91% of households, outlays 0.054% of gross income (~5bp funding levy if closed with a proportional tax); receipt 0.9142% of households (renter 0.5531%, owner 0.3611%; ages 18-33 0.5620%, 66+ 0.0000%); outlays 0.0539% of gross income.

## Cell 0 base_repro — full target fit (loss 8.9990, residual 9.32e-06)

receipt 0.0000% of households (renter 0.0000%, owner 0.0000%; ages 18-33 0.0000%, 66+ 0.0000%); outlays 0.0000% of gross income; entry_censored_share 0.0000

| moment | target | model | gap | weight | loss contrib |
|---|---:|---:|---:|---:|---:|
| tfr | 1.9180 | 2.0346 | +0.1166 | 20.0 | 0.272 |
| childless_rate | 0.1880 | 0.1888 | +0.0008 | 20.0 | 0.000 |
| own_rate | 0.5755 | 0.6583 | +0.0828 | 100.0 | 0.686 |
| own_family_gap | 0.1677 | 0.2177 | +0.0500 | 45.0 | 0.113 |
| housing_increment_0to1 | 0.6644 | 0.6808 | +0.0164 | 14.0 | 0.004 |
| prime30_55_childless_renter_mean_rooms | 3.8053 | 3.9681 | +0.1628 | 6.0 | 0.159 |
| prime30_55_childless_owner_share_rooms_ge6 | 0.5961 | 0.7584 | +0.1622 | 25.0 | 0.658 |
| young_childless_renter_liquid_wealth_to_annual_gross_income_2535 | 0.1792 | 0.3226 | +0.1434 | 12.0 | 0.247 |
| prime30_55_childless_owner_minus_renter_mean_rooms | 2.4188 | 2.3132 | -0.1055 | 12.0 | 0.134 |
| own_rate_2534 | 0.3412 | 0.2760 | -0.0652 | 80.0 | 0.340 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.3677 | 0.6017 | +0.2340 | 8.0 | 0.438 |
| old_total_estate_wealth_to_annual_income_median_7684 | 6.5013 | 6.4305 | -0.0708 | 18.6 | 0.093 |
| old_nonhousing_ge_1x_income_share_6575 | 0.6083 | 0.6100 | +0.0017 | 9435.2 | 0.026 |
| old_age_own_rate | 0.7643 | 0.9540 | +0.1897 | 160.0 | 5.760 |
| aggregate_mean_occupied_rooms_18_85 | 5.7800 | 5.6721 | -0.1078 | 6.0 | 0.070 |

Untargeted: estate p90/p50 (76-84) 1.751; owner_neg_liquid_share_2534 0.569.

## Cell 1 floor_only — full target fit (loss 8.9990, residual 9.32e-06)

receipt 0.0000% of households (renter 0.0000%, owner 0.0000%; ages 18-33 0.0000%, 66+ 0.0000%); outlays 0.0000% of gross income; entry_censored_share 0.0000

| moment | target | model | gap | weight | loss contrib |
|---|---:|---:|---:|---:|---:|
| tfr | 1.9180 | 2.0346 | +0.1166 | 20.0 | 0.272 |
| childless_rate | 0.1880 | 0.1888 | +0.0008 | 20.0 | 0.000 |
| own_rate | 0.5755 | 0.6583 | +0.0828 | 100.0 | 0.686 |
| own_family_gap | 0.1677 | 0.2177 | +0.0500 | 45.0 | 0.113 |
| housing_increment_0to1 | 0.6644 | 0.6808 | +0.0164 | 14.0 | 0.004 |
| prime30_55_childless_renter_mean_rooms | 3.8053 | 3.9681 | +0.1628 | 6.0 | 0.159 |
| prime30_55_childless_owner_share_rooms_ge6 | 0.5961 | 0.7584 | +0.1622 | 25.0 | 0.658 |
| young_childless_renter_liquid_wealth_to_annual_gross_income_2535 | 0.1792 | 0.3226 | +0.1434 | 12.0 | 0.247 |
| prime30_55_childless_owner_minus_renter_mean_rooms | 2.4188 | 2.3132 | -0.1055 | 12.0 | 0.134 |
| own_rate_2534 | 0.3412 | 0.2760 | -0.0652 | 80.0 | 0.340 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.3677 | 0.6017 | +0.2340 | 8.0 | 0.438 |
| old_total_estate_wealth_to_annual_income_median_7684 | 6.5013 | 6.4305 | -0.0708 | 18.6 | 0.093 |
| old_nonhousing_ge_1x_income_share_6575 | 0.6083 | 0.6100 | +0.0017 | 9435.2 | 0.026 |
| old_age_own_rate | 0.7643 | 0.9540 | +0.1897 | 160.0 | 5.760 |
| aggregate_mean_occupied_rooms_18_85 | 5.7800 | 5.6721 | -0.1078 | 6.0 | 0.070 |

Untargeted: estate p90/p50 (76-84) 1.751; owner_neg_liquid_share_2534 0.569.

## Cell 3 s12_lowc — full target fit (loss 639.8888, residual 2.50e-06)

receipt 0.0099% of households (renter 0.0000%, owner 0.0099%; ages 18-33 0.0004%, 66+ 0.0000%); outlays 0.0000% of gross income; entry_censored_share 0.0000

| moment | target | model | gap | weight | loss contrib |
|---|---:|---:|---:|---:|---:|
| tfr | 1.9180 | 2.5687 | +0.6507 | 20.0 | 8.468 |
| childless_rate | 0.1880 | 0.0922 | -0.0958 | 20.0 | 0.184 |
| own_rate | 0.5755 | 0.7499 | +0.1744 | 100.0 | 3.043 |
| own_family_gap | 0.1677 | 0.1972 | +0.0295 | 45.0 | 0.039 |
| housing_increment_0to1 | 0.6644 | 0.4961 | -0.1683 | 14.0 | 0.397 |
| prime30_55_childless_renter_mean_rooms | 3.8053 | 4.3801 | +0.5748 | 6.0 | 1.983 |
| prime30_55_childless_owner_share_rooms_ge6 | 0.5961 | 0.9386 | +0.3425 | 25.0 | 2.932 |
| young_childless_renter_liquid_wealth_to_annual_gross_income_2535 | 0.1792 | 0.6755 | +0.4962 | 12.0 | 2.955 |
| prime30_55_childless_owner_minus_renter_mean_rooms | 2.4188 | 3.3051 | +0.8863 | 12.0 | 9.427 |
| own_rate_2534 | 0.3412 | 0.4640 | +0.1229 | 80.0 | 1.207 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.3677 | 0.4386 | +0.0709 | 8.0 | 0.040 |
| old_total_estate_wealth_to_annual_income_median_7684 | 6.5013 | 7.4504 | +0.9491 | 18.6 | 16.743 |
| old_nonhousing_ge_1x_income_share_6575 | 0.6083 | 0.3608 | -0.2476 | 9435.2 | 578.234 |
| old_age_own_rate | 0.7643 | 0.9634 | +0.1992 | 160.0 | 6.346 |
| aggregate_mean_occupied_rooms_18_85 | 5.7800 | 6.9268 | +1.1469 | 6.0 | 7.892 |

Untargeted: estate p90/p50 (76-84) 2.177; owner_neg_liquid_share_2534 0.642.

## Cell 4 s20_lowc — full target fit (loss 452.5248, residual 1.88e-06)

receipt 0.9142% of households (renter 0.5531%, owner 0.3611%; ages 18-33 0.5620%, 66+ 0.0000%); outlays 0.0539% of gross income; entry_censored_share 0.0000

| moment | target | model | gap | weight | loss contrib |
|---|---:|---:|---:|---:|---:|
| tfr | 1.9180 | 2.1565 | +0.2385 | 20.0 | 1.138 |
| childless_rate | 0.1880 | 0.2768 | +0.0888 | 20.0 | 0.158 |
| own_rate | 0.5755 | 0.6833 | +0.1078 | 100.0 | 1.163 |
| own_family_gap | 0.1677 | 0.4602 | +0.2925 | 45.0 | 3.850 |
| housing_increment_0to1 | 0.6644 | 0.2970 | -0.3675 | 14.0 | 1.890 |
| prime30_55_childless_renter_mean_rooms | 3.8053 | 3.1935 | -0.6118 | 6.0 | 2.246 |
| prime30_55_childless_owner_share_rooms_ge6 | 0.5961 | 0.8386 | +0.2425 | 25.0 | 1.470 |
| young_childless_renter_liquid_wealth_to_annual_gross_income_2535 | 0.1792 | 1.3281 | +1.1489 | 12.0 | 15.840 |
| prime30_55_childless_owner_minus_renter_mean_rooms | 2.4188 | 4.4752 | +2.0564 | 12.0 | 50.745 |
| own_rate_2534 | 0.3412 | 0.4339 | +0.0927 | 80.0 | 0.687 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.3677 | 0.3390 | -0.0287 | 8.0 | 0.007 |
| old_total_estate_wealth_to_annual_income_median_7684 | 6.5013 | 7.3940 | +0.8927 | 18.6 | 14.810 |
| old_nonhousing_ge_1x_income_share_6575 | 0.6083 | 0.4151 | -0.1932 | 9435.2 | 352.199 |
| old_age_own_rate | 0.7643 | 0.9126 | +0.1484 | 160.0 | 3.521 |
| aggregate_mean_occupied_rooms_18_85 | 5.7800 | 6.4631 | +0.6831 | 6.0 | 2.800 |

Untargeted: estate p90/p50 (76-84) 3.065; owner_neg_liquid_share_2534 0.355.

## Cell 5 s12_midc — full target fit (loss 507.8810, residual 1.12e-05)

receipt 0.1089% of households (renter 0.0000%, owner 0.1089%; ages 18-33 0.0016%, 66+ 0.0000%); outlays 0.0058% of gross income; entry_censored_share 0.0000

| moment | target | model | gap | weight | loss contrib |
|---|---:|---:|---:|---:|---:|
| tfr | 1.9180 | 2.4278 | +0.5098 | 20.0 | 5.198 |
| childless_rate | 0.1880 | 0.1274 | -0.0606 | 20.0 | 0.073 |
| own_rate | 0.5755 | 0.7367 | +0.1613 | 100.0 | 2.600 |
| own_family_gap | 0.1677 | 0.2638 | +0.0961 | 45.0 | 0.416 |
| housing_increment_0to1 | 0.6644 | 0.4998 | -0.1646 | 14.0 | 0.379 |
| prime30_55_childless_renter_mean_rooms | 3.8053 | 4.0419 | +0.2366 | 6.0 | 0.336 |
| prime30_55_childless_owner_share_rooms_ge6 | 0.5961 | 0.8915 | +0.2954 | 25.0 | 2.181 |
| young_childless_renter_liquid_wealth_to_annual_gross_income_2535 | 0.1792 | 0.6063 | +0.4270 | 12.0 | 2.188 |
| prime30_55_childless_owner_minus_renter_mean_rooms | 2.4188 | 3.4250 | +1.0063 | 12.0 | 12.151 |
| own_rate_2534 | 0.3412 | 0.4385 | +0.0973 | 80.0 | 0.758 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.3677 | 0.5910 | +0.2233 | 8.0 | 0.399 |
| old_total_estate_wealth_to_annual_income_median_7684 | 6.5013 | 7.4489 | +0.9476 | 18.6 | 16.688 |
| old_nonhousing_ge_1x_income_share_6575 | 0.6083 | 0.3891 | -0.2192 | 9435.2 | 453.448 |
| old_age_own_rate | 0.7643 | 0.9614 | +0.1971 | 160.0 | 6.216 |
| aggregate_mean_occupied_rooms_18_85 | 5.7800 | 6.6789 | +0.8990 | 6.0 | 4.849 |

Untargeted: estate p90/p50 (76-84) 2.131; owner_neg_liquid_share_2534 0.602.

## Cell 2 lowc_nofloor — infeasible (stage forward_age_22, dead mass 0.000320)

- age 22, b -0.814, z 0.186, income 0.440, tenure 0, transfer 0.000, slack -0.090
- age 22, b -0.674, z 0.186, income 0.440, tenure 0, transfer 0.000, slack -0.078
- age 22, b -0.535, z 0.186, income 0.440, tenure 0, transfer 0.000, slack -0.067

## Cell 6 s20_midc — infeasible (stage forward_age_26, dead mass 0.000017)

- age 26, b -0.535, z 0.186, income 0.575, tenure 0, transfer 0.145, slack 0.053

## Cell 7 gate_check_smallG — infeasible (stage forward_age_22, dead mass 0.006057)

- age 22, b -0.535, z 0.186, income 0.440, tenure 0, transfer 0.080, slack -0.846
- age 22, b -0.535, z 0.380, income 0.900, tenure 0, transfer 0.000, slack -0.466
- age 22, b -0.535, z 0.778, income 1.841, tenure 0, transfer 0.000, slack 0.475

## Reproduce

```
PYTHONPATH=code/model code/model/.venv/bin/python3 code/model/tools/run_intergen_transfer_floor_probe.py
```

Round-1 (cash-on-hand means test) summary lines are preserved in this README's history section only as a caution: receipt 79%, outlays 18.9% of income, loss 3342 in floor_only — the leveraged-owner artifact. Full records: results.jsonl / results.json (v2 only).

