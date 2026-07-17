# Intergen Mechanics Packet

Diagnostic mechanics packet for the June 2026 one-market intergenerational model. This is not a production calibration or quantitative policy estimate.

Source: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/intergen_income_disciplined_recalibration_20260716/report/results.json`.
Target set: `candidate_replacement_income_disciplined_v1`.

## Headline Fit

| Moment | Model | Target |
|---|---:|---:|
| `tfr` | 2.035 | 1.918 |
| `childless_rate` | 0.189 | 0.188 |
| `own_rate` | 0.658 | 0.575 |
| `own_rate_2534` | 0.276 | 0.341 |
| `old_age_own_rate` | 0.954 | 0.764 |
| `prime30_55_childless_renter_mean_rooms` | 3.968 | 3.805 |
| `prime30_55_childless_owner_mean_rooms` | 6.281 |  |
| `prime30_55_childless_owner_minus_renter_mean_rooms` | 2.313 | 2.419 |

## Files

- `diagnostics/`: standard plots from `intergen_housing_fertility.diagnostics.write_diagnostics`.
- `first_look_policies_markets_on_path.png`: simpler 2-by-2 inspection panel using tenure-probability expected policies and dropping zero-mass points in each displayed childless-renter age/income slice.
- `first_look_policies_markets.png`: same simple panel on the displayed wealth support without the on-path mass filter; use it for numerical grid inspection, not first-pass economics.
- `first_look_policies_markets_full.png`: fuller version with all selected income states and ages.
- `first_look_policies_markets_total_wealth_on_path.png`: on-path first-look panel against expected total wealth after tenure choice.
- `first_look_policies_markets_total_wealth.png` and `_full.png`: same first-look policy panels with the policy x-axis set to expected liquid wealth plus net liquidated housing value after tenure choice.
- `first_look_wealth_density.png` and `.csv`: aggregate ergodic mass over liquid wealth `b`, plus non-stacked renter/owner and childless/parent density panels.
- `first_look_total_wealth_density.png` and `.csv`: same density diagnostics over total wealth `b + (1 - psi) p H` for owners and `b` for renters.
- `wealth_grid_coverage.csv`: how much of the liquid-wealth grid is economically occupied under several mass thresholds.
- `total_wealth_grid_coverage.csv`: occupied support under the total-wealth definition.
- `first_look_policy_lines.csv` and `first_look_market_summary.csv`: source data for the first-look panels, including chosen tenure, ergodic mass, and house-price/user-cost accounting.
- `ergodic_deterministic_policy_states.csv`: deterministic argmax policy rows for every occupied state in the stationary distribution.
- `ergodic_deterministic_policy_top_slice.png`: one literal deterministic policy function for a common occupied childless-renter state.
- `ergodic_deterministic_policy_slices.png`: readable deterministic policy slices for the most common occupied childless-renter age/income cells.
- `ergodic_deterministic_policy_bins.png`: mass-weighted deterministic policies over all occupied states by liquid wealth.
- `target_fit.csv` and `target_fit.md`: full target/model/gap table with loss contributions.
- `room_bin_fit_prime30_55_childless.csv` and `room_bin_shares_prime30_55_childless.png`: model versus ACS room-bin shares.
- `owner_rung_shares_all_owners.csv` and `.png`: realized owner mass across the full owner room ladder.
- `owner_rung_shares_prime30_55_childless.csv` and `.png`: owner rung concentration among prime-age childless owners.
- `age_profiles.csv` and `.png`: lifecycle ownership, fertility, housing, and liquid wealth profiles.
- `tenure_by_age.csv` and `.png`: simple renter/owner population shares by age.
- `owner_entry_thresholds.csv` and `.png`: childless-renter owner-entry probability thresholds by age and income state.
- `owner_entry_policy_childless_renter_age30_42.csv`: owner-entry probability lines used for age-30 and age-42 inspection.
- `classic_draft/`: the two stable quantitative figures used by the current paper draft, with CSV sidecars. These are written automatically for an M4 standard-bequest packet.
- `solution_cache.pkl`: local trusted Python pickle with the solved `sol`, `P`, and `p_eq` objects. Use `--refresh-plots-from-cache` to rebuild figures and CSVs from it without re-solving.
- `contact_sheet.png`: quick visual index when standard diagnostics are enabled.
- `policy_cases/`: diagnostic policy proof-of-concept cases from the same theta.

## Policy Caveat

Policy cases are fixed-theta diagnostic mechanics only. The estate-tax case is a terminal bequest wedge without revenue rebates or inheritance kernels.
