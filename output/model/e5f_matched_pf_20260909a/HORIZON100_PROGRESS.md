# 100-date perfect-foresight solution: provisional readout

This is a price-solver diagnostic at inherited parameters, not a new calibration. The tables below use the collected trial with the smallest market residual. Final reproduction and horizon stability are not certified by this document.

| Trial | Maximum market gap | Mapping checks | Seconds |
|---|---:|---|---:|
| 1 | 18.199032% | Pass | 3310.3 |
| 2 | 6.984612% | Pass | 3255.7 |
| 3 | 1.914974% | Pass | 3233.6 |

The market tolerance is **0.02%**. Parameters, empirical targets and weights are identical across these trials.

## Full target fit — trial 3

Objective at these provisional prices: **94.30055915**. Shares remain in fraction units. This objective is not a calibrated-equilibrium loss.

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Completed fertility | 1.918 | 1.7771433 | -0.14085666 | 1425.739 | 28.287515 |
| Childless share | 0.188 | 0.23358543 | 0.045585425 | 17180.744 | 35.702118 |
| Mean age at first birth | 26.044627 | 26.23329 | 0.18866241 | 44.444444 | 1.5819335 |
| First births at age 30+ (share) | 0.2603274 | 0.24409669 | -0.016230713 | 10000 | 2.6343605 |
| First-birth housing response (rooms) | 0.72024626 | 0.4119682 | -0.30827806 | 137.56527 | 13.073566 |
| Rooms gap: 3+ versus 1–2 children, ages 30–55 | 0.36769956 | 0.34806426 | -0.019635297 | 2958.515 | 1.1406404 |
| Parent ownership gap (share units) | 0.16766167 | 0.14660198 | -0.021059686 | 14229.591 | 6.3109713 |
| Ownership share | 0.575472 | 0.58771694 | 0.012244938 | 1207.8461 | 0.18110263 |
| Mean occupied rooms, ages 18–85 | 5.7799705 | 6.3250511 | 0.54508066 | 11.973159 | 3.5573803 |
| Wealth / annual gross labor earnings | 6.8731 | 7.0587524 | 0.18565243 | 6.2876694 | 0.21671599 |
| Annual bequests / wealth | 0.0088 | 0.0083540801 | -0.00044591994 | 5165289.3 | 1.0270898 |
| Wealth / income dispersion, ages 76–84 (p90/p50) | 3.4481108 | 3.3465804 | -0.10153036 | 56.959772 | 0.58716492 |

Four ACS targets pool 2012–2023 while the current model observer uses 2023; two parent/child group definitions remain unresolved. The approved childbirth event-study target is unchanged. See overnight_target_mapping_review.md for authoritative sources and dates.

## Complete parameter and restriction table

The eleven free coordinates below are inherited inputs to this price solve. No parameter has been re-estimated in this run. Near-bound flags reproduce the existing parameter receipt.

| Parameter | Value | Lower | Upper | Free coordinate | Near bound | Restriction/status |
|---|---:|---:|---:|---|---|---|
| beta_annual | 0.99527658 | 0.94 | 0.9995 | True | False | fixed or bounded pilot transition parameter |
| kappa_fert | 2.168173 | 0.02 | 50 | True | False | fixed or bounded pilot transition parameter |
| kappa_fert_continuation | 1.7707706 | 0.02 | 50 | True | False | fixed or bounded pilot transition parameter |
| chi | 1.053727 | 0.1 | 5 | True | False | fixed or bounded pilot transition parameter |
| H0 | 13.716049 | 0.2 | 80 | True | False | fixed or bounded pilot transition parameter |
| theta0 | 0.57034989 | 0 | 8 | True | False | fixed or bounded pilot transition parameter |
| theta1 | 0.10372395 | 0.02 | 16 | True | True | fixed or bounded pilot transition parameter |
| hbar_child_rooms | 0.24708491 | 0.1 | 1.8 | True | False | fixed or bounded pilot transition parameter |
| first_birth_fixed_cost | 4.6197314 | 0 | 8 | True | False | fixed or bounded pilot transition parameter |
| hbar_first_child_jump | 0.4697651 | 0 | 0.5 | True | False | fixed or bounded pilot transition parameter |
| psi_child_change_2023 | -0.3213878 | -1.5 | 0.2 | True | False | fixed or bounded pilot transition parameter |
| psi_child_2007 | 0.29005153 | — | — | False | False | externally normalized to old completed fertility |
| psi_child_2023 | -0.031336267 | — | — | False | False | derived from old intercept and transition coordinate |
| tenure_choice_kappa | 0.005 | — | — | False | False | externally fixed profile not estimated |
| housing_supply_elasticity | 0.63 | — | — | False | False | externally fixed profile not estimated |

## Terminal-distance checks for this trial

These endpoint distances do not replace a comparison of historical prices and moments across longer horizons.

| Distance | Value | Tolerance | Pass |
|---|---:|---:|---|
| asset price relative gap | 0.0017831804 | 0.01 | True |
| equal transfer relative gap | 0 | 0.01 | True |
| household heads relative gap | 0.0032804222 | 0.01 | True |
| normalized head age sex l1 | 0.00075906956 | 0.02 | True |
| normalized household distribution l1 | 0.017570167 | 0.02 | True |
| normalized person age sex l1 | 0.00093820274 | 0.02 | True |
| psi absolute gap | 0 | 0.001 | True |
| renter price relative gap | 0.012543935 | 0.01 | False |
| resident persons relative gap | 0.0032250431 | 0.01 | True |

No new matched policy path is available. The main intended property-tax comparison keeps equal household rebates in both tax regimes. Market convergence, replay, horizon stability, empirical alignment and matched re-estimation remain separate requirements.

Source receipts: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/meeting_receipts/historical_root_h100_01/sequential/evaluation_003`. Regenerate with `python3 review_horizon100_root_progress.py` after collecting complete trials. No model solve is performed.
