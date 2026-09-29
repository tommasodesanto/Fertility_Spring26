# Saved selected-candidate lifecycle comparison

These are 2007 stationary approximations. The two-birth rule is experimental; neither selected case is an adopted reference. Both separately normalize completed fertility to 2.1 and pass demographic renewal. The frozen block0506 reference is used only to validate saved measurement logic and is not plotted. There is no 2023 transition result here.

The two-birth model allows one optional additional birth attempt after a first success within a four-year cell. The original permits at most one birth per cell. The candidates were recalibrated separately, so their difference does not isolate the causal effect of this rule. Earnings, entry distributions, transfers and floors, preference forms, target values and weights are retained.

## Measurement

CPS 2004/2006 women provide pooled cross-sectional children-ever-born stocks at five-year interview ages 20–44. Counts are capped at three; motherhood is the share with at least one child ever born; children among mothers is the capped mean conditional on motherhood. The model integrates four-year pre/post birth masses over matching five-year windows, with uniform within-cell timing, then aggregates mass before dividing. The exact age-25 target is a separate single-year projection, not the 25–29 bin. Completed fertility 2.1 is a normalization and is distinct from the capped 40–44 stock.

NCHS first-birth shares pool 2003–2006 period first births. Data ages 12–21 and 42–49 are folded into the endpoint model cells. The plotted shares sum to one and reproduce each saved mean first-birth age and share of first births at age 30 or later. These data do not track the CPS women.

The housing plot uses the saved age-node `lifecycle_2023.csv` from each stationary 2007 run; its filename is inherited. ACS 2007 homeownership covers household heads in a 42-metro sample. This age-profile overlay is descriptive, not a new calibration target. Rooms are uncapped in the saved model curves; ACS rooms are capped at nine and are therefore omitted. Liquid wealth is plotted in model income units; no matched empirical age curve is asserted.

## Quantitative comparison

| Item | Original selected | Two-birth selected |
|---|---:|---:|
| Scored loss | 7.826 | 7.842 |
| Children ever born at exact age 25 | 0.530 | 0.606 |
| Mother share at exact age 25 | 0.448 | 0.438 |
| Children among mothers at exact age 25 | 1.185 | 1.385 |

CPS exact-age-25 comparison: children 0.810, mother share 0.457, children among mothers 1.770. The symmetric decomposition in `age25_symmetric_decomposition.csv` is an arithmetic identity, not causal attribution. `target_costs.csv` holds all 14 target and validation rows, including weights and loss contributions; `demographic_comparison.csv` holds normalization and renewal values. `plotted_data.csv`, `first_birth_age_cells.csv`, and `other_lifecycle_data.csv` preserve all plotted values.

`weak_direction_round_centers.csv` reports the smallest right singular vector of each saved final-round 10-by-10 Jacobian in scaled and physical parameter units. Its derivatives are at round centers, not at the selected candidates; numerical rank at relative cutoff 1e-6 is scale and step dependent and does not establish statistical identification.

Source hashes, selected SUCCESS artifact hashes, reference export hashes, and numeric replay checks are in `verification.json`. No checkpoint was opened, no model module imported, and no solve run. The existing 17 standard plots per candidate were not modified.
