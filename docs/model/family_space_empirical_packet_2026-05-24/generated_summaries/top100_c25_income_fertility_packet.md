# ACS/MMS Income-Fertility Cross-Section

Generated: 2026-05-23 23:37:37 EDT

## Design

- Source: ACS/IPUMS `extract27.dta` joined to the MMS center/periphery PUMA lookup.
- Sample: women ages 22-45, non-group-quarter households, years 2012-2023, positive household income.
- Fertility objects are ACS cross-sectional measures: recent birth from `fertyr == 2`, own children in household from `nchild`, and childlessness as no own child in household. They are not completed fertility.
- Middle PUMAs are assigned to: `center`.
- Income bins are within-year household-income quintiles. A residual-income quintile table additionally residualizes log household income on age, year, education, and marital status.

## Main Cross-Section Facts

- Matched women records: `1,884,591`.
- Recent birth rate, bottom income quintile: `0.072`; top income quintile: `0.066`.
- Parent-with-child-under-18 rate, bottom income quintile: `0.555`; top income quintile: `0.522`.
- Childless-in-household rate, bottom income quintile: `0.421`; top income quintile: `0.461`.
- Owner rate, bottom income quintile: `0.259`; top income quintile: `0.830`.
- Mean rooms, bottom income quintile: `4.712`; top income quintile: `7.136`.
- Share of recent births in income quintiles 1-5: 0.245, 0.201, 0.192, 0.188, 0.174.
- Share of childless women in income quintiles 1-5: 0.204, 0.219, 0.209, 0.194, 0.174.

## Center-Periphery Income Decomposition

- Center-minus-periphery recent-birth gap: `-0.0062`; within-income component `-0.0067`; income-composition component `0.0005`.
- Center-minus-periphery parent-with-child-under-18 gap: `-0.1006`; within-income component `-0.1035`; income-composition component `0.0028`.

## Family-Space Scarcity And Fertility Levels

- Overall center recent-birth rate: `0.062`; periphery: `0.068`.
- Overall center parent-with-child-under-18 rate: `0.431`; periphery: `0.532`.
- Overall center mean own children: `0.936`; periphery: `1.132`.
- Center recent-birth level slope on 3+ bedroom stock scarcity: `-0.0031 (SE 0.0026, p=0.242)`.
- Periphery recent-birth level slope on 3+ bedroom stock scarcity: `0.0026 (SE 0.0018, p=0.166)`.
- Metro-total recent-birth level slope on 3+ bedroom stock scarcity: `0.0014 (SE 0.0017, p=0.43)`.

- Fertility gaps are periphery minus center by metro, so positive values mean fertility is more peripheral.
- Recent-birth gap slope on 3+ bedroom stock scarcity: `0.0077 (SE 0.0027, p=0.00621)`.
- Mean-own-children gap slope on 3+ bedroom stock scarcity: `0.2389 (SE 0.0577, p=8.93e-05)`.
- Parent-with-child-under-18 gap slope on 3+ bedroom stock scarcity: `0.1203 (SE 0.0222, p=6.72e-07)`.

## Regression Reads

- Recent birth slope on standardized log household income: `-0.0075 (SE 0.0006, p=6.4e-21)`.
- Recent birth quadratic in standardized log household income: `0.0024 (SE 0.0002, p=8.09e-26)`.
- Parent-with-child-under-18 slope on standardized log household income: `-0.0393 (SE 0.0015, p=3.47e-40)`.
- Parent centrality interaction with standardized log household income: `-0.0032 (SE 0.0027, p=0.229)`.

## Outputs

- `acs_fertility_by_income_quintile.csv`
- `acs_fertility_by_income_location.csv`
- `acs_fertility_by_income_age.csv`
- `acs_fertility_by_income_tenure.csv`
- `acs_fertility_by_residual_income_quintile.csv`
- `acs_income_distribution_by_fertility_group.csv`
- `acs_income_distribution_by_parity.csv`
- `acs_income_composition_by_location_parent.csv`
- `acs_fertility_by_location_overall.csv`
- `acs_center_periphery_fertility_income_decomposition.csv`
- `acs_income_fertility_regressions.csv`
- `acs_fertility_gaps_by_metro.csv`
- `acs_fertility_gaps_by_metro_with_bedroom_supply.csv`, when bedroom-supply outputs are available.
- `acs_fertility_levels_by_metro_with_bedroom_supply.csv`, when bedroom-supply outputs are available.
- `acs_fertility_supply_gap_regressions.csv`, when bedroom-supply outputs are available.
- `acs_fertility_supply_level_regressions.csv`, when bedroom-supply outputs are available.
- PNG figures for recent births, parenthood/childlessness, location gradients, income distributions, and fertility levels/gaps by bedroom stock scarcity.

## Interpretation Guardrails

Household income is endogenous to marriage, labor supply, and fertility. The residual-income table is a type-composition diagnostic, not an instrument. ACS `nchild` measures own children in the household, not completed fertility.
