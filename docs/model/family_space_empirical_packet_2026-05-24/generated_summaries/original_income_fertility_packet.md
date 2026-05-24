# ACS/MMS Income-Fertility Cross-Section

Generated: 2026-05-23 18:27:15 EDT

## Design

- Source: ACS/IPUMS `extract27.dta` joined to the MMS center/periphery PUMA lookup.
- Sample: women ages 22-45, non-group-quarter households, years 2012-2023, positive household income.
- Fertility objects are ACS cross-sectional measures: recent birth from `fertyr == 2`, own children in household from `nchild`, and childlessness as no own child in household. They are not completed fertility.
- Middle PUMAs are assigned to: `center`.
- Income bins are within-year household-income quintiles. A residual-income quintile table additionally residualizes log household income on age, year, education, and marital status.

## Main Cross-Section Facts

- Matched women records: `1,589,272`.
- Recent birth rate, bottom income quintile: `0.071`; top income quintile: `0.066`.
- Parent-with-child-under-18 rate, bottom income quintile: `0.546`; top income quintile: `0.519`.
- Childless-in-household rate, bottom income quintile: `0.429`; top income quintile: `0.466`.
- Owner rate, bottom income quintile: `0.253`; top income quintile: `0.826`.
- Mean rooms, bottom income quintile: `4.648`; top income quintile: `7.115`.
- Share of recent births in income quintiles 1-5: 0.244, 0.199, 0.192, 0.189, 0.176.
- Share of childless women in income quintiles 1-5: 0.204, 0.221, 0.209, 0.194, 0.171.

## Center-Periphery Income Decomposition

- Center-minus-periphery recent-birth gap: `-0.0053`; within-income component `-0.0057`; income-composition component `0.0004`.
- Center-minus-periphery parent-with-child-under-18 gap: `-0.0922`; within-income component `-0.0929`; income-composition component `0.0007`.

## Family-Space Scarcity And Fertility Levels

- Overall center recent-birth rate: `0.063`; periphery: `0.068`.
- Overall center parent-with-child-under-18 rate: `0.451`; periphery: `0.543`.
- Overall center mean own children: `0.957`; periphery: `1.149`.
- Center recent-birth level slope on 3+ bedroom stock scarcity: `-0.0058 (SE 0.0034, p=0.0918)`.
- Periphery recent-birth level slope on 3+ bedroom stock scarcity: `0.0028 (SE 0.0031, p=0.367)`.
- Metro-total recent-birth level slope on 3+ bedroom stock scarcity: `-0.0007 (SE 0.0028, p=0.8)`.

- Fertility gaps are periphery minus center by metro, so positive values mean fertility is more peripheral.
- Recent-birth gap slope on 3+ bedroom stock scarcity: `0.0105 (SE 0.0036, p=0.00592)`.
- Mean-own-children gap slope on 3+ bedroom stock scarcity: `0.2954 (SE 0.0623, p=2.98e-05)`.
- Parent-with-child-under-18 gap slope on 3+ bedroom stock scarcity: `0.1293 (SE 0.0237, p=3.15e-06)`.

## Regression Reads

- Recent birth slope on standardized log household income: `-0.0069 (SE 0.0006, p=1.19e-13)`.
- Recent birth quadratic in standardized log household income: `0.0025 (SE 0.0002, p=5.67e-18)`.
- Parent-with-child-under-18 slope on standardized log household income: `-0.0391 (SE 0.0016, p=2.18e-26)`.
- Parent centrality interaction with standardized log household income: `-0.0060 (SE 0.0045, p=0.193)`.

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
