# ACS/MMS Couillard-Style Bedroom Supply Packet

Generated: 2026-05-23 15:28:18 EDT

## Design

- Source: ACS/IPUMS `extract27.dta` joined to the MMS center/periphery PUMA lookup.
- Household-head sample, years 2012-2023, non-group-quarter households.
- Bedroom bins follow Couillard's size contrast: `B_0_1`, `B_2`, and `B_3plus`.
- Middle PUMAs are assigned to: `center`.
- Supply-response diagnostic compares early 2012-2014 to late 2021-2023. It is descriptive, not causal.

## Main Facts

- Matched household-head records: `4,103,889`.
- Metros with identified 3+ bedroom price premium: `42`.
- Weighted childless-minus-parent center-share gap: `0.139`.
- Weighted parent-minus-childless owner-share gap: `0.212`.
- Weighted parent-minus-childless bedroom gap: `0.991` bedrooms.
- Weighted mean central 0-1/3+ bedroom stock scarcity: `0.914`.
- Weighted mean central 3+ bedroom price premium: `-0.020` log points.
- Weighted mean change in central 3+ bedroom stock share, early to late: `-0.029`.

## Regression Reads

- Parent centrality gap slope on 3+ bedroom stock scarcity: `0.1467 (SE 0.0256, p=1.3e-06)`.
- Parent centrality gap slope on 3+ bedroom price premium: `-0.0029 (SE 0.0688, p=0.967)`.
- Central 3+ bedroom stock-share change slope on price premium: `0.0529 (SE 0.0603, p=0.388)`.
- Central 3+ bedroom stock-share change slope on stock scarcity: `0.0160 (SE 0.0189, p=0.405)`.

## Outputs

- `acs_bedroom_supply_cells.csv`
- `acs_bedroom_supply_panel.csv`
- `acs_bedroom_location_summary.csv`
- `acs_bedroom_renter_price_cells.csv`
- `acs_bedroom_premia_by_metro.csv`
- `acs_bedroom_parent_centrality_by_metro.csv`
- `acs_bedroom_early_late_supply_response.csv`
- `acs_bedroom_supply_regressions.csv`
- `acs_bedroom_summary_stats.csv`
- PNG figures for the bedroom menu, parent centrality gap, and descriptive stock-share response.

## Interpretation Guardrails

This packet uses bedrooms, which is closer to Couillard than the room-bin packet. It is still descriptive. The supply-response exercise uses observed ACS stock changes, not an exogenous supply shifter, and should not be read as a causal elasticity.
