# ACS/MMS Couillard-Style Bedroom Supply Packet

Generated: 2026-05-23 23:33:16 EDT

## Design

- Source: ACS/IPUMS `extract27.dta` joined to the MMS center/periphery PUMA lookup.
- Household-head sample, years 2012-2023, non-group-quarter households.
- Bedroom bins follow Couillard's size contrast: `B_0_1`, `B_2`, and `B_3plus`.
- Middle PUMAs are assigned to: `center`.
- Supply-response diagnostic compares early 2012-2014 to late 2021-2023. It is descriptive, not causal.

## Main Facts

- Matched household-head records: `4,989,066`.
- Metros with identified 3+ bedroom price premium: `80`.
- Weighted childless-minus-parent center-share gap: `0.108`.
- Weighted parent-minus-childless owner-share gap: `0.209`.
- Weighted parent-minus-childless bedroom gap: `0.973` bedrooms.
- Weighted mean central 0-1/3+ bedroom stock scarcity: `1.058`.
- Weighted mean central 3+ bedroom price premium: `-0.030` log points.
- Weighted mean change in central 3+ bedroom stock share, early to late: `-0.053`.

## Regression Reads

- Parent centrality gap slope on 3+ bedroom stock scarcity: `0.0903 (SE 0.0192, p=1.14e-05)`.
- Parent centrality gap slope on 3+ bedroom price premium: `-0.0353 (SE 0.0350, p=0.316)`.
- Central 3+ bedroom stock-share change slope on price premium: `0.1476 (SE 0.0384, p=0.00037)`.
- Central 3+ bedroom stock-share change slope on stock scarcity: `0.0172 (SE 0.0142, p=0.232)`.

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
