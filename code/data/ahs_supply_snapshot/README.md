# AHS Family-Unit Supply Snapshot

This folder builds a fast American Housing Survey stock/menu snapshot for the
Couillard-style family-sized housing question.

Main script:

- `analyze_ahs_family_unit_menu.R`: downloads the 2023 AHS PUF CSV files when
  missing, reads the household file, and tabulates the housing menu by bedroom
  count, rooms, tenure, structure type, rents, child status, and metro.

Default run:

```bash
Rscript code/data/ahs_supply_snapshot/analyze_ahs_family_unit_menu.R
```

National run:

```bash
AHS_SAMPLE=national Rscript code/data/ahs_supply_snapshot/analyze_ahs_family_unit_menu.R
```

Raw AHS ZIP files are downloaded to `raw/` and are ignored by git. Generated
tables and figures are written to `output_ahs_family_unit_menu_metro/` or
`output_ahs_family_unit_menu_national/`, also ignored by git.

Key output packets:

- `output_ahs_family_unit_menu_national/AHS_2023_FAMILY_UNIT_MENU.md`
- `output_ahs_family_unit_menu_metro/AHS_2023_FAMILY_UNIT_MENU.md`

The main constructed objects are:

- Small units: `0-1` bedrooms.
- Middle units: `2` bedrooms.
- Family-sized units: `3+` bedrooms.
- Missing-middle proxy: 2-3 bedroom units in attached single-family, 2-4 unit
  multifamily, or 5-19 unit multifamily structures.
- Family-sized rental scarcity: `1 - Pr(3+ bedrooms | rental stock)` within
  metro.

The purpose is descriptive. These tables identify whether the family-capable
stock is bundled with detached ownership, rental tenure, or moderate
attached/multifamily structures. They do not by themselves identify the causal
effect of family-sized housing supply on fertility.
