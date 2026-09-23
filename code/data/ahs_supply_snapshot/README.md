# AHS Family-Unit Supply Snapshot

## 2007 national rooms target for the supply level

`build_ahs_2007_room_target.py` reproduces the author's chosen main *quantity*
moment for the 2007 reference economy. It downloads the Census 2007 national
AHS PUF v2.0 flat CSV if absent, verifies its SHA-256, and computes the
`WGT90GEO`-weighted mean of literal `ROOMS` for occupied housing units
(`STATUS=1`) whose householder is aged 18--85 (`HHAGE`). The public-use rooms
topcode is 21, not 9. The script also reports the nine-capped mean as a
within-survey diagnostic and a Fay-BRR standard error from the 160 replicate
weights. It does not change the model's scored moment or calibration contract.

```bash
python3 code/data/ahs_supply_snapshot/build_ahs_2007_room_target.py
```

The ignored receipt is `output/ahs_2007_room_target.json` in this folder. The
full-precision main estimate is 5.729434240102641 rooms (SE
0.008933865036622756; 37,793 sample units); imposing a nine-room cap on the
same AHS records gives 5.684708786004939. The one-off 2005--06 ACS target
5.607885960068579 is from a different survey/year and must not be subtracted
from the AHS estimate as a pure topcode effect. The next numerical calibration
must compare the AHS target with the model's *uncapped occupied physical rooms*,
not its existing nine-capped observer. The supply-level parameter remains free;
this quantity moment identifies it only conditional on the rest of demand and
the chosen reference-cost normalization.

Sources: [AHS 2007 PUF](https://www.census.gov/programs-surveys/ahs/data/2007/ahs-2007-public-use-file--puf-/2007-ahs-national-puf-microdata.html),
[AHS topcodes](https://www2.census.gov/programs-surveys/ahs/2007/AHS_2007_Topcodes.zip),
[HUD/Census replicate variance guide](https://www.huduser.gov/portal/datasets/ahs/AHSN_Public_Use_Replicate_Weight_abbreviated31OCT12.pdf).
The separately anchored price-and-quantity method is deferred as an extra;
no rent or user-cost mapping has been approved.

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
- `output_ahs_family_unit_menu_national/AHS_2023_FAMILY_UNIT_FIGURE_PACKET.pdf`
- `output_ahs_family_unit_menu_metro/AHS_2023_FAMILY_UNIT_FIGURE_PACKET.pdf`

The main constructed objects are:

- Small units: `0-1` bedrooms.
- Middle units: `2` bedrooms.
- Family-sized units: `3+` bedrooms.
- Missing-middle proxy: 2-3 bedroom units in attached single-family, 2-4 unit
  multifamily, or 5-19 unit multifamily structures.
- Family-sized rental scarcity: `1 - Pr(3+ bedrooms | rental stock)` within
  metro.
- Absolute stock counts: `ahs_absolute_unit_counts.csv`, reported in units and
  millions of units.

The purpose is descriptive. These tables identify whether the family-capable
stock is bundled with detached ownership, rental tenure, or moderate
attached/multifamily structures. They do not by themselves identify the causal
effect of family-sized housing supply on fertility.
