# Moreno-Maldonado / Santamaria Center-Periphery Build

This folder implements a paper-style within-metro downtown/suburb geography for the fertility project.

Baseline design:
- Use the CBD/tract system from the urban-revival data bundle.
- Keep only metro-like regions with `pop2000 >= 1,000,000`.
- Define the core as the smallest set of 2000 census tracts closest to the CBD that contains `MMS_CORE_POP_SHARE` of region population. The default is `0.10`.
- Define the outer suburb ring as the furthest tracts containing 50% of region population.
- Map the tract geography into PUMAs and classify a PUMA as:
  - `center` if at least `MMS_CENTER_PUMA_CORE_SHARE` of its assigned tract population is in core tracts
  - `periphery` if less than `MMS_PERIPHERY_PUMA_CORE_SHARE` of its assigned tract population is in core tracts
  - `middle` otherwise
- Downstream summaries can optionally absorb `middle` into `center` with `MMS_INCLUDE_MIDDLE_IN_CENTER=1`.
- Keep only metros where downtown-classified PUMAs capture at least 50% of downtown tract population.

Important implementation note:
- The original paper maps downtown tracts into PUMAs using more detailed geographic overlaps.
- Here, tract-to-PUMA assignment is done by point-in-polygon using 2000 tract centroids and the local PUMA shapefile already present on this machine.
- For the current project, that is a fast and transparent approximation; the scripts save diagnostics so the approximation can be audited.

Files:
- `../Spatial_aggregate_withmicrodata/build_extract28_origin_geo_supplement.R`: builds a filtered person-level supplement from `usa_00028` with `MIGPLAC1` and `MIGMETRO1`
- `build_mms_geography.R`: builds the downtown/suburb lookup and diagnostics
- `compute_mms_descriptives.R`: applies the lookup to `extract27.dta` and recomputes core validation facts
- `build_migpuma_origin_bridge.R`: uses the official IPUMS PUMA/MIGPUMA relationship files plus the PUMA lookup to build origin MMS shares by `(MIGPLAC1, MIGPUMA1, CBSA)`
- `analyze_mms_fertility_moves.R`: rebuilds the fertility/location and new-parent move facts under the MMS geography using `extract27` plus the filtered `extract28` origin supplement
- `analyze_family_size_supply_menu.R`: builds a fast ACS/MMS family-size housing menu and sorting packet by center/periphery, room-size bin, tenure, parent status, and income proxy
- `analyze_income_fertility_cross_section.R`: builds an ACS/MMS cross-section of recent births, parenthood, childlessness, parity, income ranks, and center/periphery decomposition by income

Outputs:
- `data/migpuma_mms_origin_bridge.csv`
- `data/migpuma_mms_origin_bridge_diagnostics.csv`
- `data/region_to_cbsa.csv`
- `data/tract2000_mms_flags.csv.gz`
- `data/puma_mms_lookup_2010.csv`
- `data/puma_mms_lookup_2020.csv`
- `data/city_coverage.csv`
- `output/mms_location_summary.csv`
- `output/mms_tenure_parent_summary.csv`
- `output/mms_age_profiles.csv`
- `output/mms_summary.md`
- `output/mms_location_profiles.png`
- `output/mms_location_by_parent_summary.csv`
- `output/mms_move_type_summary.csv`
- `output/mms_four_way_move_shares.csv`
- `output/mms_origin_transition_summary.csv`
- `output/mms_origin_regressions.csv`
- `output/mms_fertility_moves_summary.md`
- `output_family_size_supply/ACS_MMS_FAMILY_SIZE_SUPPLY_PACKET.md`
- `output_family_size_supply/acs_family_size_premia_by_metro.csv`
- `output_family_size_supply/acs_parent_centrality_by_metro.csv`
- `output_family_size_supply/acs_family_size_regressions.csv`
- `output_income_fertility_cross_section/ACS_MMS_INCOME_FERTILITY_CROSS_SECTION.md`
- `output_income_fertility_cross_section/acs_fertility_by_income_quintile.csv`
- `output_income_fertility_cross_section/acs_center_periphery_fertility_income_decomposition.csv`
