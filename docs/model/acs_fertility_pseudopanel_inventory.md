# ACS / PSID feasibility inventory

Verified September 17, 2026. Metadata and existing outputs only; no new sample counts or estimates. Paths below are relative to the project root.

## Samples and files

| Object | Verified evidence | Implication / unresolved item |
|---|---|---|
| ACS person extract | `code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta`, approximately 9.2 GB; `haven::read_dta(..., n_max=0)` verified its schema | Child and mother fields coexist with housing fields; no observation-level completeness checked |
| Extract documentation | `raw_data/extractor27.do:7–76` within the same directory; also `usa_00026.xml` and `usa_00028.xml` | XML time coverage lists 2005–2023. This does not verify extract27's observed year counts or sample mix |
| Latest compressed extract | Same directory, `usa_00028.dat.gz` with matching XML | Existing summary reports 59,046,776 rows processed and 10,202,560 retained for ages 22–45, 2012–2023, GQ 1/2. These are person-record counts, not eligible mother/twin counts |
| Origin-geography subset | `extract28_origin_geo_2012_2023_age22_45_households.csv` and `_summary.txt` | Only year/sample/serial/pernum/migration fields; adult restriction and missing child variables make this unsuitable for roster reconstruction |
| Earlier pseudo-panel | Targeted code/data, code/empirical and docs filename searches; targeted text in memory and SESSION_DIARY did not recover a specific cohort implementation | Existing MSA-year aggregation is not proof of a fixed-cohort pseudo-panel. Conversation archives were not exhaustively searched |
| PSID underlying data | `/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta` exists | Needed for roster validation; endpoint extracts alone cannot reproduce longitudinal trajectories |
| Preferred corrected PSID audit | `code/data/psid_followup_mar2026/iv_housing_reaudit_20260809.do`; its `output/iv_housing_reaudit_20260809/` has README, CSV tables, and both clean samples | Mother-level baseline-adjusted five-year audit; distinct from the panel and post-only followups |
| Existing PSID rooms samples | `iv_housing_reaudit_samples.csv`: twins 1,974 mothers / 34 instrument-positive; same-sex 1,836 / 917 | Unweighted sample counts, outcome-specific. These are historical corrected samples, not new validation results |

## Variables in extract27's actual schema

| Purpose | Present fields | Restriction / next check |
|---|---|---|
| IDs and links | `year sample serial pernum momloc poploc sploc relate related` | Link within year/sample/household; retain children before maternal filters. Probable links do not establish biological birth order; inspect link-rule/quality variables in an expanded extract |
| Timing and sex | `age sex eldch yngch` | No `birthqtr` in this schema. Integer-age ties are only a twin proxy; approximate event ages are not exact dates |
| Fertility | `nchild nchlt5 fertyr` | Coresident children and recent birth are not children ever born. No complete birth history or verified twin indicator in this schema |
| Housing | `rooms bedrooms ownershp ownershpd rent valueh mortgage` | Verify valid codes and top-coding by year; rooms, bedrooms, rent, and housing value are distinct outcomes |
| Cohorts / demographics | `bpl bpld race raced age citizen yrimmig` | Freeze cohort definition; national first. Birthplace grouping does not prevent immigration-related composition changes |
| Benchmark labor outcome | `empstat empstatd uhrswork wkswork1` | Check universes; replicate the benchmark's sample and treatment rather than equating different local effects |
| Weights / inference | `perwt hhwt cluster strata` | Maternal estimates use person weights; household estimates require a stated household estimand. No replicate-weight fields in schema; verify survey-design inference before production |
| Geography | `statefip puma met2013 migrate1 migpuma1 migmet131` | Harmonize PUMA vintages for spatial extensions. Migration geography is not current residence; household IDs are not longitudinal identifiers |

[IPUMS MOMLOC](https://usa.ipums.org/usa-action/variables/MOMLOC) documents probable coresident mother links. [IPUMS birth-order guidance](https://forum.ipums.org/t/birth-order-variable-s/6396) explains missing nonresident siblings and lack of complete histories. [BIRTHQTR](https://usa.ipums.org/usa-action/variables/BIRTHQTR) is offered for ACS from 2005: request it to strengthen timing, while retaining the proxy label. Quarter agreement alone does not prove twinning, and integer-age discordance near birthdays can cause false negatives.

Use annual ACS samples without overlapping multiyear extracts. Exclude 2020 from the initial standard-year comparison; [Census guidance](https://www.census.gov/programs-surveys/acs/data/experimental-data/faq.html) does not recommend comparing its experimental one-year estimates with standard estimates. Confirm actual year/sample coverage before setting the pilot years. No fresh extraction was requested here.

Verification: actual extract27 schema read with zero rows; XML coverage and variable names parsed; PSID source existence checked; corrected sample CSV read; source lines checked. R returned locale warnings but completed the schema read successfully. No household microdata were uploaded to an external model service.
