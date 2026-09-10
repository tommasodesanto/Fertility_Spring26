# Initial housing targets: bounded feasibility investigation

## Task and decision

**A defensible initial-period housing target set is feasible today under an explicit new measurement contract.** The local raw ACS contains 2005, 2006 and 2007 annual observations, and all 42 cities in the present calibration are identifiable through the existing harmonized MET2013 variable in every checked year. This reverses the stronger inference that missing active 2007 MMS cache means initial housing data are unavailable.

The viable construction preserves the 42 city codes, uses IPUMS MET2013 membership directly, and consistently observes rooms as min(rooms,9). It does **not** preserve the old admitted-PUMA footprint. Geographic membership, room measurement, corresponding model observation rules, uncertainty, and family-group decisions require a new named contract before estimation. No active target/model/parameter was changed here.

Route: bounded empirical source audit plus one selected-year raw pass;20 minute limit. Deliverable: empirical feasibility, source receipts, candidate moments and uncertainty. Lead owns identification and adoption.

## Exact proposed replacement list

All initial observations pool annual 2005–2007 by adding household-weighted numerators and denominators; no arithmetic average of year-specific ratios. The 2023 column is a comparable single-year diagnostic under the same42city rule, not the existing 2012–2023 target pool. Ownership rows retain the DUE structure restriction UNITSSTR 3:10; rooms rows retain their broader housing sample.

| Candidate empirical row | Initial2005–2007 | Metro-bootstrap SE | Same-definition2023 | Corresponding role |
|---|---:|---:|---:|---|
| Mean occupied rooms, capped at9, heads18–85 | 5.570068596577 | 0.088822794169 | 5.583722836148 | Housing supply scale H0 |
| Ownership, heads30–55, DUE structures | 0.647545995948 | 0.020782177947 | 0.587409170952 | Owner premium chi |
| Recent-parent minus no-child ownership, heads30–55 | 0.160079038798 | 0.005598405174 | 0.211160544631 | Hard overidentifying tenure/fertility-selection restriction; group mapping open |
| 3+minus1–2 resident own-child capped rooms, heads30–55 | 0.348551908260 | 0.057789871295 | 0.354798452363 | Per-child housing requirement jointly with first-child requirement |
| Ownership, heads25–34, DUE structures | 0.429682451785 | 0.021084848725 | 0.356026881044 | Optional overidentifying/validation row; no new free parameter |

Initial housing sample: 1,535,537 household records; HHWT total 159,325,005. Prime-age DUE ownership: 809,378 records. Recent-parent/control groups: 43,707 and 331,485 records. Family-room groups: 106,666 large-family and 306,068 small-family records. All 42 metro clusters represented in each of the five checked years 2005, 2006, 2007, 2012, 2023.

Bootstrap: 1,000 draws, seed 20260910; resample 42 metros with replacement, retain each selected metro's household/group totals, and form ratios and differences anew. The same draws are paired across dates. This follows the existing room-audit cluster approach and supplies covariance matrices, but is **not an official ACS survey-design/replicate-weight standard error**. No bootstrap-derived weight was adopted. Young ownership remains optional/default-off pending an explicit decision.

The reviewed PSID first-birth room response remains 0.7202462623815278 (SE 0.0852600513385958), with its original event-window measurement and source. It remains a pooled structural restriction; do not pretend it is a 2005–2007 event-study estimate, cap it at 9 silently, or rerun its reviewed regression. Using a capped ACS stock moment and an uncapped PSID response is permissible only with separate, explicit model observation operators.

## Feasible windows and geography

| Window and geography | Moments | Status today | Remaining step |
|---|---|---|---|
|2005–2007, same42MET2013city codes|All five candidate rows above, min(rooms,9)|Measured from local raw annualACS; point/metro-bootstrap covariance saved|Lead review and explicit new target contract/model observation mapping|
|2007alone, same42city codes|Same five rows|Measured with bootstrap; candidateCSV includes them|Choose single year versus short pool|
|2012and2023separately, same42city codes|Same five rows|Measured with bootstrap; usable endpoint comparisons|No claim that these two years reconstruct pooled2012–2023|
|2012–2023, same42city codes, without old PUMA admission restriction|Same five rows|Computable locally today; not yet extracted in this bounded pass|One separately authorized remaining-year pass; new source receipt|
|2005–2007 with exactly the old admitted-PUMA footprint|Same conceptual rows|Not established; PUMA definitions change|A defensible spatial concordance/area-intersection design is needed; city labels alone are insufficient|
|Four-city subset matching city-only and admitted-PUMA selection at checked2012/2023 dates|Same rows|Components available; unnecessary narrowing for main proposal|Still no guarantee of exactly unchanged early physical footprint; do not call it equivalent to42cities|
|2005/2007 PSID old-age dispersion|Living reference persons 76–84, p90/p50 of NETWORTHR/INCFAMR|RawPSIDlocal, but no already-saved initial-window point/SE|Fresh narrowly selected PSID pass/bootstrap; not run in this task|

## What the geographic check established

IPUMS MET2013 assigns 2013 MSA definitions to 2005 onward ACS using available PUMAs. Its official documentation allows omission/commission error, suppresses city codes with combined error at least 15%, and warns that PUMA changes can alter the represented population despite unchanged city labels. All 42 labels are observed in the actual checked data. This is a standard, explicitly approximate harmonized city definition, not exact constant land-area membership. [IPUMS MET2013](https://usa.ipums.org/usa-action/variables/MET2013).

The present custom MMS procedure additionally requires the exact(state,PUMA,MET2013)triple to appear in its lookup. Aggregating center,middle and periphery does not undo that selection. Among heads in the same42MET2013cities, removing this restriction adds 48.9035% of head weight in 2012 and 4.5439% in 2023. It changes coverage in 37 cities in 2012 and 16 in 2023. Four codes have no extra households at BOTH endpoint checks:12580(Baltimore),15380(Buffalo),29820(LasVegas),33100(Miami). This is empirical evidence against silently treating the current admitted-PUMA sample as the whole42-city population. The investigation did not diagnose why those historical admissions differ.

An apparent pre2012crosswalk exists at code/data/Spatial_aggregate_withmicrodata/crosswalks/geocorr_puma2000_to_cbsa2013.csv, but its contents are an HTML ApplicationError page, not a usable crosswalk. Older county-to-CBSA files do not recover exact within-PUMA household location. No new crosswalk was downloaded.

## Room observation and family definitions

ACS2005–2007codes 9 rooms as 9 or more. From 2008 the national 9 room cap ends. For the proposed new targets use z(h)=min(h,9), in data and in the model, BEFORE any household weighting or family differencing. This changes measurement, not housing preferences, supply or actual housing services. Do not cap the model's economic choices at 9. [IPUMS ROOMS](https://usa.ipums.org/usa-action/variables/ROOMS).

9.7179%of initial head weight is in the9+bin. The 2023 same-city sample has 8.3606% strictly above 9 rooms. Its uncapped mean is5.789219, versus5.583723after harmonization; the family gap is0.448516uncapped versus0.354798capped. These changes are large enough that the uncapped historical level comparison is misleading.

The new sample preserves empirical family definitions: recent parents have NCHILD>0 and ELDCH<4; controls have NCHILD=0. The family-room row compares NCHILD3+ with 1–2 conditional on YNGCH<18. Neither is automatically identical to model any-dependent versus never-parent or literal dependent counts. Resolving these mappings remains necessary. Initial pooling and MET2013harmonization do not fix that issue. The first-child housing jump and per-child housing requirement retain two housing restrictions; ownership still disciplines tenure premium and meanrooms still disciplines supply scale. Joint/local numerical identification must be checked by the lead after a complete new moment system is specified.

## PSID saving and old wealth availability

The earlier audit already established an initial 2005+2007 ratio of aggregate wealth to gross labor earnings 6.926583791073 from saved annual numerator/denominator totals. The current initial-era point is therefore available without another raw scan. A new person-bootstrap SE is still needed.

Old-age dispersion is a different object: weighted p90 divided by median of NETWORTHR/INCFAMR, living reference persons 76–84, income>1000, observed completed children, positive IW. Its existing 3.448110754 (SE .1325) uses 1984–2019. Saved outputs contain no year-specific tail components permitting a valid early quantile reconstruction: quantiles cannot be recovered by averaging the long-pool numbers. The raw 5.87 GiB PSID file is local, so remeasurement is feasible, but small-window support/precision is unverified. Retaining the long-pool value is an explicit stable structural restriction, not a measured initial level. No second large-data scan was conducted.

## Verification, provenance and files

One bounded raw ACS pass read only 2005, 2006, 2007, 2012, 2023 by fixed-record offsets, at most 250,000 rows per chunk. Header parsing used a file-backed subclass because pandas 1.5.3 otherwise copies the full 9.24 GiB file into memory; no active builder was changed. The first header initialization failed before data reading and was corrected. The completed pass wrote its components; a subsequent numpy-bool JSON serialization failure was repaired from saved components without repeating the raw pass.

The selected 2007 and 2023 national head counts, household weights, owner weights and weighted-room sums exactly reproduce the existing canonical national receipts. All four original current-geography 2023 target values also reproduce within 1e-10 after applying the existing admitted-PUMA rule. Thus both raw units and the retained target masks have independent overlap checks. The raw source was not rehashed during this bounded pass; its existing canonical SHA, current size/mtime, overlapchecks and new componenthash are explicitly distinguished in the source receipt.

Owned files (all under this report's housing/subdirectory):

- inspect_early_housing.py: reproducible bounded selected-year extraction.
- active_metros.txt: exact42metro identifiers from the pinned currentcache.
- early_housing_metro_components.csv: per-year/metro/admission numerators, denominators and counts.
- early_housing_source_receipt.json: raw support and overlap checks, source metadata and hash provenance.
- current_2023_overlap_verification.json: independently reproduces four current2023objects.
- summarize_early_housing.py: saved-component candidate/bootstrap builder; no raw access.
- early_housing_feasibility_summary.csv: literal versus capped means/gaps and sample counts.
- early_housing_target_candidates.csv: exact candidate values and uncertainty for four windows.
- metro_bootstrap_covariance_WINDOW.csv and metro_bootstrap_draws_WINDOW.csv: saved uncertainty results.
- early_housing_candidate_receipt.json: bootstrap and candidate provenance.

Authoritative existing source paths:

- /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta
- /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/moment_standard_errors/cache/acs_analysis_samples.rds
- /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/audit_ownership_targets.R
- /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/moment_standard_errors/build_active_acs_room_target_receipt.R
- /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/data/puma_mms_lookup_2010.csv
- /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/data/puma_mms_lookup_2020.csv
- /Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta
- /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/audit_intergen_bequest_family_size_targets.R

No model/calibration/target source, existing empirical output, protected manuscript, cluster job or git state was changed. Lead review of the new geographic/measurement contract and unresolved family mappings is the next decision.
