# National housing-cost inputs

Author requested national versions of the DUE depreciation and property-tax methods on September 23, 2026. These are measured candidates for the 2007 benchmark, not adopted values or edits to a frozen model. No household/equilibrium solve or fit comparison was run.

## Depreciation

`depreciation/receipt.json` records sources, hashes, definitions and arithmetic checks; `depreciation/estimates.csv` retains full precision. Candidate: **1.416% annually in 2007**.

Following the land-adjustment logic in DUE, compute structures depreciation / structures stock, then multiply by structures / total housing value. BEA tables 5.4 and 5.1, owner-occupied row 11, give $316.854 billion annual depreciation and $13.073488 trillion year-end structures stock. Federal Reserve owner-occupied real estate plus mobile homes gives $22.374424 trillion housing value, excluding separately held vacant land. Thus structures depreciate at 2.424%, the residual land share is 41.569%, and the product is 1.416%.

This uses national 2007 values, not DUE's 2016 BEA rate and later Bay Area land share. It also uses owner-occupied rows consistently instead of the broader household structures rate. Applying that broader household rate gives 1.412%, retained as a definition sensitivity. The 2006 and 2008 equivalents are 1.322% and 1.543%; these are adjacent-year checks, not confidence bounds. The land share is an accounting residual from replacement-cost structures and market-value housing, not a direct appraisal. Using this owner-occupied rate for rented homes remains a common-rate model approximation.

Reproduce offline from the retained official downloads:

```sh
/Users/tommasodesanto/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 code/data/Spatial_aggregate_withmicrodata/build_national_housing_cost_inputs.py
```

Sources: [BEA residential fixed assets tables](https://apps.bea.gov/national/FixedAssets/Release/XLS/Section5All_xls.xlsx), [Fed owner-occupied real estate](https://fred.stlouisfed.org/series/BOGZ1FL155035013A), [Fed mobile homes](https://fred.stlouisfed.org/series/BOGZ1FL155012013A), [Fed mapping of structures stock to BEA row 11](https://www.federalreserve.gov/apps/fof/SeriesAnalyzer.aspx?s=FL155012665&bc=:FL155012665&suf=A&t=). Downloads retained under depreciation/sources. Current downloads can differ from older cached web tables after source revisions.

## Property tax

Candidate: **1.060% annually**, from the **2007–2011 ACS five-year release**. This is $210,300,424,500 aggregate annual taxes divided by $19,842,731,247,500 aggregate home value, both for U.S. owner-occupied housing. It is a ratio of sums (a housing-value-weighted rate), not DUE's population-weighted average of county rates, and not a 2007-only observation. The pooled dates match the underlying Brookings source period. This is an effective self-reported burden, not a statutory marginal tax rate.

`property_tax/calculation_receipt.json` records the independent lead extraction and source hashes. The official U.S. geography record and sequence positions were checked, and the with/without-mortgage subtotals sum exactly in both tables. Compact national source rows are retained; the full 22 MB source archive remains in `tmp/national_property_tax/` and can be retrieved from the URL in the receipt. Both tables cover owner-occupied units. No ratio standard error is asserted without the numerator–denominator covariance.

## Lead recommendation

Use the calculated **1.416% annual depreciation** and **1.060% annual property tax** as the proposed national inputs. The author has authorized remeasurement; adoption of these numbers remains pending. Source geography is national, but both are owner-housing proxies if the model applies the same rates to landlords. Depreciation follows DUE's land adjustment with different, explicitly identified national source series; taxes use the national ratio of sums rather than county population weighting. Full-precision proposed values are in `national_inputs.csv`. No model, frozen objective, paper, or slide edits were made.

The independent depreciation review passes the arithmetic. Lead closes its timing question using BEA sheet row 2 ('yearend estimates') and retained FRED metadata ('Annual, End of Period'); see `lead_review.json`.

