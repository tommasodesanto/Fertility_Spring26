# National effective residential property-tax rate: bounded source review

**Date:** 2026-09-23. **Status:** recommendation only; no parameter adopted.

## Finding

The preferred national analogue to county taxes-paid divided by housing value is the ratio of national aggregate self-reported annual residential property taxes to aggregate owner-reported home value for the same owner-occupied housing universe. Using the Census Bureau's 2011 ACS five-year summary file (the 2007–2011 observation vintage), the national estimates are $210,300,424,500 in aggregate real-estate taxes (B25090_001E) and $19,842,731,247,500 in aggregate home value (B25082_001E). Their ratio is **0.01059836 = 1.060%**.

This is a ratio of sums, equivalent to a value-weighted mean of household tax/value rates for the owner-occupied housing stock. It is distinct from a simple mean of county rates or household rates and from a population-weighted county average. It is a 2007–2011 pooled ACS estimate, not a 2007-only estimate. ACS self-reported taxes and values are the measured objects; this should not be interpreted as a statutory marginal rate or as commercial property tax.

I retrieved the official 2007–2011 ACS summary-file national sequence 0103. The U.S. geography record 0000001 maps to 01000US / United States. The technical documentation places B25082 in sequence-file positions 56–58 and B25090 in positions 121–123; the first estimate in each table is the national aggregate. Mortgage-status subtotals corroborate the totals: B25082 has $14,351,943,672,500 with a mortgage plus $5,490,787,575,000 without; B25090 has $156,605,955,200 plus $53,694,469,300. The two subtotals sum to their respective aggregate totals.

The 2007 ACS one-year release has aggregate owner-occupied home value (B25082), but its aggregate-tax table B25090 was not in that release. Thus the matched 2007–2011 ACS5 ratio is the closest recovered national counterpart, not an exact single-year 2007 measurement.

The exact-2007 AHS table 3-13 search-indexed row “Annual Taxes Paid per $1,000 Value” reports a national median of $9 per $1,000 (0.9%). That is a median household rate, not an aggregate ratio; I could not retrieve its PDF directly (Census host returned HTTP 403), so it is secondary here and should not override the verified ACS aggregate calculation.

## Recommendation

For a national effective residential property tax input, recommend **1.060% as a proposed matched-period (2007–2011) national ratio-of-sums candidate**, with the five-year vintage and owner-occupied universe disclosed. It is methodologically much closer to the DUE/Brookings taxes-paid-to-home-value approach than the Bay Area 0.71% scalar or the AHS 0.9% household median. It is not a 2007-only estimate and remains a recommendation, not an adopted parameter. If the specification requires a strictly 2007 annual rate, this result does not satisfy that requirement; no exact-year national aggregate-tax ratio was verified in this review.

## Source and vintage record

- U.S. Census Bureau, *American Housing Survey for the United States: 2007*, table 3-13, PDF: https://www.census.gov/content/dam/Census/library/publications/2008/demo/h150-07.pdf ; landing page and tables: https://www.census.gov/programs-surveys/ahs/data/2007/ahs-2007-summary-tables/h150-07.html . The 0.9% is the table's national median annual tax dollars per $1,000 value.
- Census, 2007 ACS subject definitions, confirming the “real estate taxes per $1,000 value” measure is a median: https://www2.census.gov/programs-surveys/acs/tech_docs/subject_definitions/2007_ACSSubjectDefinitions.pdf .
- Official Census summary-file source directory for 2011 ACS5, U.S. sequence files: https://www2.census.gov/programs-surveys/acs/summary_file/2011/data/5_year_seq_by_state/UnitedStates/All_Geographies_Not_Tracts_Block_Groups/ ; technical documentation: https://www2.census.gov/programs-surveys/acs/summary_file/2011/documentation/5_year/ACS_2007-2011_SF_Tech_Doc.pdf .
- Census API variable groups: 2007 ACS B25082 aggregate value: https://api.census.gov/data/2007/acs/acs1/groups/B25082.html ; 2011 ACS5 B25082 aggregate value: https://api.census.gov/data/2011/acs/acs5/groups/B25082.html ; 2011 ACS5 B25090 aggregate real estate taxes: https://api.census.gov/data/2011/acs/acs5/groups/B25090.html .
- Harris and Moore (2013), *Residential Property Taxes in the United States*, Urban-Brookings Tax Policy Center, describes its county comparisons as ACS self-reported residential taxes and values and states that underlying 2007–2011 county-level statistics are county aggregates: https://www.brookings.edu/wp-content/uploads/2016/06/18-residential-property-taxes-harris.pdf . Its reported 1.15% is explicitly the mean of county-level tax/value burdens, not a national housing-value-weighted rate.
- Brookings map note: data are 2007–2011: https://www.brookings.edu/articles/3-things-we-can-learn-about-property-taxes-from-a-map/ .

- Raw files retrieved from the U.S. Census 2007–2011 ACS Summary File, United States / All Geographies Not Tracts Block Groups / sequence 0103: `tmp/national_property_tax/20115us0103000.zip` (SHA-256 `3590fd82ac6f9774b7c327ecc59b163564295807d0954352a31134c47fb06af7`); extracted estimate file `tmp/national_property_tax/e20115us0103000.txt` (SHA-256 `02c98fd8a6a1e957497383dcc05c47f96805046a244c11361e79d44a28b1759c`); geography file `tmp/national_property_tax/g20115us.csv` (SHA-256 `14f4ccbd8ab854a4bce3c8bde3d220c20e8fc2a2e50cf828a385b5c9cb438c47`). Official sequence-file documentation PDF and extracted text are in `tmp/national_property_tax/acs_tech.pdf` and `.txt`.
