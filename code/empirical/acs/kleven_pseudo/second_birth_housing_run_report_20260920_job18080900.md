# Second-birth housing diagnostic: job 18080900

Torch job `18080900` completed on 2026-09-20 in 41 seconds with 4 CPUs and 64 GB. The launcher passed the R dependency smoke and the complete source closure before reading the saved exact-support checkpoint. It used the corrected 18080591-derived checkpoint directory and wrote compact receipts to `output/second_birth_housing_20260920_job18080900/`.

The driver used the retained source identity `(YEAR, SAMPLE, SERIAL, PERNUM)` to attach housing fields from `ROOMS_RAW`, `BEDROOMS_RAW`, and `OWNERSHP_RAW`. It applied the reviewed coding rules: `ROOMS` capped at 9 with code 28 unknown, `BEDROOMS` transformed by (x-1) and capped at 5 with code 23 unknown, and `OWNERSHP` codes 1 and 2 mapped to 1 and 0. The source packet contains 2,190,987 rows from the verified Northeast 2005--2019 overlap. Used post and donor rows were checked for `SEX=2`; source-household clusters are `YEAR:SAMPLE:SERIAL`.

All three declared specifications completed all three housing outcomes. The primary and `FERTYR` event-0-yes specifications each used 5,047 gap-(geq 2) anchors and 29,546 fitted rows per outcome. The joint-negative specification used the 1,734 anchors supported at both event -2 and event -1 and 27,956 fitted rows per outcome. Each specification has a three-row `fit_status.csv` with `FIT_COMPLETE` for rooms, bedrooms, and ownership, and a three-row scalar contrast receipt.

| Specification | Anchor/support definition | Rows per outcome | Rooms +3 minus -1 | Bedrooms +3 minus -1 | Ownership +3 minus -1 |
| --- | --- | ---: | ---: | ---: | ---: |
| Primary | all gap-(geq 2) anchors; 5,047 anchors | 29,546 | -0.0493 (0.0518) | 0.0702 (0.0259) | -0.1246 (0.0136) |
| Joint negative | donors at both -2 and -1; 1,734 anchors | 27,956 | -0.1127 (0.0569) | 0.0379 (0.0285) | -0.1515 (0.0145) |
| FERTYR event-0 yes | gap-(geq 2) anchors with observed FERTYR yes | 29,546 | -0.0493 (0.0518) | 0.0702 (0.0259) | -0.1246 (0.0136) |

Numbers in parentheses are source-household clustered standard errors. Ownership estimates remain proportions in the CSV and fit objects; the local diagnostic figure converts them to percentage points for display. The figure is [second_birth_housing_event_curves.png](output/second_birth_housing_20260920_job18080900/second_birth_housing_event_curves.png).

The support receipts report 5,017 negative-donor rows at event -2 and 4,177 at event -1 under the primary and FERTYR specifications. The joint-negative specification has 4,260 rows at -2 and 3,344 at -1; post rows are unchanged. The primary `ROOMS` coding retains 107 raw code-28 observations as missing for outcome analysis. No raw `BEDROOMS=23` or `OWNERSHP=3/9` values occur in this packet. Missing or unknown outcomes remain missing and are not imputed.

The FERTYR sensitivity does not change the fitted sample because the upstream strict anchor definition already excludes observed `FERTYR=no` at event 0 while retaining unknown values. `NCHILD=2` and exactly two linked children define the strict anchor population, so mothers with observed third children are excluded; this changes the negative-composition population and is a design restriction, not evidence about the excluded group.

These are mothers/person-weighted descriptive housing fits with source-household clustered uncertainty, conditional on the constructed coresident roster proxy and exact donor matches. They do not establish historical first-two-birth links, biological parentage, an exogenous fertility shock, or a causal housing response. The machine-readable receipt is `output/second_birth_housing_20260920_job18080900/housing_run_receipt.json`; fit failures, had there been any, would be recorded explicitly rather than filled with zeros.
