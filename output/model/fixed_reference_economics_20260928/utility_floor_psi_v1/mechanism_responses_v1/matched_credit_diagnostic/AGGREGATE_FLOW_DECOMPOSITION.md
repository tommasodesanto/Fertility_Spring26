# Native saved-flow decomposition

policy first at baseline PRE; composition evaluated at credit policy; reverse corner absent. These are accounting differences, not an explanation of household mechanisms. No arrays or model calls are used.

| Flow per 1,000 households per four-year period | Policy at baseline PRE | Composition at credit policy | Total |
|---|---:|---:|---:|
| first_births | +2.775507 | -5.736835 | -2.961329 |
| second_births | +0.445821 | -2.673657 | -2.227836 |
| third_bin_entries | +0.134002 | -1.002174 | -0.868172 |
| births | +3.355330 | -9.412666 | -6.057337 |

saved cohort TFR equals topcode-adjusted birth children divided by actual entry flow; raw event flows are distinct. uniform_birth_time childless_rate_40_44 observer, not all-age childless mass.

First/second/third flows are cumulative children-count differences from PRE to POST fertility, checked to sum to native total events. Risk pools are snapshotted before any birth; fecundability multiplies attempts and no household chains two births within one period. Source identities and the three exact retained corners are in `aggregate_flow_decomposition.json`.

The credit PRE distribution was not separately saved. Granular composition requires exact native deterministic forward reconstruction from the saved post-fertility distribution and policies, with baseline identity and all retained flow checks. No such reconstruction has been executed.
