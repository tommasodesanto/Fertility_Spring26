# ACS fertility: follow Kleven's matched pseudo-event study

September 17, 2026. Revised following the author's instruction to stay very close to Kleven. Design review only; no estimates produced.

**Recommendation: reproduce Kleven's matching and PSID validation before adapting the outcome to housing.** The earlier version of this memo emphasized fixed-cohort averages and fertility IV. That was not the intended design. Kleven explicitly constructs synthetic pre-birth observations; our inability to observe actual pre-birth housing in ACS does not rule out his approach. The relevant question is whether its matching assumptions work for housing.

The reference is Henrik Kleven, *The Geography of Child Penalties and Gender Norms: A Pseudo-Event Study Approach*, November 2025. The [paper, original code ZIP, and replication guides](../literature/kleven_pseudo_event/README.md) are now local. The [data inventory](acs_fertility_pseudopanel_inventory.md) remains useful, subject to the additional ethnicity requirement below.

## Preserve the published design first

Kleven matches parents observed with an oldest child aged zero in year \(t\), age \(a\), to childless people in year \(t-n\), age \(a-n\). Those matches supply synthetic event time \(-n\); positive event times use observed parents. The baseline matches gender, education, marital status, race/ethnicity, and state. First-birth ages are 25–45; event times span −5 to +10; the omitted period is −2. Exact matches include all ties with matching weights. Separate gender regressions include age and year effects. Outcome timing distinguishes current-week from retrospective annual measures. Validation compares pseudo and actual event studies within PSID/NLSY, including a conditional-independence test. See paper sections 3–4, `matching.R`, `matching_panel.R`, and `table_validation_assumption_conditional_independence.R` in the ZIP.

These details replace the earlier proposed ages 21–35, fixed birthplace cohorts, and immediate twins/same-sex pilot. Keep the original R code read-only and adapt only a separate driver when implementation is authorized. `MASTER.R` runs the entire paper and should not be launched for this pilot.

## Minimal implementation sequence

1. **Create a code-to-variable map.** Trace `clean_acs.R`, `clean_psid.R`, `matching.R`, `matching_panel.R`, and the event-study functions. Pin the paper version, bins, survey/matching weights, event clocks, sample exclusions, unmatched-cell handling, and variance estimator. Existing extract27 has most required fields but lacks `HISPAN`; race alone cannot reproduce the author's white non-Hispanic / black non-Hispanic / Hispanic / other bins. Check all other required cleaning fields before requesting a narrow IPUMS addition. Exact reproduction uses the author's ACS/CPS coverage; our ACS-only overlap is a labeled adaptation, not an exact numerical replication.
2. **Reproduce the labor-market validation on available PSID.** Hide future birth information from the matching algorithm, retain it only for evaluation, and compare synthetic versus observed pre-birth profiles and the resulting event studies on common support. Follow the author's diagnostics before changing matching variables. PSID's biennial interviews must remain observed interviews; no invented annual observations. Check how the supplied cleaning code handles timing before implementing our mapping.
3. **Change only the outcome to housing.** Start with rooms and ownership measured at interview, then bedrooms where available. Report one household outcome per observation or explicitly account for duplicated households. Use both graphs and the within-PSID conditional-independence comparison. Preserve original matching as the benchmark; label every housing-specific adjustment and validate it separately.
4. **Move the validated design to ACS.** First national estimates, then geographic heterogeneity only if support and matching quality justify it. Report matched/unmatched counts, effective weights, pre-birth discrepancies, and uncertainty. Freeze a substantive validation tolerance before looking at housing results; failure to reject equality alone is not validation.

## What housing changes

Kleven targets a gender difference in labor-market responses. Shared housing may respond similarly for both parents, so subtracting fathers from mothers can erase the response we want. Men are not an untreated housing control. Our housing estimand therefore requires its own counterfactual and validation; success for gender gaps does not establish success for household levels.

Education, marriage, and residence are actual baseline matching variables in the paper. Preserve them for reproduction, while explicitly assessing their stability around birth for housing. Moving or changing tenure before birth can be anticipation rather than a bad pre-trend. Never match on rooms, tenure, or another housing outcome simply to force agreement. Coresident-child order remains imperfect; compare the ACS-style roster proxy against PSID histories.

Twins and same-sex instruments remain a later, separate extension. Synthetic pre-birth matching does not identify future twin or sibling-sex status among childless donors, nor establish exclusion from housing demand. The first deliverable should be the Kleven-style first-birth validation, not an IV estimate or an identified aggregate fertility shock.

## Status

Paper and 50 R source files downloaded in the unchanged code ZIP, plus both replication guides. The approximately 5 GB full archive was not retained; raw-data download and estimation were not attempted. A separate OpenCode connectivity test succeeded through the authenticated Kimi provider; it was a tool-availability test, not a review or empirical validation of this design.
