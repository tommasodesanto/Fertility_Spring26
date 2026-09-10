# Overnight target-to-model measurement review

Read-only review, September9 evening/September10 UTC. This note specifies feasible diagnostics; it does not change targets, weights, parameters, model code, the child-departure law, or the estimation objective. The current source of truth is CALIBRATION_STATUS.md, including its latest28-date equilibrium/horizon qualification; the older memory headers describe earlier experiments.

## Direct answer: what year are we targeting?

The model's measurement endpoint is **2023**, reached along a perfect-foresight historical path beginning in2007. We are not estimating twelve independent2023 empirical statistics. The twelve retained empirical rows mix approved completed-fertility cohorts, an event-time housing response, four pooled ACS housing observations, two approved pooled PSID wealth observations and one borrowed bequest normalization.

The observer computes its full target vector at global period4, corresponding to2023. Completed fertility and childlessness refer to the model cohort centered at age42; first-birth timing reconstructs that cohort's past hazards; the housing response follows a fixed2019–2023 matched branch. Other rows use the2023 cross-section. Post2023 values affect these choices through backward induction, but post2023 cross-sections do not enter the current twelve-row loss. A longer horizon therefore can change the2023 fit without changing its measurement date or re-estimating parameters.

## Complete twelve-row mapping

| Active row | Empirical source, population and vintage | Current model operator/date | Status for this review |
|---|---|---|---|
| `tfr` | June2024 CPS fertility supplement; US women40–44; weighted children ever born, with retained public-file/top-bin convention | Age42 cohort's completed-parity stock in2023 | Approved cohort anchor; retained |
| `childless_rate` | Same CPS sample; share reporting zero children ever born | Parity-zero share of that2023 age42 cohort | Approved cohort anchor; retained |
| `mean_age_first_birth` | NCHS natality files1987–2023; first births assigned to maternal cohorts1979–1984; exact ages collapsed to agreed four-year bins/midpoints | Fixed synthetic cohort: old-normalized stationary prehistory plus2007–2023 dated first-birth hazards | Approved cohorts and binning; retained |
| `share_first_births_age30plus` | Same NCHS maternal cohorts; first births at30+ divided by all observed first births | Same cohort ledger, conditional on becoming a mother | Approved cohorts; retained |
| `housing_increment_0to1` | Verified PSID Sun–Abraham occupied-rooms contrast from event−1 to+3; selected reference/spouse woman per household-year, confirmed childless controls | Equal birth/control branches from the2019 risk set;2019 policies advance,2023 policies measure continuation/housing; no Census age bridge in the matched branch | Retain reviewed estimator and current four-year mapping; no regression reassessment |
| `prime30_55_parent_3plus_minus_1to2_mean_rooms` | ACS2012–2023 pooled,42MMS metros; HHWT heads30–55 with youngest resident own child under18; mean rooms for any-age NCHILD>=3 minus NCHILD1–2 | Mean occupied rooms for dependent-count m=3 minus m=1/2, current2023 households in the model age band | Family-count mapping and pooled-date seam require explicit treatment |
| `own_family_gap` | ACS2012–2023 pooled, matched MMS, HHWT heads30–55, standard DUE structures; oldest resident own child under4 minus no resident own children | Ownership for any dependent m>0 minus never-parent n=0, current2023 households in the model age band | Parent and control mismatches; pooled-date seam |
| `own_rate` | Same pooled ACS/DUE/head/age sample, overall ownership | Owner share of2023 model age30–55 band | Pooled-date seam; no family-group contrast |
| `aggregate_mean_occupied_rooms_18_85` | ACS2012–2023 pooled,42MMS metros, HHWT heads18–85, owner/renter with positive rooms | Total occupied housing services divided by total current2023 model household mass | Pooled-date seam; no family-group contrast |
| `aggregate_wealth_to_annual_gross_labor_earnings` | PSID2005–2019; aggregate net worth, heads18–85, divided by aggregate RP/spouse gross labor earnings, ages18–65 | Beginning-of-period wealth divided by annual gross labor earnings on2023 distribution | Older pooled vintage explicitly approved; retain |
| `annual_bequest_flow_to_aggregate_wealth` | Borrowed annual US literature normalization0.0088; not a2023 microdata moment | Annualized at-death post-saving bequest flow divided by aggregate wealth,2023 | External normalization explicitly approved; retain |
| `old_total_wealth_to_annual_income_p90_p50_7684` | PSID1984–2019, living reference persons76–84; weighted p90/median of net worth/annual family income | Beginning-of-period wealth/income dispersion on corresponding2023 model ages | Older pooled vintage explicitly approved; retain |

All empirical ACS ages are exact integer-age restrictions. Current model bands retain the existing age_to_index/whole-node approximation. No rebinning is proposed here. The two ownership rows impose UNITSSTR3:10; the two rooms rows do not. None should be labeled a national household benchmark merely because the model has one market.

## What is already decided, and what is not

The July24 signed review explicitly retains the family ownership gap as a hard, overidentifying row and rejects relegating it to validation. It also approves the3+ versus1–2 rooms row. We must not describe either target as never approved. However, the review calls ownership a matched ACS object and does not explicitly accept the verified recent-parent/all-dependent-parent or no-resident-child/never-parent approximations. The earlier July audit documented the ownership discrepancy; current canonical status still marks reconciliation open.

The same signed review explicitly approves the CPS/NCHS1979–1984 cohort anchor, its seam with the modern environment, and the distinct PSID wealth vintages. Those are not new open decisions. The author has also reaffirmed the expected18-year dependency approximation. This note keeps it fixed and does not reinterpret it as exact child-age tracking.

The two unresolved family objects are precise:

1. **Ownership:** the empirical treated group has ELDCH<4, so every resident own child is younger than four. The model treated group has any m>0, regardless of duration. Empirical NCHILD0 includes empty nesters; model n=0 excludes them.
2. **Family rooms:** data classify the number of co-resident own children of any age, provided at least one is under18. Model m counts dependent units. A data household with children aged22,19 and10 is in NCHILD3+, but the model has no separate adult-child co-residence variable. This is not the old2+ versus1 bin bug: the current literal3+ convention already repaired that different issue.

## Existing date diagnostics: completed, original groups unchanged

The saved aggregate/cache exercise already reports all four original empirical objects for pooled2012–2023 and2023 alone:

| Original ACS object | Pooled2012–2023 | 2023 | Difference |
|---|---:|---:|---:|
| Recent-parent ownership gap |16.766167pp|21.303955pp|+4.537788pp|
| Overall ownership, heads30–55 |57.547241%|58.274640%|+0.727399pp|
| Large-minus-small family rooms |0.367699559|0.435891519|+0.068191961rooms|
| Mean occupied rooms, heads18–85 |5.779970482|5.771329073|−0.008641409rooms|

These are descriptive date comparisons, not adopted target replacements. The current values/weights remain unchanged; no2023 standard errors or2023 weights have been constructed. The family gaps are materially date-sensitive; mean rooms changes little. That evidence does not establish which date window the author should choose.

## Feasible diagnostic correction plan

### 1. Make the measurement contract explicit now

Produce one sidecar manifest from the current canonical provenance ledger: empirical vintage/sample/group/operator, model date/group/operator, approved approximation versus outstanding decision, and input/code hashes. Preserve the original12-row fit and11-parameter table as the main result. Correct misleading labels such as model “new parents” in the sidecar to “current dependent-child households versus never-parents.” This is reporting/provenance work and requires no new economic choice.

Preserve target set `e5_fullhistory_roomsfix_h1_20260817` and fingerprint `3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497`. **The numeric target/weight fingerprint alone does not protect measurement semantics.** Changing the model measurement operator can change the objective even if all target numbers and weights stay identical. Any diagnostic operator must have its own identifier/source hash and must not feed the production objective under the unchanged name.

### 2. Use already available2023 empirical comparisons as a separate sidecar

For a certified-at-its-given-horizon path, display the unchanged2023 model measurements beside BOTH retained pooled and diagnostic2023 empirical observations. Do not replace the loss or optimize against the diagnostic column. This immediately answers the date question using completed data work while leaving group mismatches visible. It is also feasible before a new parameter calibration; horizon qualification must remain attached to the model column.

### 3. Collect directly observable model group aggregates during an already required replay

An optional diagnostic observer can record, at2007/2011/2015/2019/2023, each group's household mass, owner mass and total occupied rooms from `evaluation.g_current`, with renter housing from the matching `hR_pol` and fixed owner-product rooms. Record groups n=0, m=0, m>0, m=1/2, m=3, using the retained age masks. This is aggregation on existing choices: no new Bellman solve, state, parameter or behavioral change is mathematically needed beyond the scheduled path evaluation itself.

In particular, compute a supplemental ownership gap m>0 versus m=0. That includes model empty nesters in the control and exposes the effect of the current never-parent control. Do not call it an exact ACS match: m=0 is no dependent units, while NCHILD0 is no resident own child of any age. The missing adult-child co-residence information remains.

The existing observer currently retains first-birth accounting at dates0–4 but only computes the complete cross-sectional target vector at date4. Saved `transition_path.csv` overall owner_rate is not automatically the target's age30–55 rate. `initial_2023.pkl.gz` stores a PRE-choice distribution, so it cannot substitute for the POST-choice2023 cross-section. Use a dated evaluation callback or the actual matching saved evaluated state; never reconstruct these diagnostics from mismatched timing.

### 4. If a pooled-model comparison is requested, define temporal aggregation rather than average gaps

The model has four-year dates, whereas ACS pools annual observations2012–2023. Exact annual interpolation is not identified by the current saved path. A declared binning or interpolation rule is needed; it must be a separate diagnostic and must reproduce constants and date endpoints. Obtain year/group HHWT denominator totals from the same pinned cache if needed, not from generic person-weighted MMS tables.

For two groups A and B, pooled gaps should be built from pooled numerators and denominators separately:

\[
\Delta^{pool}=\frac{\sum_y W_{A,y}\mu_{A,y}}{\sum_y W_{A,y}}
              -\frac{\sum_y W_{B,y}\mu_{B,y}}{\sum_y W_{B,y}}.
\]

Here W is the declared empirical group/year household-weight total and mu is the model group/year mean under the stated time mapping. An equal-weight average of annual gaps generally does not reproduce the original pooled ACS object. This standardized diagnostic also must state its use of empirical group/year weights; it is not an automatic replacement for a model-predicted population mixture. Do not fabricate annual states between model dates.

### 5. Separate missing information from easy aggregation fixes

- **Recent-parent ownership:** exact ELDCH<4 is not measurable from the current count-only state or a single checkpoint. A passive, explicitly validated cohort/recency accounting extension could describe first-birth cohorts without altering choices, but “first birth in the previous model period” is not automatically identical to oldest current resident own child under4. No exact-match claim is justified without the tracking definition and its initialization.
- **Under18 child counts:** the head-level cache's NCHILD/YNGCH/ELDCH do not generally recover the number of own children under18; neither does the model's geometric dependency count encode literal ages. Restricting data to ELDCH<18 would alter the sample and is not an automatic fix. Adult-child residence and exact ages are genuinely missing information.
- **Empty-nester control:** m=0 is directly observable in the model and can be tabulated immediately, but adopting its gap as the objective still changes measurement and requires an explicit group choice.

## Gates and decision boundary

The first overnight deliverable can be the unchanged full fit plus clearly labeled date/group diagnostics and this provenance map. No row should be dropped, demoted, reweighted or silently replaced. The baseline retains12rows and11free coordinates; row count is necessary, not proof of informative identification. Ownership chiefly restricts tenure/fertility selection jointly; mean rooms disciplines the supply scale; overall ownership disciplines the tenure premium; the family-size rooms gap disciplines the child housing requirement jointly with the first-birth response. A production remapping must preserve those identifying roles and be reviewed by the lead before any new search.

Two substantive choices remain for the author/lead: whether the ACS observations should represent a2023 endpoint or a historical pooled window, and which family/control object should define the hard comparison. Safe diagnostic aggregation can inform those choices; it cannot settle them by relabeling. No renewed decision about the reviewed Sun–Abraham estimator, approved fertility cohorts, approved wealth vintages or maintained18-year approximation is needed for this plan.

## Source anchors and reproducibility

- [Canonical latest status](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/CALIBRATION_STATUS.md:1)
- [Current twelve-row source ledger](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/e5_target_provenance.csv:1)
- [Observer chronology and2023 measurement](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/tools/e5f_matched_pf_moments.py:90)
- [Dated target extraction/wealth timing](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/tools/run_e5f_transition_calibration.py:1452)
- [Joined history and2023 callback](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/tools/run_e5f_matched_pf_history.py:71)
- [Model current-parent selector](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/solver.py:6310)
- [Model ownership parent/control measurement](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/solver.py:6486)
- [Model dependent-family room bins](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/solver.py:7515)
- [ACS ownership definitions](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/audit_ownership_targets.R:67)
- [Exact ACS room masks](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/moment_standard_errors/build_active_acs_room_target_receipt.R:131)
- [Author-approved cohort/wealth vintages](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5_target_review_20260724.md:16)
- [Author-approved ownership hard row](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5_target_review_20260724.md:119)
- [Completed ACS date comparisons/provenance](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/meeting_receipts/acs_date_diagnostic/README.md:1)

Small-file SHA256 pins inspected in this review:

- `tmp/e5f_matched_pf/code/model/tools/e5f_matched_pf_moments.py`: `4ffe5097b50775cfed88978be63804e6804c0889fe49f6609b3c5c00cd41d821`
- `tmp/e5f_matched_pf/code/model/tools/run_e5f_matched_pf_history.py`: `e394ade1e8abeb6291e1e749b029fe81b83cb5aba8eaef82e1c1071d711cd360`
- `tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/e5_target_provenance.csv`: `38e6507fdc81ae264d54ab02d9fb6c824114a3da85cfc67c901eca819bf6e78b`
- `output/model/e5f_matched_pf_20260909a/meeting_receipts/acs_date_diagnostic/acs_date_comparison.csv`: `189c0cf939428b9642a87542cccbf439e51172a0d5ba21ed3b2f058300008139`
