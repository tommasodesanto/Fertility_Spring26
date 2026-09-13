# Successive surprise history: recovered final-window candidate

Actual household history is carried through the three fitted windows into 2019. Each shock is unexpected, and households forecast its preference level to persist. The final candidate below passes all finite-path housing, pension, household and replay gates; its fertility error exceeds the 0.005 acceptance tolerance. Search was stopped by the author on September 12 at22:48EDT. The terminal-distance check fails, and no horizon certificate or production transition-policy result is claimed.

| Window | Target | Model | Model − target | Preference | Fertility fit accepted |
|---|---:|---:|---:|---:|---|
| 2007–2011 | 1.974875 | 1.974856 | -0.000019 | 0.147087157 | Yes |
| 2011–2015 | 1.861000 | 1.861069 | +0.000069 | 0.136087157 | Yes |
| 2015–2019 | 1.755375 | 1.755536 | +0.000161 | 0.125555932 | Yes |
| 2019–2023 | 1.645750 | 1.633313 | -0.012437 | 0.111555932 | No |

These are scalar shock-fitting roots, not a new weighted structural SMM calibration: the structural parameters and their original target/weight contract remain unchanged. Every historical window and its preference proposal is shown. The final-window relative miss is -0.7557%.

Reproduction: root_receipt.json records exactly zero market/fiscal/reproduction differences; source/expected_transition.csv and source/fertility.json are the saved native aggregate path and fertility measurements. The prior three accepted windows are in source/realized_fit.json. Six four-year forecast dates are used per surprise.

Cluster job17559194: candidate_path_20260911a/batches/finite_sequences_20260912/recover_saved_2019. The author stopped this and the other no-rebate runs at22:48EDT. No further trial is running; the intended baseline must rebate property-tax revenue. Standard policy diagnostics are retained in each admissible candidate’s accepted_graphs folder.


The continuation now uses the saved 2019 surprise forecast, with preference constant after its final historical change. Model period fertility is1.6262,1.6106,1.5997,1.5993,1.6287 in windows ending2027,2031,2035,2039,2043. These values have no long-horizon certificate; the final upturn must not be interpreted as a robust prediction.

## 2023 Data/Model readout

One fixed-coordinate cluster replay, job17575193, completed in248.5seconds and reproduced the entire saved path with maximum absolute difference0. No search was run. `source/readout_2023/verification.json` and `measurement_verification.json` record exact replay and successful observers. The inherited2019state is pinned by SHA256 `7787c71ae5278ef41bade3dd797b92f07db782b29ffe502277ed3d819c627a5f`; the approved initial checkpoint remains `120ffc45c0fb8756f4182f999c96b7c0236adf315cb938190ec31cd2068c87c2`.

`figures/validation_2023.csv`, `.tex`, `.pdf` and `.png` contain all13 initial-calibration moment families. Rows are untargeted comparisons at the carried2023model state. Empirical sources are CPS2024 (nearest fertility supplement), NCHS2023, ACS2023, explicitly pooled PSID samples, and the external bequest benchmark. The CSV preserves exact definitions, gaps and vintages. It does not silently label pooled wealth or event-study estimates as2023observations.

Completed fertility is children ever born at40–44, not period fertility or the initial2.1normalization. First-birth timing compares model2023–2027flows with annual2023birth counts on the same model age-cell midpoint mapping. First-birth rooms use the dated2019–2023 matched branch. National household/maternal proxies, the42-metro ACS sample, stochastic dependents as proxies for resident children, and model pension income versus PSID family income remain the same measurement limitations documented in the observer JSON; no new targets or weights were adopted.

The cluster replay driver and exact arguments are recorded in `/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/carried_readout_20260913/submission.json`. The full2023state remains there in `source/state_2023.pkl.gz`; only small readout artifacts are copied locally.

Regenerate the figures, full table and two-page PDF (no model solve):

```sh
/opt/anaconda3/bin/python -B code/model/tools/build_e5f_patch_readout.py --base output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/patch_readout_fit --sequence-base output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/recovered_sequence --pdf output/pdf/e5f_current_history_and_2023_fit.pdf
```

## Completed-fertility and housing continuation

At the author's request, one fixed-coordinate replay observed every date of the
same saved2019forecast. Job17580030 completed in2m55s (evaluation164.53s),
with exactly zero aggregate-path discrepancy. All2023 profiles, birth-flow
measurements and fertility-stock observations match the original2023packet
exactly. No parameters, targets, equilibrium prices, demographic inputs or
fiscal rules were fitted or changed. The no-rebate/horizon qualifications above
remain in force.

`source/stock_forecast/observed_dates.json` preserves the dated observations and
receipt hash; the accompanying verification and submission files identify the
remote inputs and source hash. The2023match and independent weighted-count
fertility calculation are checked by the plotting driver. The first launch,
17579982, was cancelled after11seconds to repair a NumPy-array comparison in
report validation; its scientific inputs were unchanged and its output was
never admitted.

`figures/historical_fit_stock_forecast.png` and `.pdf` add the completed-fertility
and capped-room continuations to the existing period-fertility comparison.
`figures/stock_forecast.csv` gives every dated value. Completed fertility means
children ever born at ages40–44 using the same uniform-birth-time observer as
2023:1.6951 in2023,1.6285 in2027,1.5620 in2031,1.5085 in2035 and1.4703 in2039.
It is a near-completion measure, not fertility observed at age50. Its date is
the observation year; period fertility is plotted at the end of its four-year
birth window and therefore extends to2043. The first preview began in2019; the full-history extraction below supersedes
that incomplete preview. Initial
stationary profiles must not be relabeled as historical2007observations.

The completed-fertility data are Census CPS Historical Table2, not Goldin or HFD.
The source URL, workbook hash and extraction definition are in
`source/stock_forecast/cps_source.json`; recent2022/2024 counts are capped atfive.
Housing data retain the42-metro ACS capped-nine-room definition, while the model
uses its existing household/demographic normalization. The plotted housing
measure is rooms per household, not total housing stock.

Reproduce without solving:

```sh
python3 code/model/tools/build_e5f_stock_forecast_comparison.py
```

For a fresh identical numerical replay, the existing collector now accepts
`--all-dates` with the pinned carried-state arguments recorded in
`source/stock_forecast/submission.json`. It records every observer date but
retains the original aggregate-path reproduction gate. Large states remain on
Torch. The plot is supplemental; the presentation deck has not been edited.

## Full historical stock and fertility comparison

The complete figure now includes model observations at2007,2011,2015,2019,2023,
2027,2031,2035,2039 in every panel. Three exact fixed-coordinate replays recovered
only the first, realized date of the admitted2007/2011/2015forecasts. Their later
expectations were not spliced into the realized history. Jobs17581096,17581099,
17581100 passed with exactly zero discrepancy against every saved forecast row;
evaluation times245.29,151.68,242.48seconds. The prior2019–2039observations and
2023 verification remain intact. Source hashes, admitted-vintage receipts and
pinned inherited states are in `source/historical_stock/`.

Housing is now TOTAL occupied physical rooms, capped atnine per household and
summed across households. The plot compares national ACS and model stock indices,
each2007=100; it does not validate the initial housing level. The national series
matches the geographic scope of the model's demographic conditioning. The42-metro
sample remains the initial cross-sectional calibration source; no target changed.
Both national and42-metro raw totals, household counts and definitions are retained
in `source/historical_stock/housing_data.json`. A single file-backed ACS scan took
21.36seconds and passed24 exact checks against previously saved aggregates. It did
not rerun bootstrap estimates or scan the source into memory in full.

National stock growth2007–2023 is18.2448%, model growth4.3911%. The model's
completed-fertility measure is1.8962,1.8645,1.8250,1.7611,1.6951 over the five
historical dates. These are untargeted historical comparisons conditional on the
calibrated parameters and imposed demographic inputs; they show material misses.
National ACS restricted-sample household counts and imposed CensusHH-3 counts
are not identical, so the stock gap must not be attributed wholly to endogenous
housing behavior. Vacancy is outside the model and this empirical stock measure.

### Separate constant-rate completed-fertility illustration

At the author's request, `figures/completed_fertility_constant_rates.png/pdf`
shows a separate mechanical scenario, outside the presentation. The source is
Driscoll and Hamilton (NCHS,2025),
[Table2](https://www.ncbi.nlm.nih.gov/books/NBK617829/table/nvsr74-3.t2/)
for annual age-specific rates1990–2023 and
[Table4](https://www.ncbi.nlm.nih.gov/books/NBK617829/table/nvsr74-3.t4/)
for independently published TFRs. The exact transcribed rates and TFRs are in
`source/stock_forecast/nchs_asfr_1990_2023.csv`. Reconstructed TFRs agree within
0.006 for all34years (published rates/TFRs have different rounding precision).

For cohort c, sum annual rates at ages10–49 in calendar years c+a, using the
historical rate when c+a<=2023 and the2023 age-specific rate thereafter. Rates
are uniform within each published five-year age group; the45+ group is assigned
five ages45–49, following the published TFR convention. This is an approximate
cohort reconstruction from grouped period rates, not an exact single-age Lexis
calculation or the CPS survey measure. No migration-selection or individual
birth-history response is modeled. Constant TFR alone would not identify the
scenario: the entire age schedule is held constant.

Projected completed fertility is2.1864 for the1980cohort,1.9336 for1990,
1.6610 for2000 and1.6210 for2010. The eventual level equals the2023 reconstructed
TFR1.621 by accounting identity. These estimates must not be joined to the
CPS ages40–44 line: the source, population accounting, age horizon and measure
differ. The figure's bottom axis is year of reaching50; its top axis is birth
cohort. Every cohort-year-age cell and the historical/projected decomposition
are saved in the two companion CSVs. Independent grouped sums, all-period TFR
checks, the constant-schedule identity and plotted values pass verification.

Reproduce with:
`MPLCONFIGDIR=/tmp/psid_correction_review/matplotlib python3 code/model/tools/build_e5f_completed_fertility_scenario.py`.
This is an illustrative conditional calculation, not a published forecast,
economic-model prediction or new calibration target. The deck is unchanged.

**Two-panel dotted-extension variant:** at the author's subsequent request,
`figures/fertility_introduction_with_projection.png/pdf` preserves every point
of the1980-onward introductory graph and appends dotted scenarios. Reproduce with
`python3 code/model/tools/build_e5f_fertility_introduction.py --projection`.
The default historical figure and presentation PDF are not overwritten.

The completed-fertility extension is rebuilt for ages40–44; it does not splice
the age50 birth-rate calculation above onto CPS survey means. It starts at the
exact2024 CPS value1.918 and advances the observed2024 younger-age means in
five-year cohorts. Source: Census2024 Tables1 and3a, with extracted values,
source hashes and URLs in `source/stock_forecast/cps_2024_age_profile_scenario.json`.
Ages20–39 use normalized, rounded number-of-children shares (five-plus coded5);
ages15–19 use the reported Table3a mean to avoid suppressed cells. Source
rounding and initial topcoding limit precision. Future expected births are
added to these means without reapplying topcode5: this is an illustrative
mean-stock extension, not an exact forecast of the topcoded survey statistic.
Migration, selective survival and future cohort-size changes are omitted.

The future age pattern is NCHS2023, rescaled by0.997223936 to match exactly the
existing WDI2023 period value1.6165. Thus the period dotted line has no source-
change jump. Single ages are represented at their midpoints with uniform weights
inside each five-year age band. Age exposures are integrated over the fixed
schedule; a separate half-year exposure sum reproduces every increment.
Cohorts under15 in2024 use negligible early births from the same fixed schedule.

The ages40–44 mean is1.9095 in2029,1.7538 in2034,1.6456 in2039 and converges to
1.57985. Its limit is below lifetime period fertility because it omits births
after the ages observed. The small2044–2054 variation is retained from the
initial CPS age profile, not smoothed into a forced monotone curve. Complete
scenario rows, all input hashes, exact historical-artist checks and definitions
are in the companion CSV and verification JSON. This remains a separate
illustration for review, not an economic-model result or slide replacement.

**Presentation reveal adopted:** the author subsequently requested putting the
historical graph and dotted extension on successive overlays of the same frame,
with visible sources removed. `--reveal` now exports
`fertility_introduction_reveal_history.pdf/png` and
`fertility_introduction_reveal_projection.pdf/png`. They share identical axes,
positions and historical artist values, so the second overlay only introduces
the scenario curves, endpoint labels and scenario legend. The September14 frame
Fertility in the United States uses these two images with Beamer `alt` overlays
and a fixed caption area in a top-aligned frame.
The first caption introduces the2007 initial-state approximation; the second
explains falling completed fertility across successive cohorts at fixed
age-specific rates. Source citations and construction limitations remain here
and in the saved receipts; the visible source/footer text is removed from the
presentation images and frame. Values and projection assumptions are unchanged.
Reproduce with `python3 code/model/tools/build_e5f_fertility_introduction.py --reveal`.

### Housing and population comparison, September13

**Presentation simplification, September13:** the author subsequently removed
housing from this introductory slide. The deck now uses
`figures/fertility_introduction.pdf`, two historical-data panels (annual period
fertility1980–2023 and CPS children ever born ages40–44 from1980 through2024), immediately
before Calibration and Historical Transition. It states that2007 is approximated
by an initial steady state, rather than claiming the graphs establish stationarity.
The housing and model-comparison packets below remain research diagnostics.
Reproduce the introductory figure with
`python3 code/model/tools/build_e5f_fertility_introduction.py`.
The saved WDI response and existing CPS history are its pinned inputs; the figure
receipt verifies both plotted arrays exactly. CPS2022/2024 counts are capped atfive.
Both deck PDF copies are updated from a twice-compiled, visually inspected build.
The author requested extending the start from1990 to1980 to show older cohorts.
The additional CPS observations were checked directly against Historical Table2,
rows34–42: children ever born were2.988 in1980,2.447 in1985 and2.147 in1988.
These are survey dates for women aged40–44 (1980 corresponds approximately to
birth cohorts1936–1940), not birth-cohort labels. The refreshed WDI response adds
1980–1989 with every previously plotted1990–2023 value unchanged. Both panels
share the expanded1.4–3.15 vertical range; no historical observations are clipped.

`figures/housing_population_comparison.pdf` (and PNG/CSV) provides a supplemental
historical comparison using the same national ACS housing households and all
their resident person records. It leaves the presentation deck unchanged.
Following the author's correction, every panel now compares model and data:
household counts, total occupied rooms, rooms per household, and rooms per
resident at common observed demographics. The first three include the full
saved model history and continuation2007–2039; data end in2023. The fourth ends
in2023 because the common observed denominator is unavailable beyond that date.
Top-row model aggregates retain their native Census demographic conditioning,
with its different household-count coverage disclosed in the figure. The former
data-only decomposition panels are removed. Verification checks that all four
panels contain both model and data and preserves the exact plotted native series.
Households and capped rooms use HHWT; residents use PERWT, following the
[IPUMS household](https://usa.ipums.org/usa-action/variables/HHWT) and
[person-weight](https://usa.ipums.org/usa-action/variables/PERWT) definitions.
The sample still requires a head aged18–85, valid owner/renter tenure, positive
rooms and non-group-quarter residence. It is not the entire US population.
Complete consecutive person rosters and every prior household/room total were
checked; source hashes and accounting checks are saved in
`source/historical_stock/housing_population_data.json` and the figure receipt.

Between2007 and2023, sample residents grow10.9514%, households16.4814%, total
rooms18.2448%, and rooms per resident6.5735%. The ratio of weighted residents to
weighted households falls from2.6364 to2.5112; this is distinct from directly
averaging roster size with household weights. As a weighting sensitivity,
HHWT-weighted roster population implies rooms per resident growth3.8924%.
That sensitivity is retained in the source and comparison CSV, not substituted
silently for the person-weighted population estimate.

The retained model has no separate resident-person history before2023:
`code/model/tools/e5f_successive_surprises.py` carries household distributions
through the historical dates and introduces the fixed person-population anchor
in2023. The2023 resident-person/head ratio2.6265 includes a different population
universe and must not be described as matched ACS household size. Historical
person counts cannot be reconstructed by adding a partner and dependent count.

The bottom-right model line is explicitly standardized to observed demographics:
model mean rooms per household times ACS households, divided by ACS residents.
Its2007–2023 change is−2.9494%, versus data+6.5735%; it is not a model population
prediction, a new equilibrium, or an unconditional historical fit. Model mean
rooms per household still falls7.5569%, versus data rising1.5138%, while its
level remains above the national data. At common household counts, model total
rooms grow7.6790%, rather than the native aggregate4.3911%. Initial42-metro
calibration targets and model states remain unchanged. This comparison isolates
housing intensity; it does not diagnose its economic cause. Native per-person
historical forecasts are not supplied where the necessary population is absent.

Reproduce from the repository root:

```sh
/opt/anaconda3/bin/python code/data/Spatial_aggregate_withmicrodata/build_housing_population_comparison.py
MPLCONFIGDIR=/tmp/psid_correction_review/matplotlib python3 code/model/tools/build_e5f_housing_population_comparison.py
```

All panels now use the same model dates. Period fertility is dated at the START
of its four-year birth window: e.g. the2019point summarizes births2019–2023.
This is an explicit plotting change, not a change in any numerical fertility
estimate. Every plotted period-flow number is checked against the previous
figure receipt. The separate old-state2.1normalization is no longer plotted as
an observed2007period point. The2039flow covers2039–2043; stocks are observed at
2039. The same age40–44 birth-time projection and CPS source are maintained.

The complete figure replaces the former single-panel historical-fertility image
in `latex/september_14_presentation.tex`, frame `Fertility and Housing`.
Rebuild it with the same plotting command above. The historical replay source is
`code/model/tools/collect_e5f_historical_stock_observers.py`; exact arguments are
in `source/historical_stock/submission.json`. The empirical builder is
`output/model/e5f_matched_pf_20260909a/design_research/housing/build_historical_stock.py`;
run with `/opt/anaconda3/bin/python -B` to regenerate `/tmp/full_housing_stock_history.json`.
The no-rebate baseline and unresolved horizon limitations remain unchanged.
