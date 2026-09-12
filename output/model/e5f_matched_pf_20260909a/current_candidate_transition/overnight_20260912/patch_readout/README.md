# Patch presentation figures

Reproduce the seven-page review and separate slide PDFs with:

`/opt/anaconda3/bin/python -B code/model/tools/build_e5f_patch_readout.py`

Suggested presentation order: fertility model/data history; 2023 equilibrium
(May four-panel layout); 2023 lifecycle comparison; intergenerational allocation
(May six-bar layout). Supplementary: fertility by age, aggregate validation table,
and the finite price/quantity path. All are in `figures/`. Combined review:
`output/pdf/e5f_patch_review.pdf` from repository root.

The patch consists of conditional stationary household fits for the three
four-year fertility windows ending2011,2015,2019, followed by a forward-looking
transition starting from the2019stationary household distribution. It is NOT a
fully carried2007–2023successive-surprise solution. The final preference trial
is0.09239514522037684;2019–23fertility1.5590039164103253 versus1.64575.
All four history targets/values/gaps/preferences are in `figures/historical_fit.csv`.
No data point is plotted at2007: model2.1is its initial normalization.

2023ageprofiles use the actual post-choice2023distribution of that same forward
forecast. Replay17498757 matches all saved dated price, housing, ownership,
birth and pension aggregates exactly (maxabsoluteerror0). Read `source/verification.json`.
The first readout job17498744 stopped before a solve because sequential runtime
initialization was missing; collector initialization was corrected and the exact
replay rerun. Original failure logs remain on Torch.

Housing and PAYGO budgets pass on the finite forecast. Terminal-distance/horizon
checks do not pass, and the historical final window remains unmatched. Prices
and quantities through2039are a finite forecast; the late slope may depend on
the imposed terminal boundary. Do not call it an established long-run response.
The2023person/head totals and age composition are externally conditioned inputs,
not untargeted model predictions. Consumption is in model period units and mean
rooms are physical uncapped rooms in equilibrium; model/data roomcomparisons
cap both at9. Children ever born and dependent children currently at home differ.

ACS2023cross-sectional comparisons use42active metropolitan areas, household
heads,PERNUM1,RELATE1,GQ1/2,HHWT; ownershipcode1,validpositive rooms,noROOMS99
records present. Main ageplots use18–85; aggregate validation andsixbars22–85.
Model agecells correspond to four-year ACS bins. Sixbars split modelcells at40/60
using uniform within-cell weights. Large homes mean6+rooms, owneroccupied.
Children at home means modeldependentcount>0 versusACSresidentownminors
NCHILD>0,YNGCH<18. YNGCH99 staysdenominator as noresidentminor.
This mismatch is economically substantive: model child departure is gradual,
and the late-age discrepancy should not be presented as measured resident minors.

2023profiles andlevel comparisons are NOT calibration targets. Related initial
room/ownership observations were targeted around2005–06; their full scored
contract remains unchanged. InitialcheckpointSHA120ffc45c0fb8756f4182f999c96b7c0236adf315cb938190ec31cd2068c87c2,
not the later e3a4...initialcandidate. The figures do not combine different model candidates.

Fertility-by-age is supplemental and descriptive: model2023–27births annualized
and divided by household exposure, versus2023ACSfemale recent-birth reports/PERWT.
The ACSsample includes ALLwomen18–49 in the same42metros,GQ1/2,FERTYR1/2,
FERTYR2reportingbirth, with nohead/ownership/room restriction. It reports mothers
with recentbirths rather than counting births in multiple deliveries. This is
neither the nationalCDCfour-year historical fertility target nor an exact common
female-exposure measure. Finalextractor was corrected and rerun by the lead after
review caught unintended owner restriction in children-at-home andheadrestriction
in the delegated female-birth draft. Only corrected outputs are plotted.

Sources: `source/` immutable collected aggregate/observer receipts;
`data/` bounded2023ACSextracts andmetadata;2023large-home data from sibling
`age_housing_allocation/comparison_2023/large_owner_data.csv`.
Collector `code/model/tools/collect_e5f_patch_readout.py` performs one replay at
saved coordinates on Torch (190.6seconds), no root search or parameter change.
Extractor `code/model/tools/extract_e5f_patch_validation.py` reads only2023rows.
Slide task owns `latex/september_14_presentation.tex` and incorporates these assets.
