# Completed preference-timing pilot

All three 100-date household/population calculations finished. Later timing improves the inherited twelve-target objective by 13.82%, but the historical birth-count decline remains approximately 21%, versus 11% in the data. This three-case, fixed-price comparison gives a useful direction to investigate; it is not a newly calibrated equilibrium.

| Preference decline | Inherited objective | Maximum market discrepancy | Historical birth-index mean squared gap |
|---|---:|---:|---:|
| Later | 81.423452 | 0.045556% | 0.00933069 |
| Linear baseline | 94.475576 | 0.003673% | 0.00977139 |
| Earlier | 109.029620 | 0.047158% | 0.01024835 |

The market tolerance is 0.020000%. Both changed profiles exceed it. They need a new market-clearing price path before their objective improvements can be assessed as equilibrium results. The inherited terminal unit-rent discrepancy remains 1.088827% against a 1% threshold; horizon extension is still unverified. All household budget, feasibility, distribution and population-accounting checks pass for the saved paths.

The historical birth comparison normalizes each candidate's 2008–2011 aggregate births to 100:

| Birth years | Data | Later | Linear baseline | Earlier |
|---|---:|---:|---:|---:|
| 2008–2011 | 100.00 | 100.00 | 100.00 | 100.00 |
| 2012–2015 | 97.06 | 89.43 | 88.79 | 88.14 |
| 2016–2019 | 93.93 | 83.14 | 82.75 | 82.35 |
| 2020–2023 | 89.04 | 78.78 | 79.07 | 79.36 |

Later timing raises first-block births by 1.4395% relative to the unchanged case; earlier timing lowers them by 1.4300%. Normalizing each case by its own first block removes these level differences, so the complete comparison also reports a common baseline anchor and raw flows. Later timing narrows the middle-block gaps slightly but does not fix the final decline. This is a birth-count diagnostic, not female TFR: historical household totals/age margins are externally conditioned, national births and the model housing geography differ, and additional births in the 3+ bin are imputed at entry into that bin.

The baseline job finished its 100 numerical dates but FAILED its final file comparison. Independent diagnosis found exactly four differences, all the recorded ACS file-location string inside the historical bridge audit. The recorded empirical content hash, every other CSV cell, all prices/residuals, and the complete target-fit, parameter and measurement files reproduce exactly. The original failure is preserved. The collector now verifies this exact, narrowly specified metadata difference while rejecting changed economic numbers or source-content hashes. No numerical gate or model source was relaxed, and no rerun was necessary.

All six smoke/main paths were verified against their 508 source pins and five output hashes per case. The three main paths have 300 verified date observations, zero household budget violations and all 318 recorded numerical/accounting gates passing. This does not turn their separate failed market/horizon checks into passes. Source snapshot e399c90e is unchanged. No pilot jobs remain running; the monitor is paused and no further round has been launched.

[All 36 target-fit rows](computation/all_main_target_fits.csv) include every target, model value, gap, weight and contribution. [All 15 parameter/normalization/fixed rows](computation/all_inherited_parameters.csv) include estimates, bounds and near-bound flags; every parameter was held fixed in this pilot. [Full birth comparisons](computation/all_main_birth_comparisons.csv) retain all normalizations and raw levels. [Verification receipts](computation/verified_receipts.json) and the [baseline metadata diagnosis](computation/baseline_replay_diagnosis.json) preserve the evidence.

The next useful numerical step would be to solve market-clearing prices for the later profile and compare the full tables again. Matching the historical birth decline will also require assessing the preference-change amplitude and the initial calibration; these three fixed-endpoint profiles do not settle that problem.
