# Intergen Target-Object Audit

This note audits the target objects in the active one-market intergen
calibration before any further recalibration. It is about measurement and
model-statistic definitions, not model fit.

Active reference target set:
`candidate_replacement_roomgap_14moment_tfr192_v1`.

## Main Conclusion

The wealth block has a measurement/statistic mismatch. The target value
`0.17922556` is a PSID annual-income ratio for young childless renters, but it
is currently attached to `young_liquid_wealth_to_income`, a model statistic that
divides by 4-year period after-tax income. The same period-versus-annual issue
also affects the old-age nonhousing wealth-to-income statistics.

The active target sets have not been changed yet. The solver now stores
explicit annual-gross young wealth moments, so the next target revision can be
made deliberately.

## Wealth Target

### Data Object

Source file:
`code/data/psid_followup_mar2026/output/entry_wealth_v1/entry_wealth_candidate_targets_focus_v1.csv`.

Source script:
`code/data/psid_followup_mar2026/entry_wealth_targets_v1.do`.

Relevant PSID objects:

| Sample | Ratio | Mean | Median | Note |
|---|---|---:|---:|---|
| young childless, 25--35 | `NETWORTH2R / INCFAMR` | `0.870076` | `0.141386` | annual family-income denominator |
| young childless renters, 25--35 | `NETWORTH2R / INCFAMR` | `0.179226` | `0.061550` | annual family-income denominator |
| young pre-birth, 25--35 | `NETWORTH2R / INCFAMR` | `0.823997` | `0.166442` | annual family-income denominator |
| young pre-birth renters, 25--35 | `NETWORTH2R / INCFAMR` | `0.204921` | `0.080227` | annual family-income denominator |

Important: the source script reports weighted means but unweighted medians,
because Stata `centile` was used without weights. A median target is sensible,
but it should be re-extracted as a weighted median before becoming a hard SMM
moment.

### Current Model Object

`young_liquid_wealth_to_income` is computed from childless renter states at
ages 25--35 as
\[
\frac{\sum b\,g}{\sum y^{period,aftertax}g}.
\]
With `period_years = 4` and `scale_flows_to_period = True`, this denominator is
4-year after-tax income, not annual PSID income.

### Measurement Error

The hard target attached an annual data ratio to a period-income model ratio.
At the current review point, the old model statistic is `0.353`, while the
same young-childless-renter sample is about `1.10` when expressed against
annual gross income. The PSID mean is `0.179`.

### Right Object

For the current mechanism, the right primary object is:

`young_childless_renter_liquid_wealth_to_annual_gross_income_2535`.

Store both:

- mean: use the PSID mean `0.179226` as the first candidate;
- median: use the PSID median only after re-extracting it as weighted.

Do not keep using `young_liquid_wealth_to_income` as a hard target.

## Active 14-Moment Target Ledger

| Target moment | Target | Model statistic | Data object | Status |
|---|---:|---|---|---|
| `tfr` | `1.920` | `2 * mean_completed_fertility` among post-fertility ages | Should be completed-fertility-equivalent, not period TFR | Target-choice issue. Label is misleading; model object is stationary completed fertility. |
| `childless_rate` | `0.150` | completed parity-zero share among post-fertility ages | Completed childlessness / ever-childless object | Needs source pin. Do not use ACS no-child-in-household for this. |
| `own_rate` | `0.575472` | ownership ages 30--55, `own_rate_3055` | ACS household heads, DUE-housing sample, ages 30--55 | Clean after prior extractor fix. |
| `own_family_gap` | `0.167662` | new-parent minus non-parent ownership, ages 30--55 | ACS household heads, DUE-housing sample, new parent minus no-child household, ages 30--55 | Mostly clean, but data uses current no-child household while model uses never-parent/nonparent states. Flag this distinction. |
| `housing_increment_0to1` | `0.664435` | controlled birth-cohort housing response at model horizon 0 | PSID first-birth room event-study, about 3 years post-birth | Clean after June 25 fix, given 4-year model period approximation. |
| `housing_increment_1to2` | `0.488031` | model additional-child housing proxy, one-child vs two-plus birth states | PSID quick second-birth room change, post-3 | Not a bug, but not a true sequential second-birth hazard in the model. Interpret as additional-child housing demand only. |
| `young_liquid_wealth_to_income` | `0.179226` | young childless renter liquid wealth over 4-year period after-tax income | PSID young childless renter liquid nonhousing wealth over annual income | Measurement error. Replace with explicit annual-income model statistic. |
| `old_parent_childless_nonhousing_wealth_to_income_gap_6575` | `1.007450` | parent minus childless ratio of aggregate nonhousing wealth to period income | PSID parent minus childless weighted mean of annual nonhousing wealth-income ratios | Needs fix. Denominator and mean-of-ratios versus ratio-of-sums are not identical. |
| `prime30_55_childless_renter_mean_rooms` | `3.805288` | weighted mean renter housing services for childless-in-household ages 30--55 | ACS household heads, childless renters, mean `ROOMS` | Clean. |
| `prime30_55_childless_owner_share_rooms_ge6` | `0.596131` | owner share with realized rooms/services at least 6 | ACS household heads, childless owners, `ROOMS >= 6` | Clean after Markov nonlinear-moment fix; economic miss if model fails. |
| `old_nonhousing_wealth_to_income_median_6575` | `2.230461` | weighted median of nonhousing wealth over period income | PSID weighted median of nonhousing wealth over annual income | Needs denominator fix. The model statistic should divide by annual income or target should be converted. |
| `prime30_55_childless_owner_minus_renter_mean_rooms` | `2.418762` | owner mean rooms minus renter mean rooms, childless-in-household ages 30--55 | ACS household heads, childless owner-renter mean room gap | Clean. |
| `old_age_own_rate` | `0.764261` | ownership ages 65--75 | ACS household heads, DUE-housing sample, ages 65--75 | Clean as an ownership target; still an economic/structural miss. |
| `own_rate_2534` | `0.341166` | ownership ages 25--34 | ACS household heads, DUE-housing sample, ages 25--34 | Clean. |

## Recommended Target Changes Before Recalibration

1. Replace `young_liquid_wealth_to_income` with
   `young_childless_renter_liquid_wealth_to_annual_gross_income_2535`.
2. Store the young-renter mean and median as candidate targets, but re-extract
   the data median as weighted before using it in the hard objective.
3. Add annual-income old-age wealth statistics before continuing to target
   `old_nonhousing_wealth_to_income_median_6575` or the parent-childless old
   wealth gap.
4. Rename or document `tfr` as completed-fertility-equivalent in reports. The
   target can be `1.92`, but it should not be described as a period TFR target.
5. Keep the room and ownership targets as the stable block for now; the main
   remaining problems there are economic fit and product mapping, not obvious
   measurement bugs.

## What Changed In Code

Current code now stores these annual-gross diagnostic model moments:

- `young_all_liquid_wealth_to_annual_gross_income_2530`;
- `young_childless_liquid_wealth_to_annual_gross_income_2535`;
- `young_childless_renter_liquid_wealth_to_annual_gross_income_2535`;

each with `_median` and `_mass` variants. These are exposed through
`calibration.extract_moments`. No active target set has been changed yet.
