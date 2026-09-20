# PSID comparison reference for ACS housing event diagnostics

This note records the verified saved PSID readouts that can contextualize the
ACS first-birth and coresident-second-birth diagnostics. These are external
comparators, not updates to the ACS calibration target and not causal
validation of the ACS proxy.

| Saved readout | Sample, event clock, and reference | Outcome, design, and uncertainty | Authoritative source |
|---|---|---|---|
| Corrected first birth, primary target | 49,457 estimation observations, 4,112 individuals; current women age 18+, reference person or spouse/partner, positive `IW`, single-FID dwelling, one woman per household-year. First biological birth is taken from `RELCHI1`--`RELCHI20` type/year records. Bins run from `F7` (at or before −7) through `L11` (after +10); reference is −2. | `ACTUALROOMS_` is shifted forward one observed interview within person to align the rooms item to the preceding row. Person and survey-year FE, age/education controls, `IW` probability weight, ID clustering; confirmed zero-child women are controls. Saved (+3-(-1)) contrast: 0.7202463, SE 0.0852601; normal-approximation CI [0.55314, 0.88736] is derived from the saved contrast and SE. | [`target_receipt.csv`](../../../data/psid_followup_mar2026/output/sa_rooms_first_birth_household_aligned_v1/target_receipt.csv); builder [`sa_rooms_first_birth_household_aligned_v1.do`](../../../data/psid_followup_mar2026/sa_rooms_first_birth_household_aligned_v1.do) |
| Second birth, all saved variant | 203,238 ID-year observations, 22,185 IDs. The script keeps adults with an observed first-child year, drops same-year first/second and second/third births, and retains periods after the first birth. It does not restrict to women or mother-only records; `SEX` is carried but not filtered. Second-birth event bins run from `F7` (at or before −7) through `L11` (after +10); reference is −2. | `ACTUALROOMS_` is renamed directly to `rooms` in this legacy script. Survey-year FE plus age/education controls, `IW` weighting, ID clustering; controls are IDs that remain at one observed child. At +3: 0.5658138, SE 0.1277770, derived normal CI [0.31537, 0.81626]. At +5: 0.5760105, SE 0.1262027, derived normal CI [0.32865, 0.82337]. | [`rooms_s_c_y_all_summary.csv`](../../../data/psid_followup_mar2026/output/sa_rooms_second_birth_with_onechild_controls_v1/rooms_s_c_y_all_summary.csv); builder [`sa_rooms_second_birth_with_onechild_controls_v1.do`](../../../data/psid_followup_mar2026/sa_rooms_second_birth_with_onechild_controls_v1.do) |
| Second birth, restricted diagnostic | 76,655 ID-year observations, 9,725 IDs from the same longitudinal construction, additionally requiring first-to-second gap ≥5 years and no observed third birth by +3. Reference remains −2; event bins remain at or before −7 through after +10. | Same direct `ACTUALROOMS_` outcome, controls, `IW`, clustering, and year/age/education specification as the all variant. At +3: 0.8451701, SE 0.3061890, derived normal CI [0.24504, 1.44530]. At +5: 0.2203672, SE 0.2494201, derived normal CI [−0.26850, 0.70923]. | [`rooms_s_c_y_no_third_by3_gap5_summary.csv`](../../../data/psid_followup_mar2026/output/sa_rooms_second_birth_with_onechild_controls_v1/rooms_s_c_y_no_third_by3_gap5_summary.csv); same builder as above |

The first-birth readout is closer to the ACS first-birth comparison because it
uses a single-family-unit household-year and explicitly aligns the rooms
measure to the interview. The second-birth readouts are only directional
context for the ACS coresident second-birth proxy: they use longitudinal
second-child birth-year records, one-child controls, an ID-year sample that is
not mother-only, and a legacy direct rooms measure. The ACS proxy instead links
coresident children in repeated cross-sections and infers event time from the
second-oldest linked child's age. The ACS rooms outcome requires a predeclared
harmonization after the extract27 code audit because the raw source contains
values above 9; literal comparisons to uncapped PSID rooms bins must use that
documented outcome definition.
The five-year-gap PSID restriction is relevant to ACS negative-time support,
but does not make the populations or estimands identical. The saved PSID
second-birth tables do not establish a causal effect or validate the ACS proxy.

All confidence intervals above are normal-approximation intervals calculated
from the saved coefficient and standard error; the second-birth summary CSVs
save coefficients and SEs but not CIs. The first-birth target receipt likewise
saves the contrast and SE; its component event-study table saves coefficient
CIs, not a contrast CI.
