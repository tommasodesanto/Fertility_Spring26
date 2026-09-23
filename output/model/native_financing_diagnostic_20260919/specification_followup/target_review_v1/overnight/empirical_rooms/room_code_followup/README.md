# PSID room-code follow-up (bounded read-only check)

## Finding

The 100 aligned room-code zeros in the completed Sun–Abraham `e(sample)` are not, by codebook definition, generic missing values. PSID codebooks describe zero as “None; [family unit] shares room.” It is a substantive report of no room assigned exclusively to the family unit. Treating these rows as missing would alter the outcome/sample and would need a separately specified estimand. The single surviving aligned value 98 must be interpreted by interview wave: the codebook treats 98 as an actual count in older waves but as “DK” from 1994 onward. The current builder recodes 98 from interview year 1994 onward and keeps 98 in 1985–1993, matching those definitions.

The builder shifts `ACTUALROOMS_` from the preceding PSID shelf row forward to the interview year it represents, then applies year-specific missing-code rules and requires a 1- or 2-year row gap. That makes the aligned interview year the applicable codebook wave for interpretation; the preceding row’s calendar year is only the field-storage row. This check did not produce a cell-level cross-tab of surviving code 0/98 by interview year, source year, event time, reporter role, or cohort, so it does not establish the specific year/cohort pattern of the 100 zero values or identify the lone 98 observation.

## Saved numerical evidence and limitation

The completed replay’s `room_code_diagnostics.csv` reports exact `e(sample)` aggregate counts: 49,457 rows, of which 100 have aligned raw room value 0, one has value 98, and none has value 99. A prior singleton-pruning reconstruction matches aggregate sample counts (49,457 rows, 4,112 people, 1,654 confirmed childless controls; 415 person singletons removed and no singleton years), but the individual `e(sample)` marker is not saved. Counts alone cannot establish row identity.

I made one targeted, selected-column local data-read attempt to produce the requested grouped counts, but it stopped during in-memory first-birth construction before producing any table. I did not retry the raw read. Therefore there is no new cell-level grouped evidence here, and no codebook-driven recode was applied to the target or any regression.

## Measurement proposal for a separately authorized test

Retain code 0 as a valid outcome value for the current “rooms available to the family unit” measurement, subject to the author deciding that shared-room observations belong in that outcome. Keep code 98 handling wave-specific: preserve it only in waves where the official codebook defines it as an actual count; set it missing where codebooks label it DK. Any separately proposed exclusion of zero-room reports should be described as a sample/estimand change and evaluated in a new, explicitly labeled empirical test. No such test or target replacement was run here.

## Sources

- Repository builder: `code/data/psid_followup_mar2026/sa_rooms_first_birth_household_aligned_v1.do` (shift, year-specific recodes, eligibility and row-gap restrictions).
- Local wave/variable crosswalk: `code/data/psid_followup_mar2026/output/first_birth_correction_review/all_wave_variable_crosswalk.csv` (1968–2019 room-variable names and labels).
- Official PSID [1968 Family Codebook](https://psidonline.isr.umich.edu/documents/psid/codebook/fam1968_codebook.pdf), variable V102 (“HOW MANY ROOMS”): 0 is “None, shares room”; 1–8 are counts; 9 is NA/DK.
- Official PSID [1985 Family Codebook](https://psidonline.isr.umich.edu/documents/psid/codebook/FAM1985_codebook.pdf), variable V11614: 0 is “None; FU shares room”; 1–98 are actual room counts; 99 is NA/DK.
- Official PSID [1993 Family Codebook](https://psidonline.isr.umich.edu/documents/psid/codebook/fam1993_codebook.pdf), variable V22425: 0 is “None; FU shares room”; 1–98 actual count; 99 NA/DK.
- Official PSID [1994 Family Codebook](https://psidonline.isr.umich.edu/documents/psid/codebook/fam1994er_codebook.pdf), variable ER2029: 0 is “None; FU shares room”; 1–97 actual count; 98 DK; 99 NA/refused.
- Official PSID [1997 Family Codebook](https://psidonline.isr.umich.edu/documents/psid/codebook/fam1997er_codebook.pdf) and [1999 Family Codebook](https://psidonline.isr.umich.edu/documents/psid/codebook/fam1999er_codebook.pdf) retain 0 as shared-room/no-room and 98 as DK.

## Follow-up v3 saved-artifact review

The bounded v3 Stata census completed its sample assertions and wrote both requested aggregate exports before stopping during the later metadata receipt export. Its log verifies 49,457 retained rows and 4,112 retained people in the iteratively singleton-pruned reconstruction, plus 100 aligned raw-zero rows, one aligned raw-98 row, and zero raw-99 rows. These are aggregate-equivalent singleton-pruned support counts; the individual regression `e(sample)` row marker is unavailable, so this review does not claim row-level identity with `e(sample)`.

I checked the saved cell CSV against the totals CSV: the 63 grouped cells sum to 100 rows for code 0 and one row for code 98, and those sums match `suspect_code_totals.csv`. No code-99 group is present, consistent with the logged zero count. The synthetic grouping smoke artifacts cover expected codes 0, 98, and 99. The run then failed at `singleton_pruned_confirmed_controls` with Stata `r(198)` because that metadata variable name exceeds Stata's 32-character limit. Therefore only the two aggregate exports are complete; the full execution and metadata receipt are failed/incomplete. See `reviewed_completion_receipt.json` for artifact hashes and the precise boundary.

The v3 script hash matches the recorded correction plan. Its diff from v2 changes the output directory and counts grouped observations from `ID`, which addresses the prior collapse failure. No regression or estimator command is present in the reviewed census script. The aggregate exports are usable for the requested descriptive follow-up, subject to the row-identity limitation above. This review makes no economic interpretation or code-recoding decision; the supplied project instruction assigns those decisions to the lead.
