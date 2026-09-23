# ACS resident-child mapping check

This bounded diagnostic measures how the existing ACS `family_rooms` target changes when the child groups count only linked resident own children under age 18. It uses the same 2005/2006 pooled ACS head sample, household weights, room cap, and 42 MET2013 metros as `housing_profiles_v1`; it does not change a live target, weight, model, or geographic choice.

## Evidence and method

The canonical source is `extract27.dta`, identified by the housing-profile provenance receipt. The local housing reader and `target_recomputed.json` already reproduce the four housing moments. This script independently rereads the two relevant year ranges with a memory-mapped file and preserves the original filters:

- General housing sample: sample code `year*100+1`, group quarters in codes 1 or 2, person number 1 and relationship-to-head code 1, positive `HHWT`, head age 18–85, owner/renter status 1 or 2, positive rooms, and membership in the same 42 metros.
- Ownership and recent-parent moments additionally require head age 30–55 and `UNITSSTR` 3–10. The recent-parent gap is owners among `NCHILD>0` and `ELDCH<4` minus owners among `NCHILD=0`.
- `family_rooms` additionally requires head age 30–55, `NCHILD>0`, and `YNGCH<18`; it has no `UNITSSTR` restriction. `HHWT` weights all means and group contrasts, and rooms are capped at 9 before averaging.

The candidate count uses the child records' ACS parent pointers: count a resident child once if either `MOMLOC` or `POPLOC` equals the head's `PERNUM` (1); split that count at age 18. The all-age linked-child count, capped at 3, matches the householder's `NCHILD` bin for all 275,009 households in the `family_rooms` sample. This verifies capped bins only; because `NCHILD` is top-coded at 3, exact all-age counts above 3 cannot be validated from it. The `RELATE=3` relationship category is a second check; it differs from `NCHILD` in 1,739 households (0.779% of sample weight), so it is not used for the candidate.

## Results

All four existing ACS moments reproduce to machine precision: `mean_rooms=5.561097`, `ownership_30_55=0.648334`, `family_rooms=0.347067`, and `recent_parent_ownership=0.162896`. The source housing receipt's active42 values are the comparison values.

Among the 275,009 relevant family heads (total household weight 29,286,406), 40,450 have at least one parent-linked resident child age 18 or older. Their weighted share is 14.952%. The existing all-age `NCHILD` grouping yields 0.347067 rooms for 3+ minus 1–2 children. Replacing the grouping with 3+ versus 1–2 parent-linked children under 18 yields 0.336220 rooms, a paired shift of −0.010847 rooms. The candidate groups contain 55,266 and 219,743 households, respectively; no households in this sample fall outside the two groups.

The estimate is a new target candidate only. There is no bootstrap in this pass, and the small point shift does not establish that either grouping is adequate. In the model, `family_rooms` maps current dependent count `m` (subject to `m<=n`) to the empirical group. This diagnostic does not reconstruct the model's child ages, departure process, or age-specific dependent stock, so it measures an ACS grouping approximation without validating that model-dependent interpretation.

## Files

- `recompute_child_mapping.py` — deterministic raw-source reader and calculation.
- `housing_profiles_v1_child_mapping_recomputed.json` — aggregate sample receipt and point estimates.
- `target_reproduction.csv` — four exact reproduction checks and the candidate-only row.
- `family_rooms_group_aggregates.csv` — household counts, weighted totals, and capped-room means for the four old/new groups.
- `input_hashes.json` — raw-source identity reused from the canonical receipt, upstream reader/target hashes, and this script's hash.

No person- or household-level records are written. The raw file was not rehashed in this pass; its size and modification time match the canonical housing receipt, whose stored SHA-256 is recorded in `input_hashes.json`.
