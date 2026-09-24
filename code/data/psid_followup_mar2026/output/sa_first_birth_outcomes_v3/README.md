# First-birth event studies for ownership and moving (version 3 sample)

Date: 2026-09-24. Same specification as the rooms headline
(`../sa_rooms_first_birth_v2/README.md`): Sun–Abraham in two-year interview
windows, baseline −3/−2, person and survey-year fixed effects, age and
education covariates, clustered by person, PSID weight IW, first biological
birth, confirmed-childless controls, entry at first adult observation. Two
designs: H (women who are reference person or spouse, one per single-family-unit
household-year) and A2h (all adults who were reference person or spouse in the
baseline window). Torch jobs 18428887 (smoke) and 18428888 (eight fits), all
pass; covariances symmetric positive definite; receipts reproduce coefficients.

Items rebuilt from the official year-specific PSID family variables via
`../../psid_family_item_crosswalk.csv` (41 waves; verification JSON and the
verbatim code lists are in `../sa_rooms_first_birth_v2/`). The shelf's
MOVEDFREF_ and WHYMOVED1_ were verified to carry the same one-interview
displacement as rooms and are not used; HOMEOWN was verified contemporaneous
and agrees with the rebuilt ownership on 826,417 rows with zero mismatches.

Definitions. `own`: owns the dwelling (1) vs rents (0); "neither" missing.
`moved`: moved since the previous interview (1 yes / 5 no; the reference
period is "since spring of last year" through 2001 and "since January 1 of the
prior year" from 2003). `moved_space`: moved with first-mentioned reason
"expansion of housing: more space, more rent, better place" (code 3); zero if
did not move or moved for another coded reason; missing if reason DK/NA.
`moved_nbhd`: same for reason "neighbourhood-related" (code 6). Each outcome
uses its own complete-case sample.

| Outcome, design | Baseline mean | −1/0 | +1/+2 | +3/+4 | +5/+6 | Rows |
|---|---|---|---|---|---|---|
| Owns, H | 0.41 | +0.11 (0.01) | +0.16 (0.02) | +0.16 (0.02) | +0.14 (0.02) | 61,229 |
| Owns, A2h | 0.43 | +0.16 (0.01) | +0.25 (0.01) | +0.27 (0.01) | +0.26 (0.01) | 111,272 |
| Moved since last interview, H | 0.58 | −0.04 (0.02) | −0.12 (0.02) | −0.16 (0.02) | −0.16 (0.02) | 63,817 |
| Moved since last interview, A2h | 0.54 | −0.08 (0.01) | −0.17 (0.01) | −0.21 (0.01) | −0.21 (0.01) | 118,875 |
| Moved for more space, H | 0.070 | +0.014 (0.009) | +0.034 (0.009) | +0.020 (0.010) | +0.007 (0.010) | 62,813 |
| Moved for more space, A2h | 0.064 | +0.015 (0.006) | +0.021 (0.006) | +0.010 (0.006) | +0.009 (0.006) | 116,697 |
| Moved for neighbourhood, H | 0.032 | −0.004 (0.006) | −0.005 (0.006) | −0.010 (0.005) | +0.002 (0.006) | 62,813 |
| Moved for neighbourhood, A2h | 0.030 | −0.009 (0.004) | −0.008 (0.004) | −0.009 (0.004) | −0.003 (0.004) | 116,697 |

Full paths in `window_path.csv`; figure `first_birth_outcomes_v3.png/pdf`.

Reading. Ownership rises by 16 (H) to 27 (A2h) percentage points on a base of
about 42%, with the A2h pre-period flat (0.00, 0.00, −0.02) and the jump
concentrated in the birth window and the two years after. The probability of
having moved since the last interview rises into the baseline window (A2h:
−0.13, −0.10, −0.03 before it), then falls by 16 to 21 points after the birth:
the household moves around the birth and then stays put. Moves motivated by
more space rise by 1.5 to 3.4 points on a 6–7% base in the birth window and
the two years after, then return to baseline. Neighbourhood-motivated moves do
not respond. Together with the rooms path this is one event: a move to a
larger, owned dwelling between two years before and two years after the first
birth, followed by a decade of stability.
