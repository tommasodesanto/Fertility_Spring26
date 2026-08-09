# PSID Twins and Same-Sex Housing Re-audit

## Bottom line

The old PSID IV files should not be quoted directly. Their specifications mix
different samples and horizons, and two outcome constructions are defective.
A clean mother-level re-estimation produces imprecise but mostly supportive
directional evidence: the same-birth-year proxy points toward more
space-related moves, more rooms, and more renter-to-owner transitions; the
same-sex design points toward more space-related moves and more ownership and
move-to-owner transitions. None of these weighted reduced forms or 2SLS
estimates is individually statistically significant.

These results are triangulation, not a separate causal pillar. In particular,
same-sex children may share rooms more easily, so sex composition can affect
housing demand directly and is not a clean excluded instrument for a housing
outcome.

## What was re-audited

The historical evidence came from four scripts:

- `twins_and_gender_iv_v1.do`: pre-first-birth renters, outcomes through
  $+3$, indexed at the first birth;
- `secondborn_gender_design_v1.do`: the same renter sample, using mixed sex of
  the first two children;
- `moved_for_size_iv_alltenure_v1.do`: all baseline tenures, moved-for-size
  through $+3$;
- `iv_power_recovery_v1.do`: all baseline tenures, outcomes through $+5$.

These are variants of the classic twins design of Rosenzweig and Wolpin (1980)
and the sex-composition design of Angrist and Evans (1998), applied here to
housing outcomes.

The re-audit uses one observation per mother. The twins design is indexed at
the first birth; the same-sex design is indexed at the second birth, when sex
composition becomes known. Treatments are also horizon-aligned:

\[
D_i^{2+,5}=1\{\text{second child arrives by first birth}+5\},
\qquad
D_i^{3+,5}=1\{\text{third child arrives by second birth}+5\}.
\]

The sample requires a complete pre-birth baseline no more than four years
before the event and an interview at event time $+4$ or $+5$. Birth cohorts
must permit the full five-year calendar window by 2019. The same-sex design
excludes same-year multiple-birth proxies at the first and second birth. The
regressions use the baseline PSID individual weight and robust standard errors,
and control flexibly for age, education, birth cohort, baseline tenure, and
baseline rooms; the same-sex design also controls for first-child sex.

The available twins variable is not a verified clinical indicator. It is a
same-birth-year sibling proxy, so it can include rare closely spaced births.

## Corrected estimates

The reduced form is the coefficient of the instrument on the housing outcome.
It is the safest directional object when the first stage is weak. Binary-outcome
coefficients below are probability changes; for example, `0.043` means 4.3
percentage points.

| Design | Outcome | $N$ | $Z=1$ count | First stage | $F$ | Reduced form (SE) | $p$ | 2SLS (SE) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Same-birth-year first-birth proxy | Moved for size by $+5$ | 1,974 | 34 | 0.367 | 61.8 | 0.034 (0.093) | 0.717 | 0.092 (0.252) |
| Same-birth-year first-birth proxy | Change in rooms by $+5$ | 1,974 | 34 | 0.367 | 61.8 | 0.212 (0.316) | 0.503 | 0.578 (0.852) |
| Same-birth-year first-birth proxy | Ever owner by $+5$, baseline renters | 974 | 17 | 0.390 | 36.1 | 0.038 (0.122) | 0.755 | 0.098 (0.303) |
| Same-birth-year first-birth proxy | Moved to own by $+5$, baseline renters | 923 | 16 | 0.374 | 36.5 | 0.019 (0.108) | 0.859 | 0.051 (0.282) |
| Same sex of first two | Moved for size after second birth | 1,835 | 916 same-sex | 0.072 | 7.2 | 0.007 (0.025) | 0.769 | 0.104 (0.347) |
| Same sex of first two | Change in rooms after second birth | 1,836 | 917 same-sex | 0.073 | 7.4 | -0.089 (0.069) | 0.197 | -1.229 (1.053) |
| Same sex of first two | Ever owner after second birth, baseline renters | 961 | 465 same-sex | 0.123 | 9.3 | 0.034 (0.040) | 0.396 | 0.275 (0.326) |
| Same sex of first two | Moved to own after second birth, baseline renters | 912 | 446 same-sex | 0.129 | 9.9 | 0.043 (0.030) | 0.146 | 0.335 (0.249) |

For the twins design, the moved-for-size, rooms, and move-to-own reduced-form
signs remain positive without weights; the ownership sign flips. For the
same-sex design, all four signs survive the unweighted sensitivity: positive
for moved-for-size, ownership, and moved-to-own, and negative for rooms.

The first stage is strong for the same-birth-year proxy, but there are only 34
such mothers in the broad housing sample and 16--17 among baseline renters.
That rarity explains why a high first-stage $F$ does not deliver precise
housing reduced forms. The same-sex first stage is moderate, with weighted
$F$-statistics of roughly 7--10 depending on outcome availability.

## Problems in the historical files

1. `twins_and_gender_iv_v1.do` defines `own_post3_i` as a Boolean expression
   even outside the event window. In the matched raw panel, all 18,619 parents
   without an observed ownership value in $0{:}3$ are consequently assigned
   zero rather than missing.
2. The move-to-owner code uses Stata's calendar lag `L.own`. The PSID becomes
   biennial after 1997, so this is not the previous observed interview. Among
   11,313 comparable parents, the calendar-lag rule finds 210 transitions; the
   previous-observation rule finds 634 and changes 424 classifications. The
   historical $+5$ twins estimate of `-0.029` (SE `0.011`) for move-to-own is
   therefore not a valid result to cite.
3. `moved_for_size_iv_alltenure_v1.do` uses the reason-for-move code without
   confirming that a move occurred. This changes 93 of 15,064 comparable
   classifications.
4. The old same-sex outcome window begins at the first birth, before the sex of
   the second child is realized. The corrected design starts at the second
   birth.
5. The old files include both parents for household-level outcomes and use
   final recorded family size against short housing windows. The corrected
   design uses mothers and births occurring within the matching five-year
   window.
6. `ACTUALROOMS_` contains non-room codes `0` and `99`. The re-audit treats
   values outside 1--20 as missing; 6,712 of 704,258 matched parent-years have
   such codes. This coding issue should be checked separately in every active
   rooms-event-study builder before another rooms target is treated as final.

## Defensible interview language

> I also use two traditional fertility instruments as a triangulation exercise.
> A same-birth-year first-birth proxy has a strong first stage, and same-sex
> first two children increases the probability of a third child. The housing
> estimates are noisy, but the corrected reduced forms generally point toward
> additional children increasing space-related moves and transitions into
> ownership. I do not use their magnitudes as targets or present them as a
> separate causal pillar: twins are rare in the PSID, the same-sex first stage
> is only moderate, and sex composition can directly affect room-sharing needs.

If asked for one number, the most favorable corrected directional result is
that same-sex first two children predicts a `4.3` percentage-point increase in
moving from renting to owning within five years after the second birth (SE
`3.0` points, $p=0.146$). It should be described as suggestive and
statistically imprecise.

## Reproduction

Run:

```bash
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -b do \
  /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/iv_housing_reaudit_20260809.do
```

Machine-readable outputs are `iv_housing_reaudit_estimates.csv`,
`iv_housing_reaudit_samples.csv`, `iv_housing_reaudit_balance.csv`, and
`iv_housing_reaudit_construction.csv` in this directory.
