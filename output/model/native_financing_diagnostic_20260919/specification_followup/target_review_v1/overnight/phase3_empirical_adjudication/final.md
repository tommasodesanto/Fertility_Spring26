# Measurement choices for the calibration targets after the empirical checks

This memo reads the four worker packets (fertility, entry wealth, old-age wealth, first-birth rooms), the target catalogue and the sources verified earlier. It changes nothing: every recommendation is a paper choice for the author, and tonight's experiments stay frozen.

"B" means the saved selected case from the earnings/entry battery, a finite-search diagnostic, not a calibration. Numbers marked "derived" are my arithmetic on saved receipts.

## 1. Recommended package

### Choices that can be made now

1. **Fertility level: keep 2.1 as an explicit normalization, not a target, and use the matching vintage for the 3+ group.**
   - 2.1 is the stationary benchmark's replacement closure. It pins ψ (the fertility utility scale) and the births-to-entrants renewal law, and it approximates the 2003–06 period rate of 2.06.
   - The representative number of children in the 3+ group should come from the same survey as the stock rows: 3.4828, the June 2004/06 capped mean, instead of the June 2024 value 3.6024.
   - Add two untargeted validation rows on late fertility (see §2, row 0). The economic cost of this choice is stated there.
2. **First-birth timing: keep the current convention.** Data and model are binned into the same four-year model cells, with births under 18 assigned to the first cell and births over 45 to the last. The mean-age target stays 25.976 and the share at 30+ stays 0.2493. The paper must say this is a model-cell mean, not the NCHS published mean (25.0 in 2006).
3. **Old-age dispersion: switch to the dispersion of wealth levels,** $Q_{.90}(W)/Q_{.50}(W)=4.9194$ (SE 0.5798), on the identical PSID 2003/05 sample. This is a new estimand needing a new contract name.
4. **First-birth rooms: keep the statistic, fix the model window.** Keep $\hat\beta_{+3}-\hat\beta_{-1}=0.7202$. Do not substitute $\hat\beta_{+3}=1.2017$, which is measured against the omitted −2 year. Align the model observer to a single housing decision (§5).
5. **Bequests: keep 0.0088 as a borrowed external restriction,** stating its Gale–Scholz construction. Do not rescale it.
6. **Aggregate wealth / gross earnings: keep** the July definitions.
7. **Mean rooms, ownership 30–55, family-room slope and recent-parent ownership gap: keep the current objects.** The family-room slope and the recent-parent gap are descriptive associations, not causal effects.
8. **Entry wealth: replace both current conversions with one common earnings scale and an empirical wealth–earnings association.** The two being replaced are the July "ratio × own model income" rule and B15's preserved wealth levels. Keep the July 18–24 childless-renter sample (§4).
9. **Current-income purchase eligibility (Y/R) stays the stated timing assumption.** The accounting is reconciled; no change is proposed here.

### Choices conditional on a comparison or an author decision

- **A. Economic population for the four ACS rows.** My default is national, unless the paper's claims are about large-metro housing markets.
- **B. Primary first-birth window.** My default is the one-decision window. It could be reversed if the author wants the target to include early second births.
- **C. Keep the CPS cohort stock rows under the 2.1 normalization, or move to period-consistent stock measures.** This depends on the late-fertility diagnostic.

## 2. Decision table for all 13 rows

| # | Current object | Exact model counterpart | Demonstrated issue | Recommended paper treatment | Parameters and identifying information |
|---|---|---|---|---|---|
| 0 | **2.1 normalization (unscored).** Author-chosen. Completed fertility over all households aged 46+ (population-weighted), with the 3+ group valued at 3.602359422009 (June 2024 CPS capped mean). | The same statistic. ψ is solved so it equals 2.1 within 5e-4. The same 2.1 divides births into entrants in the renewal law. | Two issues. (i) The 3+ value comes from a different vintage than the stock rows (3.4828 in 2004/06). (ii) The two sides differ in age window, not arithmetic. In B, the 40–44 projection has mean 1.883 (at T = 3.602; CPS at the same T gives 1.891). The last fertile cell (42–45) then raises the mean from 1.863 to 2.100, adding 0.237. It holds 2.5% of B's first births, against 0.54% of NCHS first births at 42+. About 10% of women enter the 3+ group at 42–45 (derived). | Keep 2.1, labelled as a stationary replacement closure. Use T = 3.4828, which requires slightly more 3+ mass to reach 2.1. Report two validation rows: the 40–44 distribution against its value at completion, and the share of first births at 42–45 (model vs 0.488%). | ψ is fixed. The late fertility tail is not directly targeted. |
| 1 | Childless women 40–44 = 0.198279. CPS June 2004+2006, supplement weights, n = 10,872 (reproduced exactly). | Stock of the model's household reproductive member, projected to ages 40–44 assuming births fall uniformly within cells. B: 0.1689. | None in measurement. It is a 1960–66 cohort stock used in a stationary model (declared). | Keep, presented as the cohort stock for the early economy. | κ_fert (entry noise) and first-birth cost; ψ enters through the normalization. |
| 2 | Share with exactly one child among mothers = 0.213655, same sample. | Same projection. B: 0.2553. | None. It is a stock share, not a second-birth hazard. | Keep. | κ_fert continuation and first-birth cost. |
| 3 | Mean age at first birth = 25.976264. NCHS 2003–06 counts (6.61m first births) mapped to four-year cell midpoints, with boundary collapse. | Stationary first-birth flow weighted by cell midpoints (20, 24, …, 44). B: 27.296. | None; it is a convention. The raw mean is 25.161, or 25.661 with +0.5; the teen collapse adds 0.250. The published 2006 mean of 25.0 is a different definition. | Keep. Footnote the published mean. On support: the first model cell must hold 34.2% of first births, where B places 24.2% (derived). | κ_fert and first-birth cost, conditional on the external fertility schedule. |
| 4 | Share of first births at 30+ = 0.249278 (exact age ≥ 30 over all first births aged 12–49). | Share of first-birth flow in cells starting at 30 or later. B: 0.30195. | None. The teen collapse does not affect it. Dropping births under 18 instead would move it to ≈0.270, and the CPS stock rows would then need the same exclusion, so that option is rejected. | Keep. | The timing tail; same parameter block as row 3. |
| 5 | Aggregate wealth / annual gross earnings = 6.145861 (bootstrap SE 0.3629). PSID 2003/05: family net worth of reference persons 18–85 over reference person + spouse earnings at 18–65 (EARNINDRRC), ratio of sums. | Beginning-of-period Σ(b + pH) over Σ annualized gross earnings in working cells. B: 6.059. | None demonstrated. The shelf labels EARNINDRRC as combined tax-year earnings; the builder treats it as gross labor earnings. PSID excludes employer DC accounts, so asset coverage is a separate question. | Keep. Disclose the tax-year/survey-year dating and the DC exclusion. | β primarily (July), with θ₀ and the housing block through pH. |
| 6 | Annual bequests / wealth = 0.0088. External: Gale–Scholz IRP DP 1019-93, Table 4, via De Nardi–Yang (2014), Table 2. Working scale 0.00044 (synthetic 5%). | Σ over deaths of max(b′ + pH, 0), divided by 4 years, over beginning-of-period aggregate wealth; death certain in the final cell. B: 0.00633. | The concepts differ, and this is not a bug. Gale–Scholz count only the children's share of mortality-imputed estates (25% for married heads, 75% for single heads or both spouses), exclude pensions and trusts, and treat life insurance separately. The model counts every household estate, with no marital structure. | Keep as an external restriction in the table of borrowed targets. Describe the construction in one sentence. Label the scale a tolerance, not a standard error. No rescaling. | θ₀ (July), with θ₁ and β. |
| 7 | Old-age dispersion: $Q_{.90}/Q_{.50}$ of wealth / family income for living PSID reference persons aged 76–84, 2003/05 = 3.515935 (SE 0.3069). | p90/p50 of (b + pH) divided by retirement income. In B, retirement income is a constant 0.51159 across all selected cells, so the statistic is exactly the wealth-level p90/p50, 3.72771. | Demonstrated mapping mismatch (not a builder bug). The model denominator is constant; the data denominator varies. On the same 644 observations the level statistic is 4.919 (SE 0.580). | Replace with the level statistic, 4.919429 (SE 0.579804), under a new contract name. B would then miss by 1.19, which is 2.05 SE. The model counterpart needs no new solve. | θ₁ (July), with θ₀. Earnings dispersion and entry wealth enter as inputs, not free parameters. |
| 8 | Mean occupied rooms, capped at 9 = 5.561097 (42 metros; metro bootstrap SE 0.0884). National value 5.607886. | Realized rooms capped at 9 before aggregation, with fractional age exposure (verified). B: 6.245. | None. Geography is a choice. | Keep the object; population per decision A. | H₀, the supply scale (July). |
| 9 | Ownership, heads 30–55 = 0.648334 (42 metros; structure types 3–10). National value 0.676260. | Owner mass at 30–55, counting half of the 54–57 cell. B: 0.4533. | None. 2005/06 is near the historical peak; Sommer–Sullivan target a long-run 0.65 instead. That precedent is secondary next to a 19.5 pp miss. | Keep; descriptive level; population per decision A. | χ, the owner premium (July), with H₀ and β through prices. |
| 10 | First-birth room response: $\hat\beta_{+3}-\hat\beta_{-1}=0.7202462624$. Covariance-based, person-clustered SE 0.0852600513. Sun–Abraham estimator, confirmed never-treated controls, person and year fixed effects, weight IW. | Current: clone each pre-choice state into birth and held-childless branches; compare realized rooms after the next date's decisions, with second births allowed. B: 1.3415. | The window is misaligned by construction: the model difference spans two housing decisions (birth date and next date). Counterfactuals and weights differ (declared). No coding bug. | Keep the statistic. Primary observer: the difference after the birth-date decision only. Report the current two-decision value and the second-birth shares as diagnostics (§5). | Floor arm: h_P. Share arm: δ_jump + δ_a (July logic), with first-birth cost, χ and H₀. |
| 11 | Rooms, 3+ versus 1–2 resident own children (heads 30–55 with youngest child under 18) = 0.347067 (SE 0.0597). NCHILD counts own children of any age, including step- and adopted children. | Current dependents 3+ versus 1–2, capped at 9. | A grouping approximation; no bias has been shown. | Keep as a descriptive cross-sectional association. Improvement I2 below. | Share arm: δ_a (July: the jump cancels among parents). Floor arm: no dedicated parameter (per-child rooms fixed at 0). |
| 12 | Recent-parent ownership gap = 0.162896: households whose oldest child is under 4 minus households with no resident children (42 metros; SE 0.00608). National value 0.127608. | Current births from previously dependent-free homes minus current empty homes, including former parents; synchronized snapshot at 30–55. | None. The repaired observer stands; it approximates an interview stock with a flow snapshot. | Keep the definition; descriptive, not causal; population per decision A. | No dedicated parameter. An overidentifying restriction on how fertility and tenure interact (first-birth cost, the κ's, χ, h_P or δ_jump). |

## 3. Why these calls

### Normalization (rows 0–4)

- **What the data do and do not show.** B does not over-predict the 3+ group at 40–44 (27.0% against CPS 28.6%). B's own 40–44 mean at T = 3.602 (1.883) is close to the CPS mean at the same T (1.891).
- **Where the 2.1 comes from in B.** The 0.217 gap between 1.883 and 2.1 accrues after the 40–44 window, mostly in the 42–45 cell. The NCHS share of first births at 42+ is 0.54%. The receipts contain no counterpart for higher-order births.
- **Economic cost.** A stationary benchmark at replacement, combined with a cohort that completed about 1.88, has to be reconciled somewhere. B reconciles it through a late-fertility tail that no scored row penalizes directly.
- **Identification cost.** ψ is fixed by the normalization. The stock and timing rows then fall on κ_fert, κ_fert continuation and the first-birth cost. The late tail is disciplined only indirectly, through the mean age and the 30+ share.
- **Why not a cohort-consistent level (≈1.86–1.89) instead.** It would break replacement stationarity and the period-based 2007–2023 path, so I keep it as a robustness check only.
- **The honest disclosure** is that the benchmark's completed fertility exceeds the stock cohort's by design (the period–cohort gap), shown with the two validation rows.

### Timing convention (rows 3–4)

- The binning is shared by data and model, which is what makes the comparison valid. Excluding births under 18 would be less coherent, because those women are in the CPS stock of mothers.
- The support fact that matters is shape: the first cell must carry 34.2% of first births (B: 24.2%) and the last cell 0.54% (B: 2.5%). B's timing miss spans the whole support, not an artifact of the collapse.

### Old-age dispersion (row 7)

- The July rationale for a scale-free statistic carries over unchanged to the level ratio.
- The ratio version asks a model with constant pensions to match a statistic compressed by variation in family income. The data drop from 4.92 to 3.52 when divided by family income; which income component causes this is not identified.
- The cost is precision: SE 0.580 against 0.307.
- The alternative would be pensions that depend on earnings history. That is a model change and is not proposed here.

### Bequests (row 6)

Without marital structure, the direction of the concept gap is ambiguous. The model's full estate at household termination exceeds Gale–Scholz's 75% for single decedents, but Gale–Scholz also count 25% of married decedents' estates, which have no model event. Disclosure beats an invented adjustment.

## 4. Entry-wealth construction (an input, not a scored row)

**What the checks established.**
- The July 18–24 sample and its five ratio nodes reproduce exactly.
- Dividing by each household's own earnings is unstable: on common rows the mean of ratios jumps from 0.286 (family income) to 0.912 (gross earnings) because of small denominators. Medians (0.108 against 0.135) and ratios of sums (0.372 against 0.483) are different objects.
- B15 preserves wealth levels exactly, but its implied median ratio is 0 (the empirical median is 0.0985). Its imposed wealth–income rank correlation is 0.098, against 0.280–0.282 for gross earnings in the PSID and 0.330 for family income.

**Recommendation: one economy-wide scale plus a measured association.**
1. For each PSID row, divide nonhousing net worth by $\bar E_t$, the wave-t weighted mean reference person + spouse gross earnings at ages 18–65. This is the same concept as the row 5 denominator, applied per wave to remove real growth.
2. Rank entrants into weighted earnings terciles, with zero earners included.
3. Within each tercile, compute weighted quintile-bin means of that normalized wealth: 15 nodes.
4. Model side: multiply the nodes by the model's mean annual gross earnings at working ages, and assign terciles to entrant income states by rank overlap.

**What this buys.**
- Wealth and earnings are measured in the same units as the aggregate wealth target.
- No household-specific small-denominator ratios.
- No rows are dropped: all 1,835 are used, where requiring earnings above $1,000 drops 93.
- The observed association is kept (about 116 family-years per cell).
- No zero-entry fallback and no transfer. Mass that is infeasible under the model's borrowing limits goes through the unchanged frontier-censor rule and is reported.

**Caveats.** It changes the mean entry wealth in model units, with sign and size unknown until computed. This is an improvement over declared approximations, not a bug fix.

## 5. First-birth rooms: statistic versus observer

**Why keep the statistic.**
- The saved coefficients are $\hat\beta_{+2}$ = 0.726 (derived from the saved mean, $\hat\beta_{+3}$ and $\hat\beta_{+4}$), $\hat\beta_{+3}$ = 1.202 and $\hat\beta_{+4}$ = 0.831. The +3 coefficient sits about 0.42 above its even-year neighbours; the correction review reports an average odd-year peak of 0.420.
- In the biennial PSID era, odd and even event years are observed for different birth-year cohorts. So $\hat\beta_{+3}$ measured against the omitted −2 compares cohort sets: 3,297 of 9,309 observations at +3 come from cohorts never observed at −2.
- Both four-year contrasts within the same odd/even phase agree: +3 against −1 gives 0.720 and +2 against −2 gives 0.726.
- The omitted −2 is the estimator's normalization, needed to identify the event-time coefficients. It does not decide which contrast is scored.

**Why change the observer window.**
- The current observer differences housing after two decisions: the birth-date choice and the next date's choice. The next-date choice also includes that date's aging, income changes and second births.
- The primary observer should be the treated-minus-control realized housing after the birth-date decision only: the housing carried into the next date, before its decisions. That is one model period and one decision, matching a four-year change around the birth.

**What stays different.**
- The data window contains a second child for 50.4% of mothers by +3 (60.7% by +4). The one-decision model window cannot contain second births.
- Control groups and weights differ.
- The current observer's inclusion of next-date second births may partly offset the first point, which is why the choice of primary window remains conditional (question B).

## 6. Utility arms and identification

**Parameter counts.** The floor arm searches nine structural coordinates (β, κ_fert, κ_fert continuation, χ, H₀, θ₀, θ₁, first-birth cost, h_P); the share arm searches ten (δ_jump and δ_a replace h_P). ψ is normalized in both.

**What the rooms rows can inform.**
- On July's logic, first-birth rooms loads on δ_jump + δ_a, the family slope on δ_a, and mean rooms on H₀. So in the share arm the rooms rows can inform both share terms.
- They are not the only parameters these rows respond to. Prices (H₀), tenure (χ), who has children through the κ's and the first-birth cost, the fixed equivalence scale, and the discreteness of housing products all move the same rows.
- In the floor arm, the family slope has no dedicated parameter.

**What cannot be claimed.** Local rank remains unverified without a Jacobian at each arm's own point. The share arm has one more free coordinate, so a lower loss is partly mechanical. Adaptive search paths diverge after the common anchors, which confounds attributing fit differences to the utility form. Compare the arms row by row; no moment is proposed for dropping.

## 7. Implementation order (no solves needed except where noted)

### Fixes: the model counterpart currently measures a different object

**F1. Old-age level row.**
- Input: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_old_wealth/wealth_results.csv`.
- Observer: weighted p90/p50 of beginning-of-period b + pH (b for renters) over the 76–84 overlap weights 0.5/1.0/0.75, with no income division, in `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/e5f_initial_housing_observer.py`.
- Output: a new contract row, 4.919429 with SE 0.579804.
- Check: B returns exactly 3.7277115.

**F2. First-birth window.**
- Code: in `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_transition_calibration.py`, realize both branches' current cross-section at the birth date with the origin-date policies, using the same mass gates.
- Output: the one-decision response; the current two-decision response; the treated second-birth mass at the next date as a share of origin mass, shown next to the PSID 50.4%/60.7% (49.5%/59.5% for strictly later birth years).
- Checks: the two-decision value reproduces 1.341457 at B, and treated and control masses are equal.
- Needs a saved evaluation (policies and pre-fertility distribution); no equilibrium solve.

### Improvements: the current mapping is an acceptable approximation

**I1. Entry construction (§4).**
- Inputs: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_entry/replicate_entry_units.R` (sample) and the per-wave 18–65 earnings totals from the aggregate builder.
- Output: 15 nodes with person-cluster bootstrap (499 draws, seed 20260715), plus the model mapping produced through `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/e5f_earnings_wealth_contract.py`.
- Checks: the pooled node mean equals the data ratio of sums; the weighted rank correlation matches (descriptive); report the censored mass.

**I2. Family-room child count.**
- If `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/housing/inspect_early_housing.py`'s extract carries a household identifier, count the householder's children under 18 (RELATE = child, AGE < 18). Otherwise keep NCHILD.
- Check: first reproduce 0.3470669 with NCHILD.

**I3. Late-fertility validation rows.**
- From `.../empirical_fertility/empirical_reproduction.json` (B's post-fertility shares by age) and `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/nchs_natality_timing/first_birth_counts_year_age.csv`, tabulate first births by model cell for model and data.
- If the local natality cache keeps all live-birth orders, add the share of all births at ages 42–45.
- Check: B's cell-42 mean reproduces 2.0999.

**I4. Top-group value and national rows (author decision).**
- Set T = 3.4828 (a normalized solve is needed later).
- If national is chosen, compute the national family-room slope with the same filters, and decide how the four ACS rows are weighted.

## 8. Claims withdrawn or narrowed

- **Withdrawn and contradicted by data:** that B's 3+ group is overpredicted at 40–44 and that part of the CPS-row loss is forced. B's 3+ share at 40–44 is 27.0% against the CPS 28.6%.
- **Withdrawn:** "β(−1) ≈ 0.056, little anticipation". That came from the May curve. The current $\hat\beta_{-1}$ is 0.4815, most plausibly the odd/even cohort-set alternation; this is an inference, not proven.
- **Narrowed:** "the old-age denominator mismatch is not binding". Under the level statistic, B misses by 2.05 SE.
- **Narrowed:** "asset income compresses the data ratio". The data show compression from 4.92 to 3.52; which component of family income causes it is not identified.
- **Narrowed:** "B's entry ratios differ from July's". The pooled mean is close (0.254 against 0.259); the median (0 against 0.0985) and the rank association (0.098 against 0.280) differ.
- **Withdrawn:** my phase-1 suggestion to build entry nodes as each household's wealth/earnings ratio. The long upper tail (mean 0.912) rules it out.
- **Still withdrawn:** the "0–8-year" window, the claim that the range of coefficients bounds the model estimand, and the claim that the housing floor causes the first-birth room gap.

## 9. Open questions that could change a recommendation

1. **A.** Does the paper claim to represent the U.S. population or large-metro housing markets? A metro claim would keep the 42-metro rows and relabel the fertility and wealth rows as national approximations.
2. **B.** Should the first-birth target include second births within three years? If the author says yes, and the model's second-birth share at the next date matches the PSID's 50–61%, the current two-decision window becomes the better counterpart.
3. **C.** Is the benchmark's late-fertility tail acceptable? If all-order data show negligible fertility at 42–45 while B's tail persists, the period–cohort gap is better handled with period-consistent stock measures than with 2.1 plus CPS cohort stocks.
