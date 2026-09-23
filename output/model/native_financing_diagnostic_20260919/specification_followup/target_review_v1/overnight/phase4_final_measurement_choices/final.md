# Final measurement choices for the calibration targets

Nothing below is adopted, and all running targets and bundles stay as they are. "B" is the saved finite-search point from the earnings/entry battery, not a calibration.

## 1. Recommendation to the author

**Supported now, with no new contract needed:**
- **2.1 fertility normalization.** Keep it for the current tests. The production contract requires completed fertility to equal 2.1 and uses 2.1 to convert births into entrants. Disclose it as a replacement-level closure, not an empirical estimate. That is a property of this contract, not a theorem that a stationary model cannot use another level.
- **First-birth timing rows.** Keep the matched four-year cell convention: mean age 25.976, share of first births at 30+ 0.2493. Say explicitly that 25.976 is not the raw NCHS mean.
- **Wealth / gross earnings 6.146.** Keep the July definitions and the exact current fiscal and wealth units.
- **Bequests / wealth 0.0088.** Label it a borrowed restriction. The concept gap is real: Gale–Scholz count only children's shares of mortality-imputed estates, while the model counts every positive household estate. Its 5% scale is a synthetic tolerance, not a measured match.
- **Mean rooms, ownership 30–55, recent-parent ownership gap.** Keep the current objects. The population they describe is decision 1 below, not a fit choice.

**Changes that would need a named new contract:**
- **Old-age wealth dispersion.** Replace the ratio with p90/p50 of wealth levels: 4.919 (SE 0.580) on the same sample.
- **Family-room slope.** Candidate that counts only linked resident children under 18: 0.336220 (the current all-age grouping gives 0.347067). It has no bootstrap yet.
- **First-birth room response.** A contrast that holds the set of cohorts fixed. It is contract-eligible only after a standard error exists (§3, item 1).
- **Entry wealth (an input law, not a target).** A common-scale construction (§3, item 4).
- **Value assigned to the 3+ group.** 3.4828 in place of 3.602 is conditional on decision 3. It is an economic change, because ψ and the birth flows move when the model is recalibrated. It is not a metadata cleanup.

**Choices that depend on the paper's population or estimand:**
1. **Population for the four ACS housing rows: national or 42 metros.** This follows the paper's claim about whom it represents, never which fits better.
2. **The first-birth room estimand.** Choose the reference window and cohort weights (§2, row 10). Stay with the author's −2 reference unless the author explicitly moves.
3. **The fertility vintage and population.** If the early economy is meant to be the 2004/06 CPS cohorts, the 3+ value 3.4828 follows and the late-fertility diagnostic becomes part of the disclosure. Otherwise keep 3.602 and state the vintage gap.

## 2. Disposition of all 13 rows

| # | Row | Disposition | Rationale (one sentence) | Parameters informed |
|---|---|---|---|---|
| 0 | 2.1 normalization (unscored) | **Retain as a normalization**; 3+ value conditional on decision 3 | Required by the production contract; B's late-fertility tail is a finite-point diagnostic (below), not proof of impossibility. | ψ, and the entrant conversion |
| 1 | Childless share, women 40–44 | **Retain** | Reproduced exactly (n = 10,872); a cohort stock under the declared 40–44 uniform-birth-time observer. | κ_fert and first-birth cost; ψ through row 0 |
| 2 | Exactly one child, among mothers | **Retain** | Stock share, not a hazard. | κ_fert continuation and first-birth cost |
| 3 | Mean age at first birth | **Retain** | A matched binning convention: births under 18 add 0.250 years, and the raw mean is 25.161. | Timing block: κ_fert, first-birth cost, with the external fecundity schedule |
| 4 | Share of first births at 30+ | **Retain** | The threshold sits on a cell boundary and the teen collapse does not affect it. | Same as row 3 |
| 5 | Wealth / gross earnings | **Retain** | July definitions preserved; PSID's source-coverage note alone does not show a mismatch with a stated model asset concept. | β, with housing wealth through pH |
| 6 | Bequests / wealth | **Retain, labelled borrowed** | Concept gap and synthetic tolerance disclosed; no rescaling. | θ₀, with θ₁ and β |
| 7 | Old-age p90/p50 | **Proposed replacement** | The model's retirement income is constant (0.51159), so its statistic already is a level dispersion (3.7277). The data ratio divides by family income that varies across households; on identical observations the level statistic is 4.919. | θ₁ and θ₀; inherited earnings dispersion enters as an input |
| 8 | Mean rooms (capped at 9) | **Retain; population open** | Reproduced; cap and age aggregation verified. | H₀ |
| 9 | Ownership, heads 30–55 | **Retain; population open** | Descriptive level. | χ, with H₀ and β through prices |
| 10 | First-birth room response | **Measurement not settled** | 0.7202 depends on cohort normalization constants: a witness normalization gives 0.7588 with identical fitted values. Common-cohort, fixed-weight candidates are 0.811 and 0.798, with no SE. | Floor arm: h_P. Share arm: δ_jump, δ_a. Both: first-birth cost, χ, H₀ |
| 11 | Rooms, 3+ versus 1–2 children | **Retain; candidate replacement** | The all-age grouping includes the 14.95% of heads with a child aged 18+; counting only children under 18 shifts the slope by −0.0108. Neither grouping is shown equivalent to the model's dependent-exit process. | Share arm: δ_a directly, δ_jump indirectly. Floor arm: no dedicated parameter |
| 12 | Recent-parent ownership gap | **Retain; population open** | The repaired observer stands; a descriptive association, not a causal effect. | Overidentifying restriction on how fertility and tenure interact |

**Late-fertility tail in the saved B point (row 0).** Within the age-42 cell, the mean number of children rises from 1.8635 to 2.1000. Third births account for 9.98% of that cell's mass. First births at 42+ are 2.53% of model first births, against 0.536% in NCHS 2003–06. No local all-order NCHS series exists to compare higher-order births.

**Row 10: why the estimand is unsettled.** Each birth-year cohort's event coefficients are identified only up to a cohort-specific constant. Because cohort support differs across endpoints (34 cohorts at −1, 32 at +3), the pooled contrast inherits those constants. Differences taken within 27 common cohorts cancel them. That repair changes composition and weights, and the choice of pre-birth or post-birth weights matters. Moving to the −2 reference is itself a change of estimand, because cohort support at −2 differs.

**Identification, direct versus indirect.** In the share arm, δ_a raises desired rooms for each additional dependent, so it moves the family slope directly. δ_jump enters every parent's utility. But chosen rooms are endogenous, and parents of 1–2 and of 3+ children differ in income, wealth, age, tenure, product constraints and selection into family size. A common jump therefore shifts their rooms by different amounts and does not cancel. Prices and h_P, which binds for some parents and not others, add further indirect channels. The two utility classes are not nested, so an extra coordinate does not guarantee lower loss. None of this is rank identification; that needs Jacobians at each arm's own point.

## 3. Next-work order

**1. Covariance for the first-birth room contrasts.**
- *Input:* `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/sa_rooms_first_birth_household_aligned_v1.do`, on the frozen sample.
- *Operation:* re-run this existing data regression, with authorization, exporting the full vector and covariance of the cohort-by-event coefficients. The saved folder holds only aggregated estimates plus one covariance term; the Torch 20260912 covariance files belong to a different, original specification.
- *Output:* the current aggregate; the common-cohort −1→+3 contrasts under pre-birth and post-birth weights; a common-cohort contrast referenced to −2 for cohorts with −2 support. Each gets a delta-method SE.
- *Pass:* exact reproduction of 0.7202462624 and 0.0852600513, and of 0.811397 and 0.797759.

**2. First-birth model timing diagnostics (no solve).**
- *Input:* the saved stationary evaluation for B and for each overnight arm, when collected.
- *Operation:* using the existing branch functions, record (a) treated-minus-control realized housing at the birth date and (b) the current destination measure, plus the treated branch's second-birth mass at the destination.
- *What each represents:*
  - The model has no within-period birth date. The pre-state is housing chosen at t−1, occupied until t.
  - The birth is realized at the start of t, before the date-t housing choice.
  - (a) is the first post-birth housing choice, occupied over [t, t+4), with no second births possible.
  - (b) is housing chosen at t+1, occupied over [t+4, t+8). It follows the treated branch's possible second birth, with the control held childless.
  - The data's −1 and +3 are one year before and three years after the birth. About 50% of mothers have a documented second child by +3 (50.37%; 60.72% by +4).
  - Mapping (a) or (b) onto event years requires an explicit dating convention. Neither can match exactly.
- *Disposition:* keep (b), the current advanced branch, as the scored counterpart and report (a) as a diagnostic. Change only if evidence supports it.
- *Pass:* (b) reproduces 1.341457 at B; treated and control masses are equal.

**3. Old-age level row.**
- *Input:* `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_old_wealth/wealth_results.csv`.
- *Operation:* add a named candidate row and an observer computing weighted p90/p50 of beginning-of-period wealth, with no income division.
- *Pass:* the empirical value reproduces 4.9194287748 (SE 0.5798039598); B returns 3.7277115, which equals its current ratio statistic.

**4. Entry common-scale construction.**
- *Definition:* for July-sample row i in wave t, $\omega_i=W_i/\bar E_t$. Here $\bar E_t$ is the weighted mean of reference person + spouse gross earnings among reference persons aged 18–65 in wave t. Build weighted quintile-bin means of ω within weighted earnings terciles, with zeros included.
- *Model mapping:* $b=\omega\,\bar E^{m}$, where $\bar E^{m}$ is model annual gross labor earnings per working household.
- *What it is not:* $E[\omega]$ is not the pooled ratio of sums, and discretizing the bins need not preserve the rank correlation. Report both, without thresholds.
- *Pass:*
  - The node table reproduces the sample's weighted mean of ω exactly.
  - All entry mass is feasible under the unchanged gates, with no censoring.
  - Any infeasible mass is reported as a failure. It is not censored, transferred or set to zero.

**5. Family-room candidate SE.**
- *Input:* `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_acs_child_mapping/recompute_child_mapping.py`.
- *Operation:* add the existing 1,000-draw metro bootstrap.
- *Pass:* the same run reproduces both 0.347067 and 0.336220.

**6. All-order late fertility (data acquisition).**
- *Input:* the archived NCHS 2003–06 natality public-use files, keeping the mother's-age (`mager`) and live-birth-order (`lbo_rec`) fields.
- *Operation:* tabulate births by maternal age and birth order, stating the unknown-order policy.
- *Pass:* first births reproduce 35,414 at 42+ and 6,611,269 in total.

## 4. Withdrawn or narrowed phase-3 claims

- **Withdrawn:** "keep 0.7202 as settled." Its similarity to the +2-minus-−2 contrast certifies nothing, because the normalization witness changes it without changing fitted observations.
- **Withdrawn:** "the omitted −2 is only a display normalization." Cohort support differs across reference periods, so switching references changes the estimand.
- **Withdrawn:** "birth-date-only housing should be primary." The data include later births, and a single housing decision is not a criterion for choosing the counterpart.
- **Withdrawn:** "the jump cancels within parents, so the family slope identifies δ_a alone," and "the extra coordinate makes lower loss partly mechanical."
- **Withdrawn:** "report censored entry mass under the frontier rule," and "the pooled node mean equals the ratio of sums."
- **Narrowed:** "stationarity requires 2.1." It holds only as the production contract's gate and entrant conversion. Replacing the 3+ value is an economic change.
- **Narrowed:** the late-fertility tail is established for the saved B point only, and without an all-order data comparison.
- **Narrowed:** PSID's omission of employer DC accounts is a source-coverage statement, not a demonstrated model–data mismatch.
- **Narrowed:** the family-room grouping issue is now quantified at −0.0108 rooms. It is an acknowledged approximation, not a bug.

No cluster outcomes are claimed; collection still awaits the login renewal.
