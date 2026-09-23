# Decision memo: integrating the room, entry, family-room and old-age evidence

**Scope.** This pass was read-only: no model runs, no data changes, and nothing adopted. The 40 submitted cluster jobs keep their frozen targets and parameters. Everything below is a candidate measurement choice for the author.

## 1. What the new evidence changes

- **First-birth room response (the PSID event-study row).**
  - An exact replay of the original regression reproduces the current target: 0.7202462624, SE 0.0852600513, on 49,457 observations of 4,112 women. It also exports the full 900×900 covariance matrix of the cohort-by-event coefficients.
  - That covariance gives standard errors for contrasts computed within a fixed set of common birth-year cohorts. These SEs hold the chosen cohort weights fixed; they exclude the uncertainty from estimating those weights.
  - The earlier normalization witness still stands: the pooled contrast depends on additive constants for the specific cohorts the witness identified. Taking differences within each cohort removes that dependence.
- **Room codes in the PSID are not a bug.** The 100 answers coded 0 mean "shared room / no exclusive room" and are substantive. The single code 98 belongs to a never-treated control woman in the 1991 wave, where 98 is a permitted count. No recode is justified.
- **Family-room slope.** Counting only resident own children under 18 changes the gap by −0.0108, with a paired SE of 0.0081. The case for the minor-only definition is conceptual, not the number: it is closer to the model's idea of a dependent child.
- **Entry wealth (common-scale candidate).**
  - The candidate's mean, 0.1858, is essentially the inherited level marginal's 0.1865. What changes is the shape of the distribution and its association with earnings, not the level.
  - It passes the current-period resource check with strictly positive slack for all entrant mass, using the actual frozen parameters.
- **Old-age wealth dispersion.** Nothing new beyond the earlier finding: the model's constant pension makes its ratio and level statistics identical (3.7277).

## 2. Room response: target and observer (Q1)

### The timing algebra

This is algebra conditional on three facts about the model, not a behavioral claim:
1. In the model, housing changes only at four-year dates: $h_{t-1}$ is held over $[t-4,t)$, $h_t$ over $[t,t+4)$ and $h_{t+1}$ over $[t+4,t+8)$.
2. The first birth is realized at date t, before the date-t housing choice. Its calendar date $t+u$, with $u\in[0,4)$, is not modeled.
3. The data's rooms at event time $k$ are the housing held at calendar time $t+u+k$.

Define $\Delta_s=E[h^T_s-h^C_s]$: the treated-minus-control difference in housing at date $s$, starting from identical pre-states. The current observer computes $\Delta_{t+1}$.

| Data window | Model-implied counterpart |
|---|---|
| $-4\to+4$ | $\Delta_{t+1}$ exactly, for **any** distribution of $u$ |
| $-1\to+3$ (current target) | $\Pr(u<1)\,\Delta_t+\Pr(u\ge1)\,(\Delta_{t+1}-\Delta_t)$. This equals $\Delta_{t+1}$ under no dating convention unless $\Delta_t=0$ |
| $-2\to+2$ | $\Pr(u<2)\,\Delta_t+\Pr(u\ge2)(\Delta_{t+1}-\Delta_t)$. Under uniform dating within cells, which the CPS observer already assumes, this is $\tfrac12\Delta_{t+1}$ |
| Symmetric $[-a,a]$, uniform dating | $(a/4)\,\Delta_{t+1}$ |

### What this implies

**I challenge the lead's framing here.** The current observer does not correspond to a four-year data window. Its exact data counterpart is the eight-year window $-4\to+4$. The frozen pairing of the pooled −1→+3 target with $\Delta_{t+1}$ is mismatched under every possible within-period birth date.

Choosing −1→+3 because it "is four years" is therefore not justified. A matched two-year window (−1→+1) maps to $\tfrac14\Delta_{t+1}$ under uniform dating: it depends on the convention and has no advantage.

### Recommendation

- **Model side:** keep the current observer unchanged. No new model implementation is needed.
- **Data side:** replace the target with the common-cohort within-cohort contrast
  $\sum_g q_g[\beta_g(+4)-\beta_g(-4)]$.
  - Here $\beta_g(k)$ is cohort g's event-time coefficient, and $q_g$ are fixed pre-birth weights (the IW weights at event −4).
  - Report post-birth (+4) weights as a sensitivity check.
  - −2 stays the estimator's omitted reference period, so the author's reference is preserved inside the estimation.
- **Why this window:**
  - It is the only data window whose model counterpart does not depend on the unobservable birth date within the period.
  - Taking differences within each cohort removes the additive cohort constants behind the witness.
  - The saved $e(V)$ supplies its conditional SE.
- **Why pre-birth weights:** the model weights by the stationary flow of first births from pre-birth states, which the pre-birth endpoint approximates. Weights must not be chosen by fit.

### Who the two sides describe

- **Data:** women who are reference persons or spouses, in single-family-unit dwellings, with a first birth in cohorts observed at both −4 and +4. That means annual-era cohorts plus odd-year biennial cohorts. Controls are confirmed never-treated women. Rooms are measured four years before and four years after the birth. About 61% have a second child by +4 (60.72%; 59.52% counting only strictly later birth years).
- **Model:** households with a first birth realized at date t from childless, settled states at ages 18–42, weighted by the stationary first-birth flow. They are compared at t+1, after that date's fertility and housing decisions, against a same-state copy held childless. Treated households then have one or two current dependents (two after a second birth at t+1); controls have none.

### Subsequent births

Second births sit in both objects by construction. Removing them would condition on an outcome, so neither side should. As a diagnostic only, report the model's share of treated households with a second birth at t+1 next to the PSID share by +4. They are not the same timing object: under uniform dating, a second birth at date t+1 falls between event times 0 and 8.

### What still blocks a final choice

1. The −4/+4 cohort support, the estimate and its SE have not been computed. They come from the saved coefficients and covariance; no new regression is needed.
2. The longer window is more exposed to housing changes before the birth that the model does not represent, such as couples moving in together. The within-cohort change from −4 to −2 measures that exposure and uses the same files.
   - If that change is material, the author chooses between convention-free −4→+4 and −2→+2 paired with $\tfrac12\Delta_{t+1}$ (which depends on the uniform-dating convention).
   - The choice must be made before looking at the arms' fits.

**Note on the existing candidates.** The already-computed −1→+3 (0.811/0.798) and −2→+2 (0.735/0.704) contrasts use different biennial cohort sets: even-year births for the first, odd-year births for the second. Differences between them mix window and composition. Neither should be paired with $\Delta_{t+1}$ in final inference. The frozen 0.7202 stays only in the running experiment.

## 3. Entry wealth (Q2)

The candidate defines each entrant's wealth as $\omega$ = nonhousing net worth ÷ the same-wave mean gross earnings of working reference-person households. Taking its pieces separately:

| Piece | Assessment |
|---|---|
| **Denominator** | **Preferable.** It uses the same concept as the aggregate wealth/earnings target. The model's mean annual working-age gross earnings is exactly 1, so model entry wealth equals $\omega$. It never divides by the entrant's own earnings, and it no longer depends on the retired income process that generated the inherited wealth levels. |
| **Wealth–earnings association** | **Better than inherited.** The imposed rank correlation is 0.241, against 0.2825 in the data and 0.098 in the inherited coupling. The within-tercile coupling puts 0.28% of mass on nine income-state/wealth combinations that the finer microdata ranking never produces; disclose it. |
| **Discretization (three earnings terciles × five wealth bins)** | **Main open issue.** The 15 cell means leave 79.1% of the variance of $\omega$ inside cells. The nodes span [−0.547, 1.567], while the inherited reference points spanned about [−8.2, 11.5]. The down-payment channel is a threshold, so tail mass matters. A compression that preserves the mean can move the eligible share either way. |
| **Age mismatch** | Unchanged July assumption: an 18–24 sample stands in for entrants at 18. The common scale does not fix it. |
| **Sample and vintage** | The earnings scale is extended back to the 1984–2003 waves with the aggregate builder's rules plus an alive-in-wave filter. It differs by at most 0.024% in overlapping waves. Document it as a support extension. |
| **Numerical validity** | The current-period renter resource check passes for all mass: minimum slack 0.082 at the nodes, 0.072 on interpolation support. There is no Bellman, continuation or equilibrium validation, and no censoring. |

**Verdict.** The candidate is economically preferable in units and in its association with earnings. It is not yet adequate as the paper's baseline until its tail representation is checked against the purchase threshold.

**Smallest next steps, after the author agrees:**
1. **No-solve check.** Using B's saved price and the model's eligibility rule ($b+Y/R\ge(1-\phi)pH$ for each owner size at age 18), compare eligible shares under the row-level $\omega$ distribution (mapped through earnings-rank overlap), the 3×5 nodes and the inherited marginal.
   - If the row-level and discretized shares differ in a way the author judges material, split the lowest and highest $\omega$ quintiles within each tercile and repeat.
2. **Only then:** one stationary solve at B's parameters under the candidate law, reporting all 13 rows and gates. Label it an input-law sensitivity.

## 4. Corrections versus choices, and a plan for tomorrow (Q3)

**Measurement corrections (same population, mapping shown algebraically or in code):**
- **Room response:** the frozen observer/target pairing is algebraically mismatched. The −4→+4 pairing is the candidate replacement, pending the support and pre-trend computation.
- **Old-age dispersion:** the model's operator is a level dispersion, so the matched statistic is the same-sample level ratio, 4.9194 (SE 0.5798).
  - B's fit is worse under it (2.05 SE), which argues against any suspicion that it was chosen for fit.
- **Family-room slope:** switch to linked children under 18, 0.3362 (SE 0.0584). This is a small definitional alignment. It does not validate the model's dependent-exit process.

**Not corrections:** the room codes; aggregate wealth; the fertility timing and CPS stock rows; ownership; mean rooms; the recent-parent gap.

**Substantive author choices:**
- National versus 42-metro geography.
- Top-group vintage (3.4828 versus 3.6024), with 2.1 retained as the production normalization.
- Adopting the entry law.
- Room-contrast weights, and the window fallback if the pre-birth change is material.
- Disclosure wording for the borrowed bequest restriction.

**Plan (safe for the frozen run; at most five tasks):**
1. **Room windows from the existing replay (no regression).**
   - *Input:* `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_rooms/covariance_replay/`.
   - *Operation:* compute −4→+4 with pre- and post-birth weights, plus the within-cohort −4→−2 change, with conditional SEs.
   - *Done when:* the four existing contrasts reproduce to 1e-12, and cohort lists, retained mass, estimates and SEs are reported.
2. **Candidate contract v2, separate from the frozen contract.**
   - *Content:* the old-age level row, the minor-only family row, and the −4→+4 room row with the unchanged observer; every other row as now; source hash and SE type recorded per row.
   - *Done when:* re-scoring B's saved moments (no solve) reproduces its model values (old-age 3.727712 through a level observer), and the frozen fingerprint is untouched.
3. **Entry threshold check** as described in §3.
   - *Done when:* the eligibility and tail table is complete for all three distributions; any tail refinement follows the author's judgment, with no preset cutoff.
4. **When collection resumes.**
   - *Operation:* for each arm, extract the saved $\Delta_{t+1}$ and the treated second-birth share at t+1 (no solve). Report both arms under the frozen contract and under v2, as re-scorings of identical model moments.
   - *Done when:* both tables exist for every collected arm, with fingerprints.
5. **Author decision sheet.**
   - *Content:* geography, top-group vintage, entry law, room weights and window fallback, each with alternative values already computed and their uncertainty type.
   - *Done when:* every entry cites a hashed source; nothing is adopted.

## 5. What the overnight comparison can and cannot show (Q4)

**What it can show.** It is a useful controlled diagnostic: both utility classes face identical frozen targets, weights and observers. Legitimate conclusions:
- how each class moves each of the 12 scored rows at its searched point;
- whether their housing-by-children patterns differ qualitatively;
- feasibility and the standard diagnostic graphs.

**What it cannot support:** final claimed SMM inference, or a preference between the two utility classes. The reasons:
1. The weights are heterogeneous working scales, not efficient or covariance-based, so no J-test is possible.
2. Two scored rows have demonstrated mapping mismatches: the room pairing and the old-age denominator. These affect both arms, possibly unequally, so a loss difference may partly reflect mismeasured rows.
3. Population and vintage choices are still open.
4. The two utility classes are non-nested, the adaptive search paths diverge, and there are no Jacobians. Moment counts are not identification.
5. Numerical adequacy is unproven: the seven-state income grid and the wealth grid have not been shown to converge.

**Reporting rule.** Present the arms row by row under the frozen contract, and separately as a re-scoring under the candidate contract once it exists. Draw no parameter estimates, standard errors or policy conclusions from either.
