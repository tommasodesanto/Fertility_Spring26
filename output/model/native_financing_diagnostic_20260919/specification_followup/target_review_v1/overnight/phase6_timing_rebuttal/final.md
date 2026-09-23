# Rebuttal on the first-birth room observer (phase 5)

## 1. The lead's corrections

It helps to separate three objects:

- **What the code does (exact).** The observer copies identical pre-choice childless states at the origin date, t. One copy gets a first birth; the other stays childless. Both advance one model period. Only the treated copy may have a second birth at t+1. The observer then compares realized rooms at t+1, weighting by the stationary flow of first births. There is no calendar birth date $U$ and no regression.
- **A proposed observation convention (not in code).** Housing is constant within each four-year period. A birth decided at t falls at calendar time $t+U$, $U\in[0,4)$. Observed rooms equal the housing held at that calendar time.
- **An empirically identified estimand.** The PSID Sun–Abraham coefficients are cohort-specific effects on the treated, relative to −2 and to never-treated women. They are identified only under parallel trends and no anticipation relative to −2.

**Corrections accepted:**

1. **"−4/+4 is the proven unique match" is withdrawn on two grounds.**
   - The identity $E[\text{contrast}]=D_1$ holds only under the convention above. Here $D_0$ and $D_1$ are the treated-minus-control room differences at t and t+1. The convention also needs the model's built-in zero effect before the origin date (the copies share their history) and compatible cohort weights. None of this is code.
   - Even under those assumptions, uniqueness is false. Every window $[k_1,+4]$ with $k_1\le-4$ also maps to $D_1$ for any law of $U$, because the effect before the origin date is zero by construction.
2. **The "−1/+3 never equals $D_1$" claim was false.** With $a=\Pr(U<1)$, the contrast is $aD_0+(1-a)(D_1-D_0)$. It equals $D_1$ exactly when $aD_1=(2a-1)D_0$. For example, $a=1$ with $D_0=D_1$ satisfies this, as the lead showed. The correct statement is narrower: no value of $a$ makes the equality hold for all $(D_0,D_1)$, since that would require $a=0$ and $a=\tfrac12$ at once. Equality at particular effect values validates nothing.
3. **The scalar formulas also need independence.** The timing indicator must be independent of, or uncorrelated with, household-level effects. Otherwise conditional means replace $D_0$ and $D_1$. So $-2/+2=\tfrac12 D_1$ is not a rescaling rule.
4. **Two further gaps in my phase-5 algebra.**
   - Dating the model's origin-date birth at $t+U$ reinterprets a child who is present in period t's utility and budget. With $U>0$, the date-t housing choice is observed before the calendar birth. That is "anticipation" created by the convention, and it has to be declared.
   - The model's forced same-state control is not the regression's never-treated comparison. Never-mothers are selected on income, wealth and tastes, so model-world parallel trends generally fail.

**Other lead corrections, accepted without dispute:**
- The primary arms use 15 income states, not 7.
- Diagnostic parameter estimates and a qualified finite-search working preference should be reported; I withdraw my blanket ban.
- The old-age wealth-level and under-18 dependent alternatives are author choices, not automatic fixes.

## 2. Proposed protocol: run the PSID regression on model-simulated panels

This matches the same estimator on model-generated panels, in the style of indirect inference. It is neither exact nor standard by assertion. It changes no behavior, utility or income process.

- **A. Histories.** Using the frozen stationary policies, transition kernels and entry law, simulate household histories from age 18 to 85 at four-year model dates. Record age, children ever born, current dependents, tenure, uncapped rooms, income and survival. All subsequent births are kept.
- **B. Declared timing convention.**
  - Baseline: $U\sim\text{Uniform}[0,4)$, independent of the household state, with housing held constant between model dates.
  - Report two bounding conventions alongside: $U\equiv0$ and $U\equiv2$.
- **C. Match the PSID observation pattern.**
  - Assign simulated treated households to PSID first-birth cohorts and interview-year patterns by resampling the fitted cohort × event-time support, with its IW mass. That support table is already exported by the replay (`estimation_cohort_event_support.csv`).
  - Controls are simulated households still childless through their fertile ages. They get PSID control observation patterns.
  - This reproduces the annual/biennial composition of the data.
- **D. Same estimator.**
  - Run the identical cohort × event interaction regression with never-treated controls: person and year fixed effects, the same event dummies and omitted periods, and age controls.
  - The model has no education, so that control is dropped; disclose this.
  - Compute exactly the same aggregated or common-cohort contrast used for the data.
- **E. Separating pre-birth dynamics from the intervention effect (model world only).**
  - The simulated pre-period coefficients (−4, −3, −1 relative to −2) measure the model's own selection plus convention-induced anticipation.
  - Also report the intervention contrasts: $D_1$ from the existing observer, and $D_0$ from the same copied branches measured at the birth date (a small measurement addition).
  - The gap between the simulated regression contrast and $D_1$ bounds how much the regression mixes selection, anticipation, weighting and convention into the effect.
- **Output.** The model moment is the simulated regression contrast, with a Monte Carlo SE across seeds. The pre-period coefficients and $D_0$/$D_1$ are diagnostics.

**Feasibility.** Existing saved output cannot support this: it holds only aggregate branch means and the continuation-birth mass. New measurement-only code is required, either a history simulator or a forward distribution tagged with event time. No re-solve is needed if the frozen checkpoints' policies are retrievable; they are currently on /scratch, behind the expired SSH login.

## 3. Recommended morning decision

- **Frozen.** Tonight's comparison keeps 0.7202 and the native observer. The room row is reported as a diagnostic whose estimand is unresolved. No rescaling and no v2 contract.
- **Decide next: the measurement principle.** The model moment for the room response should be the same PSID contrast computed on simulated panels under the declared convention. It should not be the intervention contrast compared directly with the regression.
- **Window: none is justified yet.** Two facts are needed first:
  1. The within-cohort pre-birth contrasts in the data (−4→−2 and −2→−1), with SEs. These are being computed mechanically elsewhere.
  2. How sensitive the simulated contrast is to $U$ (uniform versus 0 versus 2).
- **Conditional rule.**
  - If pre-birth changes are small relative to their SEs and to the post-birth change, use the common-cohort −2→+2 contrast. It respects the author's −2 reference and has a conditional SE.
  - If pre-birth changes are material, start the window before the rise (e.g., −4), and require the simulated panel to reproduce the pre-period coefficients.
- **What would change this.** If the simulated contrast moves materially across the $U$ conventions, the four-year model cannot discipline annual event-time contrasts. I would then recommend a timing-robust statistic instead: mean rooms over event years +1 to +4 minus mean rooms over −4 to −1, computed identically on both sides. That also needs author approval.
