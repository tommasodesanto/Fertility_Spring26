# PSID Empirical Tests for the Spatial-Fertility Theory Block

## Project context

This is a structural spatial-lifecycle model with endogenous fertility,
housing tenure, and location choice. Two locations (Center $C$ and
Periphery $P$) with $p_C > p_P$ (rents) and $\mathcal E_C > \mathcal E_P$
(amenity premium). Common wages across locations. Households have
Stone–Geary preferences nested in a CRRA wrapper:

$$u(c,h;m) = \frac{\bigl[(c-\bar c(m))^\alpha (h-\bar h(m))^{1-\alpha}\bigr]^{1-\sigma}}{1-\sigma}$$

with calibration $\sigma=2$, $\alpha=0.70$, and parity-dependent
subsistence: $\bar h(0) < \bar h(1)$. The discrete jump $z_1 = \bar h(1) - \bar h(0)$
is the key first-birth housing-requirement shifter.

The theory block delivers four sharp predictions for the data. Three
are **untested** in the existing PSID work and need new estimation.
This document briefs you on what to run and where the existing
infrastructure lives.

---

## Repo orientation

**Existing PSID empirical infrastructure** (do not duplicate; extend):

- `code/data/psid_followup_mar2026/fertility_wealth_v1.do` through
  `v12` — the wealth-fertility analysis pipeline. v1 produces the
  "within-wealth-gap" chart (32pp ownership-transition gap by quintile).
  v12 is the simplest variant.
- `code/data/psid_followup_mar2026/sa_*.do` — Sun–Abraham staggered-DiD
  event studies for first and second births.
- `code/data/psid_followup_mar2026/moved_for_size_iv_alltenure_v1.do`
  — IV for moving for housing space, using twin-firstbirth and
  same-sex-first-two-kids as instruments.
- `code/data/psid_followup_mar2026/output/` — output subdirectories for
  each analysis. Mirror the naming convention.

**Required external data merges**:

- HUD Fair Market Rents (FMR) by year and MSA, public:
  https://www.huduser.gov/portal/datasets/fmr.html
- Or Census ACS median gross rent by MSA, B25064.
- PSID has MSA codes (variable typically `er<year>5004` or similar
  state/county/MSA identifiers). Verify variable names in the existing
  do-files before merging.

**Output naming convention**: each test gets its own subdirectory
`output/<test_name>_v1/` with a `.log`, a `.csv` summary table, a
`.tex` regression table, and a PDF figure when applicable. Match the
style of `rooms_first_birth_one_vs_two_horizon_v1/`.

---

## What's already tested (do not redo)

| Prediction | Status | Existing output |
|---|---|---|
| Cross-sectional sorting (parents less central) | ACS | mentioned in main slide deck; see `latex/april_20_project_presentation.tex` |
| First-birth raises rooms / sqft / moves-for-space | PSID | `rooms_first_birth_one_vs_two_horizon_v1`, `sa_rooms_first_birth_*`, `moved_for_size_iv_alltenure_v1` |
| First-vs-second-birth housing-step heterogeneity | PSID | `rooms_first_birth_one_vs_two_horizon_v1` (1 kid: +0.26 rooms; 2+: +1.37 rooms) |
| Within-wealth ownership-fertility gap | PSID | `fertility_wealth_v1` (32pp gap, flat across wealth quintiles) |
| Hump-shaped wealth-fertility | PSID | `lpoly_parenthood_v12` |

---

## Tests to implement (priority order)

### TEST 1 (HIGH PRIORITY) — Directional flight + wealth gradient

**Theory.** Combines two predictions:
- **P1b**: First birth triggers relocation toward lower-rent areas.
- **P3a**: Flight rate is decreasing in pre-birth wealth (because the
  sorting threshold $\bar h^\ast(M)$ is increasing in $M$, so high-wealth
  parents can absorb $\bar h(1)$ in $C$).

**Design.**

1. **Merge HUD FMR (or ACS B25064 median rent)** with PSID by year and
   MSA. Construct `rent_msa_{i,t}` for each household-year.

2. **Define directional flight indicator**:
   ```
   flight_{i,t} = 1{moved_{i,t} == 1 AND rent_msa_{i,t} < rent_msa_{i,t-1}}
   ```
   Two flavors: strict (rent strictly lower) and quantile (moved to
   bottom-tercile rent MSA from a higher tercile). Report both.

3. **Sun–Abraham event study around first birth** (mirror
   `sa_rooms_first_birth_variants_v1`):
   ```
   flight_{i,t} = sum_k beta_k · 1{event_time_{i,t} = k}
                  + alpha_i + gamma_t + X_{i,t}'phi + e_{i,t}
   ```
   k ranges over event-time bins (e.g., −5 to +5 years), with k=−2 as
   reference. Sample: childless household-years prior to first observed
   birth, plus post-birth follow-up window.

4. **Interaction with pre-birth wealth quintile**. Define wealth
   quintile from net worth or income at age 25–30 (mirror
   `fertility_wealth_v1` definition). Interact event-time bins with
   quintile dummies, or run the event study separately within each
   quintile.

5. **Outputs**:
   - Pooled event-study coefficient plot (figure).
   - 5-panel figure: event-study by wealth quintile.
   - Summary table: post-birth average flight rate by quintile, with
     standard errors. Predicted: monotone decreasing in wealth quintile.
   - One-line headline number: "X% of low-wealth parents moved to a
     lower-rent MSA within 2 years of first birth, vs Y% of high-wealth
     parents."

**Why this matters**: this is the headline test for the spatial-sorting
mechanism. Together, P1b and P3a put a sharp wedge between the
housing-flight story and any pure-amenity sorting story.

---

### TEST 2 (HIGH PRIORITY) — Housing-elasticity attenuation at parenthood

**Theory.**
- **P2a**: Income elasticity of housing demand is smaller for parents
  than childless. Stone–Geary closed-form: $\eta_M(n) = (1-\alpha)(M-\bar c)/(p\,h^\ast(n))$,
  decreasing in $\bar h(n)$.
- **P2b**: Same for price elasticity, $\eta_p(n) = -\eta_M(n)$.

**Design.**

1. **Outcome**: housing services. Best measure available is rent paid
   (renters) plus owner-equivalent rent imputed from house value
   (owners) — same measure used in your existing rooms work where
   applicable. Secondary: rooms.

2. **Income-elasticity regression** (P2a):
   ```
   log(housing_{i,t}) = beta_y · log(income_{i,t})
                       + beta_int · I(any_kid_{i,t}) · log(income_{i,t})
                       + delta · I(any_kid_{i,t})
                       + alpha_i + gamma_t + X_{i,t}'phi + e_{i,t}
   ```
   Household FE absorbs time-invariant heterogeneity; identification
   from within-household income variation. **Predicted**: $\beta_y > 0$,
   $\beta_{\text{int}} < 0$, with $\beta_y + \beta_{\text{int}}$ smaller
   than $\beta_y$ for parents.

3. **Price-elasticity regression** (P2b): replace `log(income)` with
   `log(rent_msa)` (from the FMR/ACS merge in Test 1). Identification
   from movers across MSAs and from MSA-level rent shocks. Predicted:
   $\beta_p < 0$ and parents' coefficient closer to zero in magnitude.

4. **Robustness**: alternative kid definitions (any kid under 18, any
   kid under 5, kid count). Run also as continuous-in-$\bar h$ proxy
   using kid-count.

5. **Outputs**:
   - Two regression tables (income elasticity, price elasticity) —
     side-by-side parents vs childless coefficients.
   - One coefficient comparison figure (bar chart of $\hat\beta_y$ by
     parental status, with 95% CI).
   - Headline number: "Parents' housing income elasticity is X
     [precision]; childless is Y [precision]; difference statistically
     significant at p<0.05."
   - Implied $z_1 = \bar h(1) - \bar h(0)$ from the elasticity gap, given
     $\alpha = 0.70$ — for cross-checking with the calibrated structural
     value (currently $z_1 \approx 0.66$ from the rooms event study).

---

### TEST 3 (MEDIUM PRIORITY) — Slack-dependent first-birth hazard

**Theory.**
- **L3b** (Lemma 3): Constrained households on a fixed housing branch
  face a *power-law* divergent cost of accommodating the first child as
  pre-birth slack $s^0 = h - \bar h(0)$ approaches the parent-subsistence
  increment $z_1$. Closed-form CEV:
  $$\mathrm{CEV} = (c-\bar c)\bigl[(s^0/(s^0-z_1))^{(1-\alpha)/\alpha} - 1\bigr]$$
  Interior renters and owners pay flat (linear) cost; capped renters
  pay the super-linear cost.

**Design.**

1. **Construct pre-first-birth slack**:
   - Estimate $\widehat{\bar h(0)}$ from childless renters' average rooms
     (or use a quantile, e.g., 10th percentile) — the "subsistence
     floor" for childless. From the rooms event-study work,
     pre-first-birth mean is $\sim 6$ rooms, so $\bar h(0)$ is somewhere
     around 4–5 rooms.
   - Slack: `slack_{i,t} = rooms_{i,t} - bar_h_0_hat`.
   - Construct in the year prior to first birth (or moving 12-month
     window).

2. **Hazard model**:
   ```
   Pr(first_birth_{i,t+1} = 1 | childless at t)
     = Phi(g(slack_{i,t}) + h(income_{i,t}) + X_{i,t}'phi + alpha_t)
   ```
   With `g` a flexible function — spline or 5-bin step function. Run
   parametric (probit) and nonparametric (binned local-linear).

3. **Predicted shape**: hazard strongly decreasing in slack at low
   slack (super-linear cost), flattening at high slack (interior
   household, constraint slack).

4. **Robustness — IV for slack**: pre-birth slack is endogenous
   (people who plan to have kids choose larger units). Candidate
   instruments:
   - MSA housing supply elasticity (Saiz 2010) interacted with national
     rent inflation: shifts available unit sizes exogenously.
   - Local rent control / vacancy decontrol shocks.
   - Inherited / family-provided housing (if PSID has).
   First-stage: predict slack from instrument(s); second-stage hazard
   uses predicted slack.

5. **Tenure cut**: run separately for renters and owners. Predicted:
   the super-linear shape is concentrated in renters (they have the
   binding cap $\bar h_R$); owners with multiple rungs face flatter
   cost.

6. **Outputs**:
   - Hazard plot: predicted first-birth probability vs slack
     (nonparametric), with confidence band.
   - Regression table: slack coefficients (parametric).
   - Comparison renters vs owners — should look qualitatively
     different (steeper for renters).
   - Headline: "Renters in the bottom slack tercile have a $X$pp
     lower first-birth hazard than renters in the top slack tercile,
     after conditioning on income and household FE."

---

### TEST 4 (LOWER PRIORITY but cheap) — Pre-vs-post rent-share jump by wealth

**Theory.** Stone–Geary closed-form for the housing-expenditure share
jump at first birth (under (B2) $\bar c(0) = \bar c(1)$):
$$\Delta(p h^\ast / M) = \alpha\, p\, z_1 / M$$
**Declining in $M$** — low-wealth parents take a bigger budget hit
from family housing.

**Design.** Compute pre-birth and post-birth housing-expenditure share
in PSID. Run event study on `housing_share_{i,t}` interacted with
pre-birth wealth quintile.

**Predicted**: post-birth jump in housing share is largest for
low-wealth quintile, smallest for high-wealth quintile.

This is a *very* cheap test using existing PSID variables and the
existing event-study infrastructure. One additional outcome variable;
one additional figure.

---

## Sun–Abraham specification (use throughout for event studies)

Mirror the existing `sa_rooms_first_birth_variants_v1` spec. Key points:
- Treatment = year of first birth (or first observed birth in the
  panel, with appropriate sample restriction).
- Cohort by year-of-treatment.
- Control = never-treated by end of panel (or last-treated, depending
  on variant — report both).
- Bins: $k \in \{-5, \ldots, +5\}$ with $k = -2$ as reference.
- Cluster SEs at household level.

Run with the `eventstudyinteract` Stata package if available; else
manual SA via dummy interactions.

---

## Coding standards

- Stata `.do` files mirror the existing naming convention
  (`<analysis_name>_v1.do`).
- Each script writes a `.log`, a `.csv` summary, and saves figures to
  `output/<analysis_name>_v1/`.
- Headline numbers go into a one-line summary CSV that I can paste
  into the slide deck without reformatting.
- If the test design needs a structural decision (e.g., how to
  construct the wealth quintile, how to define "first birth", how to
  handle missing MSA codes), follow the convention in
  `fertility_wealth_v1.do` unless deviating gives a clearly better
  identification, in which case flag the deviation in the log.

---

## Reporting back

For each test, produce:

1. A short summary (3–5 sentences) of the design choices made.
2. The headline result (one number with one CI).
3. The figure file path.
4. The regression table file path.
5. A flag if the result *contradicts* the model's prediction — these
   are useful too, and need to be reported up rather than buried.

Run tests in priority order (1 → 2 → 3 → 4). Stop and report after
each one before proceeding to the next, so calibration adjustments
can be made if needed.