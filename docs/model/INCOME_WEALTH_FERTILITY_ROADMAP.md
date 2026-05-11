# Income, Wealth, Fertility, And Non-Homotheticity Roadmap

Date: `2026-04-25`

Purpose: organize the next research block around the income and wealth
distribution of fertility. The question is not only how to model
non-homothetic housing demand, but what joint facts the model must reproduce:
fertility by income, fertility by liquid wealth, fertility by housing wealth,
and the way these gradients differ by tenure, location, and age.

This note is a living roadmap. It should be updated when new PSID or model
diagnostics are added.

## 1. Core Research Question

The current model uses Stone-Geary housing and consumption needs,
\[
u(c,h;m)=
\frac{\left[(c-\bar c(m))^\alpha (h-\bar h(m))^{1-\alpha}\right]^{1-\sigma}}
{1-\sigma}
+\psi_{\text{child}}m,
\]
with child-dependent needs
\[
\bar c(m)=\bar c_0+\bar c_1m,\qquad
\bar h(m)=\bar h_0+\bar h_{\mathrm{jump}}\mathbf 1\{m\geq 1\}+\bar h_1m.
\]

This architecture mechanically makes family housing requirements more binding
for households close to subsistence. It therefore predicts a wealth gradient in
fertility even if all households have the same income profile and the same child
taste. That is the right object to study, but it creates a discipline problem:
the same non-homotheticity that repairs housing-by-parity can also generate the
wrong fertility-by-wealth or fertility-by-income gradient.

The work should answer three questions:

1. What do the data say about fertility over the income, liquid-wealth, total
   wealth, and housing-wealth distributions?
2. What does the current model imply for those same objects?
3. Where does the model fail, and can the failure be reconciled by changing
   preferences, constraints, heterogeneity, or measurement?

## 2. Literature Map

### A. Canonical Fertility Modeling

Becker (1960) and Becker--Lewis (1973) make children durable household goods
and emphasize quantity-quality tradeoffs. Their lesson for this project is that
``children are normal goods'' is not enough: the shadow price of children changes
with child quality, parental time, goods costs, and housing needs.

Jones, Schoonbroodt, and Tertilt (2008, 2011) show that standard theories do
not generically produce a negative fertility-income relationship without special
assumptions. That is a direct warning against forcing a monotone negative
gradient into the model without checking the modern U.S. data.

Sommer (2016) models fertility choice in a lifecycle setting with earnings and
infertility risk. Children behave like commitments, so higher risk can delay and
reduce fertility. This is relevant because our model currently has deterministic
age-income profiles but no idiosyncratic income risk.

Doepke and Kindermann (2019) show that fertility preferences and bargaining
within couples matter quantitatively. This is not in the current model, but it is
a candidate residual if wealth and housing constraints cannot explain
childlessness at the top of the distribution.

Bar, Hazan, Leukhina, Weiss, and Zoabi (2018) show that the U.S.
income-fertility relationship flattened between 1980 and 2010 as high-income
families increased fertility, and they rationalize this with marketization of
childcare. This is especially relevant because high-income households may buy
out of time costs even while housing costs bind lower-wealth households.

Doepke, Hannusch, Kindermann, and Tertilt (2022) summarize the new fertility
era: old facts such as a universally negative income-fertility gradient no
longer hold cleanly in high-income countries.

### B. Housing, Wealth, And Fertility

Lovenheim and Mumford (2013) use PSID housing wealth changes and find that
housing wealth shocks raise fertility among homeowners. Dettling and Kearney
(2014) distinguish the price effect from the home-equity effect: rising house
prices reduce births among non-owners but increase births among owners.

These papers imply that the model must not report a single housing-cost effect.
It must split:

- renters / future first-time buyers: higher prices are a child-cost shock;
- owners: higher prices are also a wealth shock;
- constrained renters: liquid wealth and down-payment constraints should matter
  most around first birth.

Kulu (2013) finds that housing conditions explain a substantial part of the
urban-rural first-birth gap but less of higher-parity variation. This lines up
with the model's first-child housing jump, but it also cautions that the second
and third child margins may need different objects than just higher
\(\bar h(m)\).

### C. Housing Demand And Non-Homotheticity

Albouy, Ehrlich, and Liu (2016) estimate housing demand with non-homothetic
features and emphasize that housing demand is income and price inelastic, with
large distributional cost-of-living implications for the poor. This matters for
our Stone-Geary block because a strong subsistence floor can over-amplify
housing sensitivity at the bottom and erase enough room for high-income fertility
mechanisms.

Davis and Ortalo-Magne (2011) document roughly stable housing expenditure
shares across cities and time. That is useful as a counterweight: if the
Stone-Geary housing block creates too much variation in housing shares by income
or location, it will conflict with this evidence.

### D. Spatial Equilibrium

Rosen (1979) and Roback (1982) remain the sorting backbone: rents, wages, and
amenities jointly rationalize location choice. In this model wages are common
across center and periphery, so fertility gradients should not be interpreted as
cross-labor-market income effects.

Baum-Snow and Han (2024) and Saiz (2010) are still the relevant supply anchors.
The supply side should not be allowed to absorb fertility-income failures unless
there is a direct housing-market reason.

## 3. What The Repository Already Tells Us

### A. PSID Wealth-Fertility Evidence

Existing infrastructure:

- `code/data/psid_followup_mar2026/fertility_wealth_v1.do` through `v12`
- `code/data/psid_followup_mar2026/output/fertility_wealth_v1/`
- `docs/prompts/PSID_TEST_PROMPT.md`

Current usable fact from `v12`:

- parenthood by age 35 over young-adult net-worth quintiles is not monotone:
  `0.703`, `0.741`, `0.739`, `0.761`, `0.642`;
- parenthood by age 45 is flatter but still lower in the top quintile:
  `0.729`, `0.766`, `0.775`, `0.796`, `0.683`.

Interpretation:

- this is not evidence for a simple monotone positive wealth-fertility gradient;
- the middle and upper-middle of the wealth distribution look most fertile;
- the top quintile may reflect delayed fertility, preferences/career costs,
  sample selection, or measurement of wealth at ages `25-30`;
- some older output files in the same folder show all-zero summaries and should
  not be treated as evidence without script-level validation.

Existing PSID note also claims:

- a within-wealth ownership-fertility gap of about `32pp`;
- a hump-shaped wealth-fertility relation;
- birth-related moves and housing-size increases.

Those claims are plausible and model-relevant, but they need one consolidated,
audited output table before becoming calibration or validation targets.

### B. PSID Housing Event Studies

Live targets already use:

- first birth `+3` rooms: `0.664`;
- second birth `+3` rooms, stay-at-one controls: `0.566`;
- second-birth last-treated robustness: `0.704`.

These identify birth-related housing responses, not completed-fertility
parity-transition probabilities. The current one-shot fertility architecture
therefore cannot turn them into a sequential second-birth hazard target.

### C. ACS/MMS Spatial Facts

The current live MMS target system implies:

- center population share target: `0.450`;
- center/periphery unit-rent ratio: `1.140`;
- fertility gradient: `0.133`;
- ownership rate: `0.627`;
- ownership gradient: `0.170`;
- new-parent ownership gap: `0.110`;
- center share nonparents: `0.494`;
- center share new parents: `0.416`.

These are spatial and tenure facts, not income-gradient facts. The missing
empirical object is a common sample that cross-tabs fertility by:

- income or permanent income;
- liquid wealth;
- total wealth;
- housing wealth / ownership;
- center/periphery status;
- age and completed-vs-flow fertility concept.

## 4. What The Current Model Tells Us

Current benchmark from `CALIBRATION_STATUS.md`:

- strict loss: `21.33`;
- TFR: `1.969` vs target `1.700`;
- childlessness: `0.123` vs target `0.150`;
- mean age first birth: `34.26` vs target `26.00`;
- fertility gradient: `0.205` vs target `0.133`;
- ownership: `0.902` vs target `0.627`;
- ownership gradient: `0.060` vs target `0.170`;
- young liquid wealth/income: `0.461` vs target `0.600`;
- old-age ownership: `0.974` vs target `0.863`;
- old-age parent-childless ownership gap: `0.026` vs target `0.070`.

Direct implications:

- the live benchmark is too fertile, too old at first birth, and over-owns;
- the model's fertility gradient is currently spatial, not an income gradient;
- it has no persistent income heterogeneity, only a deterministic age-income
  profile and wealth heterogeneity from lifecycle saving, tenure, location, and
  fertility choices;
- therefore it cannot yet speak cleanly to the income distribution of fertility
  except through age and endogenous wealth sorting.

Existing diagnostic plot:

- `output/model/figures_current_candidate/diag09_fertility_by_wealth.png`

Read from the current policy plot:

- for renters, the probability of staying childless generally falls with liquid
  wealth at older fertile ages;
- the probability of choosing one child is hump-shaped in liquid wealth for many
  ages;
- this is a policy-function object, not a population moment, and must be
  converted into distribution-weighted wealth-bin moments before comparison to
  PSID.

Important limitation:

- a model with common wages and no idiosyncratic earnings states can match
  fertility by wealth only through asset accumulation, constraints, and tenure;
- if the data show a distinct fertility gradient by income conditional on wealth
  and tenure, the current state space is missing a major dimension.

## 5. Likely Failure Modes

### A. Wrong Monotonicity

The current Stone-Geary schedule may predict that fertility rises monotonically
with wealth among renters because richer households can clear child housing
needs. The PSID evidence points toward a hump shape: low wealth constrains
fertility, middle wealth supports parenthood, and top wealth may delay or reduce
fertility.

Candidate reconciliation:

- add income/time opportunity cost heterogeneity;
- add childcare marketization so high-income households can substitute money for
  time;
- separate child taste heterogeneity from budget constraints;
- discipline with parenthood by age 35 and completed fertility by age 45.

### B. Wealth Versus Income Confounding

The model currently treats income as age-deterministic. Data will separate:

- current income;
- permanent income;
- liquid wealth;
- housing wealth;
- homeownership.

If fertility is hump-shaped in wealth but negative in female earnings or
positive in housing wealth among owners, the model must not collapse these into
one scalar budget state.

Candidate reconciliation:

- add persistent earnings heterogeneity only if data show an income gradient
  conditional on liquid wealth and tenure;
- otherwise keep income external and let wealth/tenure constraints do the work.

### C. Housing Non-Homotheticity Too Strong

The room audit shows the current model bunches renters too heavily and has a
compressed renter upper tail. If \(\bar h(m)\) is doing too much, the model may
match birth-related room increments while distorting fertility-by-wealth.

Candidate reconciliation:

- discipline baseline childless room distributions before family increments;
- estimate housing demand elasticity by parental status in PSID/ACS;
- allow smoother or heterogeneous family-space needs rather than a single hard
  \(\bar h_{\mathrm{jump}}\).

### D. Tenure Wealth Channel Misread

Housing wealth papers imply opposite signs for owners and renters after house
price shocks. A single price effect in the model is insufficient.

Candidate reconciliation:

- compute model fertility by tenure and housing equity;
- run PSID event studies by pre-birth tenure and local rent/price growth;
- keep housing wealth separate from liquid wealth in empirical targets.

## 6. Immediate Model Diagnostics

Add a diagnostic table from each saved benchmark:

1. Distribution-weighted completed fertility by liquid-wealth quintile at ages
   `35`, `40`, `45`, and `50`.
2. Childlessness / ever-parent status by liquid-wealth quintile at the same ages.
3. Fertility by tenure and liquid-wealth quintile.
4. Fertility by location and liquid-wealth quintile.
5. Fertility policy by liquid wealth for renters and owners, separately.
6. Wealth distributions by completed fertility:
   \[
   E[b\mid n^\ast=0],\quad E[b\mid n^\ast=1],\quad E[b\mid n^\ast\geq 2].
   \]
7. Ownership and housing-equity distributions by completed fertility at older
   ages.

The critical distinction is policy versus distribution:

- policy plot: what a household would choose at a given state;
- distribution moment: what the model predicts in the stationary population.

Only the latter should be compared to PSID/ACS facts.

Implementation added:

- `calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/plotting/export_income_wealth_fertility_diagnostics.m`
- `docs/archive/model_status_2026-05-07/INCOME_WEALTH_FERTILITY_MODEL_DIAGNOSTIC_2026-04-25.md`

This exporter loads a saved model anchor and writes a long CSV by age bin,
wealth concept, group, and wealth-bin count. It reports parent share,
childlessness, mean completed fertility, ownership, location, income, and
housing rooms. It distinguishes liquid wealth, sale-net wealth, and accounting
net worth.

First run result: the current benchmark's fertility-by-liquid-wealth profile is
not a clean monotone budget gradient. It is strongly shaped by tenure, mortgage
debt, and center-periphery sorting. This makes tenure/location splits mandatory
before interpreting non-homotheticity.

## 7. Immediate Data Roadmap

### A. Consolidate Existing PSID Wealth-Fertility Outputs

Build one audited CSV with rows:

- wealth concept: liquid net worth/income, total net worth, total net worth
  excluding housing, housing wealth, current income, permanent income proxy;
- age outcome: parent by `35`, parent by `40`, parent by `45`, children ever
  born at `45`;
- group: all, women, men/household head if appropriate, renters at `25-30`,
  owners at `25-30`;
- bins: quintiles and deciles;
- controls: none, cohort, education, marital status, metro.

This can reuse the `fertility_wealth_v12_simple.do` definitions but should
write one clean `income_wealth_fertility_master_v1.csv`.

Implementation added:

- `code/data/psid_followup_mar2026/income_wealth_fertility_master_v1.do`

This script builds the master PSID sample and writes long-format bin summaries
for parenthood, childlessness, children-count outcomes, and ownership by young
adult income and wealth concepts.

### B. PSID Local-Housing Interaction Tests

Use the `docs/prompts/PSID_TEST_PROMPT.md` priority tests:

1. First-birth flight to lower-rent areas by pre-birth wealth.
2. Housing demand income elasticity by parental status.
3. Price elasticity by parental status.
4. First-birth hazard by pre-birth housing slack.
5. Housing expenditure share jump at first birth by wealth.

The local-housing merge should use HUD FMR by year and MSA if feasible; ACS
median gross rent is the fallback. The key outcome is whether low-wealth parents
move, delay, or compress housing when family-space needs bind.

### C. ACS / CPS / Census Cross-Section

Use ACS/CPS for broad, high-powered facts:

- children in household by income, age, tenure, location;
- completed-fertility or children-ever-born from Census/CPS fertility
  supplements where available;
- room distributions by income, parent status, tenure, and MMS location;
- rent burden by income and parent status.

Caution:

- ACS `NCHILD==0` is not completed childlessness;
- use it for current household children and housing demand, not completed
  fertility;
- for childlessness and completed fertility, use children-ever-born concepts.

### D. PSID Versus ACS Division Of Labor

Use PSID for:

- wealth;
- lifecycle fertility timing;
- event studies around first and second birth;
- local housing conditions around births;
- ownership transitions.

Use ACS/CPS/Census for:

- high-powered cross-sectional gradients;
- room distributions;
- center/periphery location shares;
- income and rent burden by parent status;
- completed-fertility benchmarks where children-ever-born is available.

## 8. Structural Options To Evaluate

Do not change the model before the diagnostics above are available. Once they
are, evaluate these options in order.

### Option 1. Keep Stone-Geary, Recalibrate Needs

Use if the data show a mostly positive wealth-fertility gradient among
constrained renters and no strong high-income reversal after controls.

Change:

- discipline \(\bar h_0\), \(\bar h_{\mathrm{jump}}\), and \(\bar h_1\) with
  baseline room distributions and event-study room increments;
- possibly smooth the child housing threshold to reduce hard kinks.

### Option 2. Add Persistent Income Heterogeneity

Use if fertility varies by income conditional on wealth, tenure, and age.

Change:

- add an earnings state \(z\) with transition matrix;
- allow child costs or opportunity costs to depend on \(z\);
- match fertility by income and wealth jointly.

### Option 3. Add Childcare / Time-Cost Marketization

Use if the data show middle-wealth fertility high but top-income fertility lower
or flatter, especially for women with high earnings.

Change:

- separate child goods costs from parental time costs;
- allow high-income households to buy childcare;
- preserve the housing-space mechanism for renters while preventing a purely
  monotone wealth gradient.

### Option 4. Add Fertility-Taste Heterogeneity

Use only after budget and time-cost channels fail.

Change:

- introduce persistent taste heterogeneity in child utility;
- match childlessness and high-wealth low-fertility facts without forcing
  housing costs to explain all extensive-margin variation.

## 9. Working Acceptance Criteria

A model non-homotheticity revision should not be accepted unless it passes:

1. aggregate fertility, childlessness, and age-first-birth targets;
2. birth-related room increments \(H01\) and \(H12\);
3. baseline room distributions by tenure and parent status;
4. fertility by liquid-wealth bins in PSID;
5. fertility by income bins in ACS/CPS/PSID, with the right fertility concept;
6. tenure-specific housing price effects: negative for renters/non-owners,
   positive housing-wealth channel for owners;
7. old-age parent-childless ownership and wealth gaps.

This list is deliberately stricter than the current calibration target vector.
The goal is not to add all of these as hard moments immediately, but to prevent a
non-homotheticity fix from improving one margin by breaking the distributional
mechanism of the paper.

## 10. Source Checklist

External sources used for this first pass:

- Becker (1960), "An Economic Analysis of Fertility."
- Becker and Lewis (1973), "On the Interaction between the Quantity and Quality
  of Children." NBER chapter:
  `https://www.nber.org/books-and-chapters/economics-family-marriage-children-and-human-capital/interaction-between-quantity-and-quality-children`
- Jones, Schoonbroodt, and Tertilt (2008/2011), "Fertility Theories: Can They
  Explain the Negative Fertility-Income Relationship?"
  `https://www.nber.org/papers/w14266`
- Sommer (2016), "Fertility Choice in a Life Cycle Model with Idiosyncratic
  Uninsurable Earnings Risk."
  `https://www.sciencedirect.com/science/article/pii/S0304393216300745`
- Doepke and Kindermann (2019), "Bargaining over Babies."
  `https://www.aeaweb.org/articles?id=10.1257/aer.20160328`
- Doepke, Hannusch, Kindermann, and Tertilt (2022), "The Economics of
  Fertility: A New Era."
  `https://www.nber.org/papers/w29948`
- Bar, Hazan, Leukhina, Weiss, and Zoabi (2018), "Why Did Rich Families Increase
  Their Fertility?"
  `https://doi.org/10.1007/s10887-018-9160-8`
- Lovenheim and Mumford (2013), "Do Family Wealth Shocks Affect Fertility
  Choices?"
  `https://doi.org/10.1162/REST_a_00266`
- Dettling and Kearney (2014), "House Prices and Birth Rates."
  `https://www.nber.org/papers/w17485`
- Kulu (2013), "Why Do Fertility Levels Vary between Urban and Rural Areas?"
  `https://doi.org/10.1080/00343404.2011.581276`
- Albouy, Ehrlich, and Liu (2016), "Housing Demand, Cost-of-Living Inequality,
  and the Affordability Crisis."
  `https://www.nber.org/papers/w22816`
- Davis and Ortalo-Magne (2011), "Household Expenditures, Wages, Rents."
  `https://doi.org/10.1016/j.red.2010.12.003`
- Census Bureau, "Fertility of Women in the United States: 2020."
  `https://www.census.gov/data/tables/2020/demo/fertility/women-fertility.html`
- Pew Research Center (2015), "Childlessness Falls, Family Size Grows Among
  Highly Educated Women."
  `https://www.pewresearch.org/social-trends/2015/05/07/childlessness-falls-family-size-grows-among-highly-educated-women/`
- Rosen (1979), Roback (1982), Saiz (2010), Baum-Snow and Han (2024).

Local sources used:

- `CALIBRATION_STATUS.md`
- `code/model/README.md`
- `docs/model/LIVE_HOUSING_SIZE_AUDIT.md`
- `docs/prompts/PSID_TEST_PROMPT.md`
- `code/data/psid_followup_mar2026/fertility_wealth_v12_simple.do`
- `code/data/psid_followup_mar2026/output/fertility_wealth_v1/parenthood_by_quintile_v12.csv`
- `output/model/repair_candidate_bridge_best_main_h0down15_uniform_true_strict_strict_eq.txt`
- `output/model/figures_current_candidate/diag09_fertility_by_wealth.png`
