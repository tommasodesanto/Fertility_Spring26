# Reviewer memo, phase 1: calibration targets, entry law and income/purchase timing

**Scope.** This phase was read-only. For startup I read the top of `AGENT_MEMORY.md`, including the September 22 clarification, and the top of `CALIBRATION_STATUS.md`. There are no files under `memory/daily/`. I also read every starting-evidence file you listed, the scored observers and adapters in code, the saved "cell B" fit and parameter tables, and the entry audits. "Cell B" is the selected case from the September 22 earnings/entry battery: one persistent earnings process, the inherited heterogeneous entry wealth, and a finite search. It is not a calibration.

I ran no microdata calculations. Numbers marked "my arithmetic" are computed from saved receipts. Primary sources I could not open are labelled as such.

## 1. Bottom line

1. **No coding error found in the 12 scored observers or the frozen values.** The frozen contract, the scorer and the saved receipts agree. The open problems are mapping and joint coherence, plus weights that read as far more precise than the targets' definitions allow.
2. **The 2.1 normalization and the two CPS rows cannot all hold at the observed family-size mix.** This follows from arithmetic alone and is the most consequential coherence finding (§4.1). The design report already notes that 2.1 differs from the same women's mean of 1.857 children. What it does not state is the consequence: hitting 2.1 while matching childlessness and the one-child share requires about 66% of women with two or more children to be in the model's 3+ group. The same CPS sample has about 46%. Cell B sits almost exactly at that forced corner (≈0.66).
3. **Cell B's entry law is not the July convention.** July converts each entrant's wealth/income ratio into wealth using the entrant's own annual gross income. B instead keeps the reference model's wealth levels and pairs them with the new income ranks. Those levels came from the old 15-state income process, whose age-18 annual income runs from 0.080 to 3.698, a 46-fold range. So B reproduces the July ratios only if the new process gives the same age-18 income at each paired rank. That is checkable tonight without a new solve (Task 1).
4. **The July aggregate-wealth definitions survive the income/timing change unchanged.** The old-age dispersion row keeps its July object, but its denominator mismatch (PSID family income versus the model's pension-plus-transfer income) now depends on the new pension rule. It deserves a level-based robustness check; B fits it (3.73 against 3.52).
5. **The first-birth room target and its model counterpart are not the same quantity** (§4.4). They differ in horizon, starting point, second births and weighting. Still, B's miss (1.34 against 0.72) is larger than the whole range of saved empirical specifications (0.43–0.81). The miss is economic, driven by the parenthood housing floor sitting at its bound, not a measurement artifact.
6. **Several working scales are far tighter than the definitional choices around them.** This holds for the recent-parent gap, the bequest flow, first-birth rooms and, arguably, first-birth timing. Pending author decisions (geography, the first-birth room estimator, the normalization) move the loss by tens of units, which is the size of plausible arm differences. Read tonight's utility comparison row by row.

## 2. What comparable papers measure and match (verified in this session)

| Block | Source | Location | Convention actually used | Bearing on our row |
|---|---|---|---|---|
| Wealth/earnings, bequests | De Nardi & Yang, "Bequests and heterogeneity in retirement wealth," *European Economic Review* 72 (2014) | Table 2 and text, p. 187; timing pp. 185–186 | Wealth / **after-tax** earnings 6.90 (Hendricks 2007a, PSID); bequests / wealth 0.0088 (Gale–Scholz 1994); 90th percentile of bequests / income 4.53 for **single decedents** (Hurd–Smith). Five-year periods, entry at 20, no death before 60 (p. 185). | Verifies 0.0088. The 6.90 is an after-tax object. Their third moment concerns decedents, not the living old. |
| Wealth/earnings, entry, timing, ownership | Kaplan, Mitman & Violante, "The Housing Boom and Bust: Model Meets Evidence," *JPE* 128 (2020) | Table 1 p. 3303; Table 2 p. 3304; pp. 3293–3295, 3305 | Aggregate net worth / aggregate labor income **5.5** (median ratio 1.2); median net worth at 75 / at 50 = 1.51; ownership 0.66, under-35 0.39; owned/rented house size 1.5. Two-year periods. Loan-to-value cap at origination (eq. 5); payment-to-income cap on income at purchase (eq. 6). Initial wealth drawn to mimic "financial assets and its correlation with earnings at age 21" (p. 3305). Earnings AR(1) 0.97/0.20 with **initial SD 0.42**, below the stationary SD. | Closest precedent for an aggregate ratio of sums. They **measure** the entry wealth–earnings correlation; B only assumes one. Entry income dispersion starts below its stationary level. |
| Ownership target | Sommer & Sullivan, "Implications of US Tax Policy for House Prices, Rents, and Homeownership," *AER* 108 (2018) | pp. 256–257; Table 5 p. 258 | Four targets: ownership 0.65 (landlord share, rent-to-wage ratio, share of owners with a mortgage). Ownership is set to the **long-run** 0.65, explicitly not the 2006 peak of 69%. The discount factor targets borrowing behaviour, not wealth. | Contrary precedent to treating 2005/06 as a steady state. Second-order here, because B misses ownership by 19.5 pp. |
| Old-age wealth | De Nardi, French & Jones, "Why Do the Elderly Save? The Role of Medical Expenses," *JPE* 118 (2010) | p. 47, eq. (18); p. 48 | Match **median asset levels** by cohort, age and permanent-income quintile; medians because of cell sizes; cohort and mortality bias reproduced in simulation. | Precedents use levels or medians by income group. None I verified uses a living-old p90/p50 of wealth over family income. |
| Fertility level and timing | Sommer, FEDS WP 2014-32 (**working-paper version** of *JME* 83, 2016) | Table 4, p. 21; Table 2, p. 20 | NLSY79 **same-cohort** completed fertility 1.90 and mean age at first birth 25.5. Earnings persistent AR(1) 0.95 with innovation SD 0.21, plus transitory SD 0.17. | The precedent is internally consistent by cohort. Ours mixes a 2.1 level, 1960–66 cohort stocks and 2003–06 period timing. |
| Childlessness | Baudin, de la Croix & Gobbi, "Fertility and Childlessness in the United States," *AER* 105 (2015) | Abstract-level only | Identified from 1990 Census cohort facts, including the U-shape of childlessness in education. | Cohort-stock convention. Not verified beyond the abstract. |
| Published timing statistic | NCHS Data Brief No. 21 (Mathews & Hamilton 2009) — **agency report** | Main text | Mean age of first-time mothers **25.0 in 2006**; 21% of first births to mothers under 20; 1 in 12 to mothers 35+. | Our 25.98 is a model-cell convention and must not be described as the NCHS mean. |
| Bequest flow, contrary evidence | Feiveson & Sabelhaus, FEDS Note (June 2018) — **report** | Main text | SCF-reported inheritances plus inter vivos gifts ≈ $350B per year (2016$), 1995–2016, ≈3% of disposable personal income. | At a net worth / disposable income multiple of about 6.5 (my assumption, not verified here), that is ≈0.5% of wealth per year. It is lower than 0.88% even though it includes gifts. Survey receipts are underreported. |
| PSID wealth coverage | Cooper, Dynan & Rhodenhiser, Boston Fed **WP** 19-6 (2019) | Abstract page | The standard PSID wealth summary excludes employer defined-contribution accounts such as 401(k)s; adding them moves the median household much closer to the SCF. | The 6.146 target excludes wealth that the model's single asset contains. Direction known, size not measured. |
| Child counts in the ACS | IPUMS NCHILD documentation | Variable page | Counts own children "of any age," including step- and adopted children, living in the household. | The 3+ versus 1–2 row includes adult and step-children. |
| Housing responses to births | Bergsvik, Cools & Hart, *European Journal of Population* (2023) — **demography journal** | PMC text | Norwegian registers; twins and sibling sex mix around second births. Moves concentrate among dwellings of four rooms or fewer. They report moves before first births in other contexts. | Responses depend on the initial dwelling, so the weighting of treated households matters. Supports a baseline before −1. |
| Gale–Scholz primary | *JEP* 8(4) (1994) | **Not retrieved** (AEA 403; mirror refused) | A search-result summary attributes $105.0B = 0.88% of 1986 SCF net worth. **Unverified.** | Whether it counts spousal transfers, and whether it is gross of estate tax, remains open. |

Not verified in this session: Couillard (2025), Greaney–Parkhomenko–Van Nieuwerburgh, Lovenheim–Mumford, Dettling–Kearney, Hacamo, and the Census P20-555 tables (located, not read).

## 3. Cell B under the frozen objective

These rows come from the saved B `target_fit.csv`. The last two model values are rounded from the earlier review table; I did not re-read those two rows.

| Row | Target | B model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| 2.1 normalization | 2.1 | 2.1000029 | 3e-6 | — | unscored |
| Childless women 40–44 | 0.198279 | 0.168899 | −0.02938 | 35,532 | 30.67 |
| Exactly one child, among mothers | 0.213655 | 0.255261 | +0.04161 | 26,953 | 46.66 |
| Mean age at first birth | 25.9763 | 27.2961 | +1.3198 | 139.83 | **243.56** |
| Share of first births at 30+ | 0.249278 | 0.301953 | +0.05268 | 13,866 | 38.47 |
| Wealth / gross earnings | 6.14586 | 6.05925 | −0.0866 | 7.595 | 0.06 |
| Bequests / wealth | 0.0088 | 0.006329 | −0.00247 | 5,165,289 | 31.54 |
| Old wealth/income p90/p50, 76–84 | 3.51594 | 3.72771 | +0.2118 | 10.62 | 0.48 |
| Mean rooms, capped at 9 | 5.56110 | 6.24542 | +0.6843 | 128.02 | 59.95 |
| Ownership 30–55 | 0.648334 | 0.453278 | −0.1951 | 2,339 | 89.01 |
| First-birth room response | 0.720246 | 1.341457 | +0.6212 | 137.57 | 53.09 |
| Rooms, 3+ versus 1–2 children | 0.347067 | ≈0.306 | ≈−0.041 | 280.5 | ≈0.48 |
| Recent-parent ownership gap | 0.162896 | ≈0.118 | ≈−0.045 | 27,056 | ≈54.9 |
| **Total** | | | | | **≈648.9** |

The first-birth timing pair accounts for 43.5% of the loss; the mean-age row alone accounts for 37.5%.

**Searched parameters** (`parameters_actual_bounds.csv`):

| Parameter | Estimate | Actual bounds | Note |
|---|---:|---|---|
| Annual discount factor | 0.9861 | [0.94, **0.99**] | Upper bound is 0.99, not 0.9995 |
| κ_fert | 0.293 | [0.02, 50] | Flagged near bound |
| κ_fert continuation | 0.454 | [0.02, 50] | Flagged near bound |
| χ | 1.009 | [0.1, 5] | |
| H₀ | 7.580 | [0.2, 80] | |
| θ₀ | 0.040 | [0, 8] | Flagged near bound |
| θ₁ | 0.092 | [0.02, 16] | Flagged near bound |
| First-birth fixed cost | 0.315 | [0, 8] | |
| Parenthood housing floor h_P | **2.3** | [0.1, **2.3**] | At upper bound |

**Fixed:** per-child rooms 0 (zero restriction); ψ = 0.171, solved from the 2.1 normalization; payroll tax 0.179; tenure κ 0.005; α = 0.733; σ = 2; supply elasticity 0.63.

## 4. Decision table

Status labels: retain; clarify; targeted repair; author decision; robustness only.

| # | Row | Data estimand (builder, sample) | Model estimand | Status | Confidence | Core reason |
|---|---|---|---|---|---|---|
| 0 | 2.1 normalization | None: author choice. Period TFR 2003–06 averages 2.0605; CPS women 40–44 in 2004/06 average 1.857 capped (1.878 uncapped). | Completed fertility over children-ever-born bins 0/1/2/3+, top bin valued 3.602 (June 2024 CPS vintage). ψ is solved to hit 2.1 exactly. | **Author decision; add a coherence diagnostic** | High | §4.1: forces the 3+ group to ≈66% of women with two or more children, against ≈46% in the data. |
| 1 | Childlessness 40–44 | CPS June 2004+2006, women 40–44, supplement weights; scale is a generalized-variance approximation of the pooled SE (0.0053). | Uniform birth-time projection onto 40–44. | Retain; clarify | High | Correct cohort stock. Part of the miss is mechanical under row 0. |
| 2 | Exactly one, among mothers | Same sample. | Same projection. | Retain; clarify | High | Not a second-birth hazard; same coherence issue. |
| 3 | Mean age at first birth | NCHS first-birth counts 2003–06, ages 12–49, mapped to the model's four-year cell midpoints with tails collapsed. Scale is the 2003–06 annual SD (0.085). | Stationary first-birth flow weighted by cell midpoints. | Retain value; clarify prose; **weight is an author decision** | High (value), medium (weight) | A convention consistent with the model's cells. NCHS published 25.0 for 2006. Also a cohort/period mix. |
| 4 | Share of first births 30+ | Same counts, exact age ≥ 30. | Share of the flow in cells from 30 upward. | Retain | High | Threshold sits on a cell boundary. |
| 5 | Wealth / gross earnings | PSID 2003/05 reference-person families 18–85: net worth over reference person + spouse gross earnings at 18–65; ratio of sums; bootstrap SE 0.363. | Beginning-of-period b + pH over annualized, grossed-up labor earnings. | Retain; clarify | High | July object, September window. PSID excludes employer DC accounts. Not the same object as De Nardi–Yang's 6.90. |
| 6 | Bequests / wealth | External: Gale–Scholz via De Nardi–Yang Table 2. Synthetic 5% scale. | Expected non-negative estates at death, annualized, over beginning-of-period wealth; terminal death imposed. | Retain as external; **tolerance is an author decision** | Medium | Table 2 verified; primary not. Survey-reported flows are lower (§2). |
| 7 | Old p90/p50 of wealth/income, 76–84 | PSID 2003/05 living reference persons, net worth / family income, family income > $1,000, children history observed; SE 0.307. | b + pH over pension + transfer, 76–84 overlap. | Retain + **robustness** | Medium | Family income includes asset income; the model denominator does not (§4.3). θ₁ is weak. |
| 8 | Mean rooms, cap 9 | ACS 2005/06, 42 metros, heads 18–85; SE 0.088. | Realized rooms, capped at 9 before aggregation. | Retain; **geography is an author decision** | High | Cap order verified. National value 5.608. |
| 9 | Ownership 30–55 | ACS 2005/06, 42 metros, structure types 3–10; SE 0.021. | Owner mass, half-weighting the 54–57 cell. | Retain; geography is an author decision | High | The miss is on the model side. 2005/06 is the historical peak (Sommer–Sullivan). |
| 10 | First-birth room response | PSID event study (Sun–Abraham): β(+3) − β(−1), baseline −2; women 18+; never-treated controls; SE 0.085. | Cloned treated/control branches advanced one period; destination rooms uncapped; treated may have a second birth. | **Targeted repair of the observer alignment**; estimand is an author decision | Medium | Horizon, baseline, second births and weights differ. The gap exceeds the full specification range. |
| 11 | Rooms, 3+ versus 1–2 children | ACS heads 30–55 with a child under 18; child count includes adult and step-children. | Current dependents 3+ versus 1–2. | Clarify (low priority) | Medium | Fits within 1 SE. No dedicated parameter in B (per-child rooms = 0). |
| 12 | Recent-parent ownership gap | ACS heads 30–55: oldest own child under 4 versus no own children in the household; 42 metros; SE 0.0061. | Current births from homes with no dependents versus current empty homes, including former parents. | Retain observer; geography is an author decision | High | National value 0.1276 would cut this row's loss from 54.9 to ≈2.5 (my arithmetic, same weight). |
| E | Entry distribution | PSID reference persons 18–24, childless renters, 1984–2019 wealth waves: nonhousing net worth / family income (> $1,000), weighted quintile means; 1,835 family-years. | July: ratio × annual gross income at the entrant's state. B: reference **levels** kept, paired with new ranks. | Keep frozen tonight; **targeted compatibility repair** | High (mechanics), unknown (size) | §5. |
| T | Income/purchase timing | — | Current four-year income, discounted (Y/R), counts toward the down payment; end-of-period debt floor b′ ≥ −φpH. | **Author decision**; robustness later | Medium | Accounting verified. Economically generous with four-year periods (§5). |

### 4.1 The normalization and the CPS rows (model-free arithmetic)

Let $c$ be the childless share, $s_1$ the share of mothers with exactly one child, and $\mu_{2+}$ the mean number of children among women with two or more. Completed fertility $\bar n$ satisfies

$$\bar n=(1-c)\,[\,s_1+(1-s_1)\,\mu_{2+}\,].$$

- **Imposing the targets.** With $\bar n=2.1$, $c=0.1983$ and $s_1=0.2137$, we need $\mu_{2+}=3.06$.
- **The same CPS women.** Their capped mean is 1.857, which gives $\mu_{2+}=2.67$ (2.71 uncapped).
- **Share of the 3+ group required.** With the model's top-bin value $T=3.602$, the 3+ share among women with two or more must be $(3.06-2)/1.602=0.66$. That is about 42% of all women.
- **What the CPS shows.** From the saved unweighted counts, the share is 0.46 in both 2004 and 2006, about 29% of all women (`fertility_availability.json`, lines 40–81 and 144–154).
- **Cell B.** B's two CPS values imply 0.655, assuming the normalization and the CPS projection refer to the same women.
- **Top-bin vintage.** The 2004/06 capped mean within the 3+ group is ≈3.47–3.50 (unweighted, my arithmetic), against 3.602 from the 2024 vintage. The design report already flags this difference (`DECISION_REPORT.md:122`).

Two consequences follow:

- Part of the CPS-row loss (77.3 units) is forced by the choice of 2.1.
- The model's families with 3+ children are counterfactually common. That feeds the family-size housing rows.

This is not an error: 2.1 is an author choice recorded on September 10. The options are:

- keep 2.1 and report the 3+ share as a visible diagnostic, or
- as robustness, normalize to the same-cohort CPS mean using a 2004/06 top-bin value, as Sommer does with 1.90.

### 4.2 Wealth and bequest rows

**Wealth / gross earnings.** The observer matches the July definition exactly (`solver.py:6817-6897`). The residual issues are small or signed:

- **Measurement date.** Data wealth is measured at the interview; the model stock is beginning-of-period. In a stationary economy, $\sum_j m_j(b'_j-b_j)=\text{estates}-\text{entrant wealth}$. So the gap is about half of one period's estates, ≈1.7% of wealth, or ≈0.3 SE.
- **Coverage.** Employer DC accounts are excluded from PSID wealth, so the target sits below the model's wealth concept (size unmeasured).
- **Cross-study range.** Kaplan–Mitman–Violante target 5.5 (SCF 1998, gross labor income). De Nardi–Yang's 6.90 becomes ≈5.3–5.5 on a gross basis at an average earnings tax of about 20–25% (my arithmetic, assumed tax rate). The project's 6.15 (2003/05) and 6.87 (2005–19) sit at the high end, which is plausible given 2003–05 house values.

**Bequests / wealth.** Retain as signed off in July. The 0.00044 scale is about 10 times tighter than the gap between 0.88% and survey-reported flows of ≈0.5%. B's 28% shortfall is within that cross-source band.

Open question: do model deaths start only after age 66? If so, estates of younger decedents are excluded while Gale–Scholz includes all decedents. De Nardi–Yang also rule out death before 60.

### 4.3 Old-age dispersion

The approximation behind the direction claim:

- Suppose family income equals non-asset income $y$ plus asset income $rW$, with a common $r$. Then $W/\text{INCFAM}=f(x)=x/(1+rx)$ with $x=W/y$.
- Because $f$ is increasing, quantiles pass through $f$. The measured ratio becomes $\frac{x_{90}}{x_{50}}\cdot\frac{1+rx_{50}}{1+rx_{90}}<\frac{x_{90}}{x_{50}}$, so the data ratio is compressed.
- Imputed rent from owner-occupied housing is not in family income. That makes the effective return lower at the (housing-heavy) median than at p90, which strengthens the compression.

The lead's point stands that no sign follows once returns and income sources are heterogeneous. Under this approximation, though, the model's pension-only denominator yields larger p90/p50 for the same wealth. B at 3.73 against 3.52 is consistent with that, but it is not proof.

Since precedents use levels (De Nardi–French–Jones medians; Kaplan–Mitman–Violante's 75/50 median ratio), the natural robustness check is p90/p50 of net-worth levels, which needs no denominator (Task 3).

### 4.4 First-birth rooms: equivalence to the data regression

Per the contract, the model observer selects successful first births from the initial distribution. It clones each pre-birth state into a treated branch (first child) and a control branch that stays childless. Both branches advance one period, and destination rooms are compared.

The two estimands coincide only if all of the following hold:

1. **Horizon.** The model's post-birth observation falls about 0–8 years after the birth, mean ≈4, against the data's +3.
2. **Starting point.** The clone starts before the joint decision, so any adjustment made in anticipation of the birth sits inside the model's response. The data target subtracts β(−1), which is small: ≈0.056 on the May curve (`first_birth_correction_review/README.md:24`).
3. **Second births.** The treated branch may have a second birth at destination; the shares must be compared with PSID second births by +3.
4. **Controls.** The model's control is an exact same-state counterfactual; the data relies on never-treated controls and parallel trends.
5. **Weighting.** The model weights by the stationary first-birth flow; Sun–Abraham weights by cohort. Responses vary with the initial dwelling (Bergsvik et al.).

Items 1 and 3 push the model's estimand above the data's. Even so:

- The saved empirical range is 0.43 (binned design) to 0.81 (August Sun–Abraham).
- Re-scoring B's 1.34 across that range gives losses of 40–115 (my arithmetic).
- **Repair the observer's horizon and second-birth treatment; do not expect a closed gap.** The floor at its 2.3 bound is the economic symptom.

## 5. Do the income and timing changes break the July conventions?

- **Aggregate wealth: no.** The numerator (b + pH) and the denominator (gross earnings, annualized once, payroll grossed up) are unchanged and are independent of the new earnings process.
- **Entry: yes, mechanically.**
  - July's convention converts ratios using each entrant's own income (`solver.py:464-484`).
  - B keeps the reference wealth levels instead (`e5f_earnings_wealth_contract.py:16-79`; `selected_B_checkpoint_audit.json:41-75`: mean 0.1865 model units, L1 gap 2e-16).
  - Those levels are the July ratios times the old process's annual gross income. It spans 0.080–3.698 across 15 states, putting entry wealth between −8.22 and 11.48 (`entry_reference_audit.json:26-49`). B's income grid spans 0.136–4.464 on seven states.
  - So the entrants' wealth/income ratios in B are not the July ratios unless age-18 income matches rank by rank.
  - Separately, the 18–24 sample remains a deliberate, sound July repair. The family income versus labor earnings denominator, and the assumed pairing of wealth and income ranks, are measurable (Task 1). Kaplan–Mitman–Violante calibrate exactly that correlation at age 21.
- **Old age: object unchanged.** The model denominator now depends on the new pension rule (`pension_period` 2.046, budget-derived). The denominator mismatch should be re-checked (Task 3).
- **Purchase timing: an author decision, not an error.**
  - Counting current-period income toward a purchase is standard in discrete time: Kaplan–Mitman–Violante (two-year periods, payment-to-income on income at purchase) and De Nardi–French–Jones (income realized, then consume and save).
  - With four-year periods it amounts to a bridge loan of up to four years of income. That weakens the (1−φ) threshold for young households and makes the entry distribution matter less.
  - It does not touch the July definitions, and its equilibrium sign is not established. Because the paper's mechanism is the down-payment constraint, a beginning-of-period-cash sensitivity belongs in the next round, not tonight.

## 6. Joint coherence, weights and identification

**Cohort and period mix.** The system combines:

- a 2.1 level,
- stocks from the 1960–66 cohorts (CPS 2004/06),
- 2003–06 period timing, mostly from the cohorts born 1978–85.

The design report discloses this mix (`DECISION_REPORT.md:37`); §4.1 quantifies what it costs.

**Weights.** The loss is a working minimum-distance criterion: diagonal weights, cross-moment covariance ignored, no efficiency and no overidentification test. Scales mix four kinds of uncertainty: generalized-variance sampling (CPS), year-to-year variation (NCHS), bootstrap SEs (PSID and ACS) and a synthetic 5% (bequests). Comparing the gap between documented alternative values (same object, different window, geography or specification) with each row's scale:

| Row | Alternative values | Gap between them relative to scale |
|---|---|---:|
| Recent-parent gap | 0.163 (42 metros) vs 0.128 (national) | 5.8 |
| Bequests / wealth | 0.88% vs ≈0.5% survey-reported | ≈9 |
| First-birth rooms | 0.43–0.81 across specifications | ≈2–4 |
| Wealth / earnings | 6.15 (2003/05) vs 6.93 (2005/07) | 2.2 |
| Ownership 30–55 | 0.648 vs 0.676 | 1.4 |

Mean age at first birth: the working scale is not a sampling or mapping uncertainty at all. Re-scoring B nationally with unchanged weights gives ownership 116.3, rooms 51.9 and recent-parent gap 2.5, a net change of about −33 (my arithmetic).

**Identification.** Twelve rows for nine searched coordinates proves nothing.

- In B, five coordinates are flagged at or near bounds, so local rank checks there say little about interior identification.
- The bequest pair sits near its lower bounds while the bequest flow is 28% short. That reflects a trade-off or the finite search, not evidence that θ₀ is identified by the flow.
- The three family-housing rows rest mostly on h_P, which sits at its bound. Per-child rooms are fixed at zero, and the recent-parent row has no dedicated parameter.

**For tonight (no job changes).** Report both arms row by row. Separately, re-score the saved model moments under the already-measured alternatives (national housing targets; first-birth rooms at 0.730, 0.805 and 0.4265) as a clearly labelled sensitivity. This tests whether the arm ranking depends on pending decisions. No targets are injected into the running jobs.

## 7. Empirical work order for tonight (at most four tasks, in priority order)

### Task 1: entry-law compatibility

This validates the July entry proxy under the new earnings process.

**Model side (no solve).**
- Inputs:
  - B's 160×7 fixed entry-wealth matrix and age-18 income (`P.income[0,0]`, z grid) from `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/selected_B_checkpoint_audit.json`, and the checkpoint it names;
  - the reference conversion in `.../earnings_entry_battery_v1/entry_reference_audit.json`.
- Calculation: implied entry wealth / B's annual gross income at 18, by B income state and pooled. Report the weighted mean, median and quintile means.

**Data side.**
- Inputs: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/estimate_intergen_income_entry_targets.R`, Block 2 (`make_entry_sample`, `quintile_bin_table`), on `PSIDSHELF_MOBILITY.dta`.
- Sample: identical to July (reference persons 18–24, childless, renters, same waves and weights).
- Calculations:
  - (a) Current nodes.
  - (b) Wealth / reference person + spouse gross earnings, requiring earnings > 1,000.
  - (c) Wealth / group-mean family income, a level normalization.
  - (d) Weighted labor-earnings share of family income (ratio of sums).
  - (e) Weighted Spearman correlation of wealth with income, plus a 3×5 table (earnings tercile × wealth quintile).
- Uncertainty: person bootstrap, 499 draws, seed 20260715.

**Output.** A candidate "entry joint" table (a new estimand, not a correction).

**How to evaluate.**
- B departs from the July convention if its implied mean or median ratio differs from 0.259 / 0.099 by more than about one SE (0.126), or if its quintiles differ by more than 20%.
- The denominator matters if (b) and (a) differ by more than 10% in the middle quintiles.
- The rank coupling overstates dependence if the Spearman correlation is small (for example, below 0.3).

### Task 2: fertility normalization coherence and timing convention

**Inputs.**
- Saved CPS 2004/2006 records via `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/parameter_target_audit/fertility/extract_fertility_availability.py`.
- `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/nchs_natality_timing/first_birth_counts_year_age.csv`.
- B's children-ever-born distribution through `code/model/tools/e5f_initial_fertility_observer.py` on the saved state.

**Calculations.**
- Weighted shares with 0/1/2/3+ children and the weighted capped mean within 3+.
- Implied completed fertility at $T=3.602$ and at the 2004/06 top-group mean.
- The 3+ share required for 2.1.
- B's 3+ share at 40–44 and at completion.
- NCHS: the plain single-year mean (age + 0.5) against the cell-midpoint convention, separating the teen-collapse contribution.

**Output.** A candidate diagnostic table; no target change.

**How to evaluate.** If B's 3+ share exceeds the CPS share by about 10 pp or more, record the tension as an author decision. The NCHS reconciliation should reproduce ≈25.0–25.2 on the published basis.

### Task 3: old-age dispersion robustness

This validates the July living-old object under the new pension rule.

**Inputs.**
- `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/wealth/build_initial_wealth.R` (old-wealth block) on the PSID shelf.
- B's saved distribution, wealth grid, pH and pension.

**Calculations** (2003/05, living reference persons 76–84):
- (i) Current p90/p50.
- (ii) p90/p50 of net-worth levels.
- (iii) Ratio (i) without the children-history filter.
- (iv) If the shelf carries any non-asset income component, the ratio with a non-asset-income denominator.
- Model counterpart of (ii): levels of b + pH at the 76–84 overlap.
- Uncertainty: 499 draws, seed 20260715.

**How to evaluate.** If data and model agree in levels as well as in ratios (within 1 SE), the July object survives. If only the ratio fits, the fit rests on the denominator mismatch, which becomes an author decision on θ₁.

### Task 4: first-birth room alignment (no new regression)

**Inputs.**
- Saved coefficients and covariance under `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/sa_rooms_first_birth_household_aligned_v1/` and `.../first_birth_correction_review/`.
- B's matched-branch outputs from `code/model/tools/e5f_initial_housing_observer.py` (`_stationary_birth_diagnostic`: treated and control means, second births, origin mass).

**Calculations.**
- Data: β(+3), β(+3) − β(−1), β(+4) and the average of β(+2..+4), all with covariance-based SEs; the event sample's share of women with a second birth by +3 and +4.
- Model: the 1.341 split by origin age and tenure; the response excluding second births; origin ages reweighted to the PSID first-birth age distribution.

**How to evaluate.** If the adjustments move the model's estimand by less than 0.2 rooms, the gap is economic and the target stays. Otherwise, align the observer's horizon and second-birth treatment. That repairs the observer, not the target.

## 8. Unresolved questions for the author, and a bounded next phase

**Author decisions:**
1. National versus 42-metro housing targets. The recent-parent row alone moves by 5.8 SE.
2. Keep 2.1 with a visible 3+ diagnostic, or add a same-cohort normalization as robustness.
3. The bequest tolerance.
4. The first-birth room estimand: baseline −2 versus −1, and the horizon.
5. The purchase-timing convention.
6. Whether entry conversion preserves July ratios or reference levels.
7. The mother's age versus the household head's age across rows. Heads are typically older than wives, which shifts the 30–55 windows relative to maternal timing.
8. Whether working scales should carry an explicit floor for definitional uncertainty. That is a weight decision, not a target change.

**Next phase (≈60 minutes).** Adjudicate Tasks 1–4 against the rules above, then verify the three items this session could not open: the Gale–Scholz primary (needs library access), the published *JME* tables of Sommer (2016), and the Couillard and Greaney–Parkhomenko–Van Nieuwerburgh ACS conventions. Candidate amendments for the September 23 target review should be drafted only as labelled candidates, each naming the parameters it affects and its replacement information.

## Sources
- [De Nardi & Yang (2014), EER](https://users.nber.org/~denardim/research/De-Nardi_Yang_EER.pdf)
- [Kaplan, Mitman & Violante (2020), JPE](https://violante.economics.princeton.edu/sites/g/files/toruqf5621/files/documents/kaplan-et-al-2020-the-housing-boom-and-bust-model-meets-evidence.pdf)
- [Sommer & Sullivan (2018), AER](https://kamilasommer.net/Taxes.pdf) ([AEA page](https://www.aeaweb.org/articles?id=10.1257/aer.20141751))
- [De Nardi, French & Jones (2010), JPE](https://users.nber.org/~denardim/research/De_Nardi_French_Jones_JPE_2010.pdf)
- [Sommer, FEDS WP 2014-32](https://www.federalreserve.gov/pubs/feds/2014/201432/201432pap.pdf) ([JME 2016 record](https://ideas.repec.org/a/eee/moneco/v83y2016icp27-38.html))
- [Baudin, de la Croix & Gobbi (2015), AER](https://www.aeaweb.org/articles?id=10.1257%2Faer.20120926)
- [NCHS Data Brief 21](https://www.cdc.gov/nchs/products/databriefs/db21.htm)
- [Gale & Scholz (1994), JEP — not retrieved](https://www.aeaweb.org/articles?id=10.1257%2Fjep.8.4.145)
- [Feiveson & Sabelhaus (2018), FEDS Note](https://www.federalreserve.gov/econres/notes/feds-notes/how-does-intergenerational-wealth-transmission-affect-wealth-concentration-20180601.html)
- [Cooper, Dynan & Rhodenhiser (2019), Boston Fed WP 19-6](https://www.bostonfed.org/publications/research-department-working-paper/2019/measuring-household-wealth-in-the-panel-study-of-income-dynamics-the-role-of-retirement-assets.aspx)
- [IPUMS NCHILD](https://usa.ipums.org/usa-action/variables/NCHILD)
- [Bergsvik, Cools & Hart (2023), European Journal of Population](https://pmc.ncbi.nlm.nih.gov/articles/PMC10115922/)
- [Census P20-555 — located, not read](https://www.census.gov/content/dam/Census/library/publications/2005/demo/p20-555.pdf)
- [PSID-SHELF user guide — income-component list not confirmed](http://fabianpfeffer.com/wp-content/uploads/DaumlerFriedmanPfeffer2025.pdf)
