# A pre-decline benchmark, then a fertility shock

Decision report for Tommaso De Santo | 10 September 2026

**I recommend calibrating an approximate pre-2007 economy first, then holding its structural parameters fixed while fitting an announced fertility-preference decline. Policy scenarios should begin from the resulting 2023 distribution.** This is a recommendation for a new empirical design; no target, parameter, model specification or cluster job was changed during this review.

The purpose is exactly the exercise you described: assume that child preferences declined, condition on the transition this produces, and assess housing policies. The fitted preference change absorbs what this maintained model needs to reproduce the selected fertility decline. It does not identify why preferences changed, or establish that preferences caused the decline.

## The design in four steps

**1. Fit the initial benchmark.** Use observations from before the announcement: NCHS 2003-2006 fertility, CPS June 2004/2006 parity, ACS 2005/2006 housing, and PSID 2003/2005 wealth. Keep the reviewed birth-housing response as a pooled structural restriction, with its original estimator and an explicit stability assumption. These windows describe an approximate early economy; none is relabeled as a literal 2007 cross-section.

**2. Assess that fit before freezing parameters.** Show current birth rates, older cohorts' completed fertility, wealth, ownership and housing together. A stationary approximation is acceptable only if its discrepancies are understood and do not drive the policy mechanism. A low optimizer loss alone cannot establish that.

**3. Fit the preference decline with perfect foresight.** Households learn the entire path in 2007. Fit one amplitude to the 2020-2023 average period TFR, 1.64575, using the calendar mapping on page 6. Retain middle windows as checks and keep the fitted path fixed across policies. The existing linear 2007-2023 path with a flat continuation is a parsimonious first specification, not an estimated path shape. [10]

**4. Solve matched policies from 2023.** The policy and baseline must inherit the same historical state and future preference path. The main property-tax comparison still needs its agreed equal-rebate treatment in both regimes, verified market clearing and horizon stability.

## What this investigation changed

**The early data are available.** We constructed housing and wealth candidates, including uncertainty. Earlier statements that these targets were unavailable were too strong: the existing production cache did not contain them. The remaining issues are measurement, the initial stationary approximation and numerical certification.

The literature provides precedents for both initial calibration followed by a transition and joint estimation along a dated path. Neither forces us to fit a 2023 steady state. My preference for the staged design is driven by your policy question, the newly established early-data feasibility and its lower computational cost. [1-5]

# The economic choice and its main limitation

Two decisions have been getting conflated. **The initial-state assumption** determines the distribution of wealth, parity and housing inherited in 2007. **The estimation design** determines whether later observations can change the common structural parameters. Choosing joint estimation does not automatically change the initial-state assumption.

| Design | What it does | Assessment |
| --- | --- | --- |
| Initial fit, then shock | Fit common parameters to early moments; fit the shock on the transition. | Recommended first route. Fast initial diagnosis and a clear separation between benchmark fit and shock fitting. |
| Joint dated transition fit | Fit common parameters and the shock together, comparing every observation with its own date/cohort/window. | Coherent alternative if later outcomes must discipline parameters. A larger estimator; the same stationary initial-state restriction remains. |
| Empirical initial state | Construct the initial joint distribution from data or explicit statistical matching; estimate behavior from flows and subsequent evolution. | Addresses missing cohort histories directly. Requires additional state construction and a new identification design. |

## The fertility discrepancy is real, but it is not a contradiction in the data

Period fertility summarizes births under current age-specific rates. Completed fertility sums births experienced by a particular cohort over its life. In the early data, the period index is **2.0605**; women aged 40-44 in the 2004/2006 CPS report **1.8566** births capped at five, or **1.8784** uncapped. They are different women and histories. [6,9; local fertility receipt]

Under matching populations, exposures and age ranges, a time-invariant fertility schedule makes the full period and cohort integrals equal:

\[TFR(t)=\int f(a,t)\,da,\qquad CF(b)=\int f(a,b+a)\,da.\]

This identity does not literally equate a capped stock at ages 40-44 with full lifetime births. But removing the cap closes only 0.0218 of the gap. At frozen early rates, even the entire remaining age-40-49 fertility integral is only 0.0479 births. We cannot reconcile the gap by capping and incomplete completion alone under that schedule. Nor have we established that postponement explains it.

**My proposed approximation is period-oriented.** Anchor the initial current fertility level; use early childlessness, one-child share and birth timing as additional discipline; always display old completed fertility and higher-parity outcomes as checks. The parity restrictions may force an excessive higher-parity tail. This is a testable risk, not a promise that the initial fit will work.

If that approximation fails on the states that matter for housing and policy, the substantive fallback is a historically heterogeneous initial state. A joint fit with the same stationary distribution changes the compromise; it cannot manufacture the missing pre-2007 histories. [3,6]

# The complete proposed initial target set

Each row states the main parameter discipline, not exclusive identification. Parentheses contain newly measured uncertainty where available. These are candidate observations, not active calibration weights. All housing values use the same 42 MET2013 city codes, a new sample definition explained on the next page.

| Parameter / restriction | Empirical moment and window | Candidate value (SE) |
| --- | --- | --- |
| Initial child preference ψ₀ | Period TFR, NCHS 2003-2006; proposed initial anchor | 2.060500; weight not set |
| First-birth cost F | Childless, women 40-44, CPS 2004/2006 | 0.198279; SE pending |
| First-birth dispersion κE | Mean first-birth age in model midpoint bins, period 2003-2006 | 25.976264; SE/scale pending |
| First-birth dispersion κE | First births at 30+, same period pool | 0.249278; SE/scale pending |
| Further-birth dispersion κC | Exactly one child among mothers 40-44, CPS 2004/2006 | 0.213655; SE pending |
| Patience β | Total net worth / gross labor earnings, PSID 2003/2005 | 6.145861 (0.362855) |
| Bequest strength θ₀ | Annual bequests / wealth; inherited external restriction | 0.008800; historical external |
| Bequest wealth shift θ₁ | Old wealth/income p90/median, ages 76-84, PSID 2003/2005 | 3.515935 (0.306911) |
| Housing supply scale H₀ | Mean min(rooms,9), heads 18-85, ACS 2005/2006 | 5.561097 (0.088381) |
| Owner preference χ | Ownership, heads 30-55, standard structures, ACS 2005/2006 | 0.648334 (0.020675) |
| First-child housing jump hJ | Reviewed first-birth PSID Sun-Abraham contrast, -1 to +3 | 0.720246 (0.085260) rooms |
| Per-child housing floor hC | 3+ minus 1-2 resident-child rooms, ACS 2005/2006 | 0.347067 (0.059705) |
| Joint housing/fertility restriction | Recent-parent minus no-resident-child ownership, ACS 2005/2006 | 0.162896 (0.006080) |

**Count:** 13 restrictions for 10 common structural parameters plus the initial preference level. If the period-rate anchor recovers ψ₀ internally, the remaining objective has 12 rows for 10 searched coordinates. This passes only the count requirement; the local weighted Jacobian must establish informative variation. The shock amplitude is fitted subsequently, not in this initial objective.

**Keep visible as initial checks:** CPS completed fertility 1.856608 capped at five; young ownership 25-34 = 0.431158 (SE 0.021023); old wealth/income median 7.285793 (0.522531). The first two diagnose fertility history and young housing access. The last is a candidate extra bequest restriction if the existing block is weak; adoption requires a new rank check.

# What the new data permit, and what must be aligned

## A genuine early housing sample, with a declared geographic change

The new ACS construction uses the same 42 city identifiers in early and later years. It does **not** preserve the current custom set of admitted geographic areas within those cities. Restoring the old filter reproduces all four current 2023 targets to numerical precision. The old filter excludes 48.9% of the full 42-city head weight in 2012, and 4.54% in 2023. The new construction therefore requires a new target contract and consistent remeasurement of later comparisons.

The IPUMS city assignment itself is approximate because public-use areas change. A fixed city code is not an exact fixed land area. ACS also removed the national nine-room topcode in 2008. For comparable room-level observations, apply min(rooms,9) to both data and model *before averaging*; this changes the observation rule, not the model's housing choices. The reviewed uncapped PSID event-study response retains its own observer. [7,8]

| Same proposed observation rule | ACS 2005/2006 | ACS 2023 |
| --- | --- | --- |
| Ownership, heads 30-55 | 64.8334% | 58.7409% |
| Ownership, heads 25-34 | 43.1158% | 35.6027% |
| Mean rooms, capped at nine | 5.561097 | 5.583723 |

These differences are descriptive, not causal estimates. They show why later ownership levels cannot simply be used as initial levels under an unnamed pooling assumption. The new metro-bootstrap errors describe resampling 42 cities; they are not official ACS replicate-weight survey errors. Full covariance and paired draws are saved.

## Wealth is measurable before the shock

The PSID 2003/2005 pool supplies wealth/earnings and old wealth dispersion with fresh person-cluster bootstrap uncertainty. Both the original long-pool point estimates and their bootstrap errors reproduce the authoritative builders. Wealth/earnings is 6.146 in that early pool, 6.452 in 2005 alone and 7.364 in 2007: the window is economically consequential. Use the saved alternatives as initial-state sensitivities. These are living-household stocks, not estates; distinguish survey-wave dates from income reference years.

## Identification changes must be explicit

The staged proposal replaces later fertility discipline with early period level, timing and parity shape for ψ₀, F, κE and κC; remeasures housing levels and family groups for the housing block; and replaces long-pool wealth/dispersion with early observations for the saving/bequest block. The reviewed childbirth response and external bequest-flow restriction remain. Parameters are jointly disciplined within these blocks. Later outcomes become stated validation under this *new* design; that has not happened to the live objective.

The Sun-Abraham coefficient is an empirical housing response, not itself a deep preference parameter. Applying it to the initial economy requires an explicit assumption that this response transfers across periods; its estimator, sample and reviewed value remain unchanged.

The supply elasticity 0.63, tenure dispersion 0.005 and child-invariant bequests remain maintained restrictions. The recent CEX preference estimate α₀ = 0.733 requires an explicit time-stability assumption; the source of 0.63 remains unresolved. The bequest wealth shift is weakly disciplined in existing evidence. No sensitivity matrix for the proposed system has yet established identification.

National fertility and PSID restrictions combined with metro housing remain a maintained geographic approximation. Choosing this new metro target sample does not by itself certify the population or geographic closure of policy experiments.

# What the primary literature actually establishes

The methodological lesson is to match an observation to the model object that generated it. There is no universal requirement that every estimate come from one calendar year. A maintained preference elasticity and an equilibrium ownership level have different reasons for being pooled. The following distinctions were checked in primary texts.

| Paper | Verified method | What it supports here |
| --- | --- | --- |
| Sommer, Sullivan & Verbrugge (2013), JME [1] | Section 3.5 calibrates four parameters to stationary moments from several vintages. Section 6 starts from an initial steady state, applies an unexpected permanent change in rates/downpayments, then solves perfect foresight. | A close housing precedent for initial calibration plus transition. Their experiment is stylized; it does not jointly estimate a complete annual historical path. |
| Greenwood, Seshadri & Vandenbroucke (2005), AER [2] | Section III.B jointly fits preference/technology parameters and dated household-technology levels to fertility along 1800-1990; market productivity is supplied. | A genuine deterministic fertility-path estimation precedent. Joint fitting is legitimate, but their causal technology story and simpler OLG structure differ from ours. |
| De Nardi, French & Jones (2010), JPE [3] | Section IV initializes from the observed 1996 state distribution, then matches cohort/age/income-specific wealth profiles using simulated moments. | Observed initial distributions are a coherent alternative to stationary initialization. This is a retiree lifecycle model, not our fertility/housing GE transition. |
| Borella, De Nardi, Pak, Russo & Yang (2023), JEEA [4] | Section 6 and Appendix C.9 estimate a cohort lifecycle model under dated tax regimes with explicit perfect foresight; 19 parameters, 448 moments. | Calendar-specific observations can enter one objective. Taxes are measured inputs; this is not joint estimation of an endogenous aggregate price path. |
| Kaplan, Mitman & Violante (2020), JPE [5] | Section III calibrates a stochastic ergodic economy and aggregate regime/belief processes; the boom-bust is a realization of a Markov process. | Useful for disciplined initial micro moments and aggregate drivers. Its expectations structure is not the deterministic announcement assumed here. |
| Kohler & Ortega (2002), Demographic Research [6] | Sections 1 and 3.1 distinguish period/cohort fertility and age-parity exposure; adjusted period rates do not identify actual completed cohorts without additional assumptions. | The old CPS stock and early period index require different observers. A generic appeal to tempo effects does not resolve our measured discrepancy. |

**My judgment:** the staged design is a credible baseline for the conditional policy question. Joint dated estimation is a possible robustness exercise if later observations must help discipline common parameters. Neither is a shortcut around initial-state specification or imperfect empirical measurement. No literature claim here certifies the current implementation.

# Implementation, numerical status and a bounded next run

## Four prerequisites before activating the proposed objective

**Fertility measurement.** The current period diagnostic divides births by adult-household mass, not by correctly aligned female exposure. The annual population bridge does not yet provide the maternal age/parity mapping needed for an exact NCHS comparison. First-birth timing must use period birth flows, while CPS parity uses the relevant completed cohorts and age window. The model's inherited representative count for the 3+ bin also differs from the early CPS count.

**Calendar mapping.** The annual bridge allocates a decision dated t to births in t+1,...,t+4. Under that convention, a 2020-2023 rate window belongs to the 2019 decision block, not the 2023 decision. The historical observer must implement and verify that convention; simply renaming its 2023 row would be wrong.

**Initial demographic closure.** The existing 2.1 normalization also sets births-to-household conversion, initial birth queues and an entry-flow gate. Replacing it with a female period TFR is not a justified one-line substitution. The new fertility anchor and household renewal law must be specified separately. The initial supply construction also inherits elasticity 1.75 while the dated restriction is 0.63; its intended role must be reconciled.

**Family groups and source contracts.** Recent-parent versus no-resident-child data are currently compared with different model groups; resident own children and model dependents are also different. Correct the observers, finish CPS/timing uncertainty, pin the new target/weight fingerprint, then test the exact candidate loop. No missing uncertainty should inherit an unrelated old standard error.

## There has been real progress on the price solver

The inherited-parameter sequential 100-date path now clears markets: maximum residual **0.003673%**, below the unchanged **0.020000%** gate. The root records a zero-difference final replay; collected bytes and fit arithmetic have been checked. The terminal unit-rent gap is still **1.088827%** against a **1%** requirement. Horizon stability, re-estimation and matched policies remain outstanding. Complete inherited fits and parameters are on the final page.

| Work unit | Observed / conditional budget | Interpretation |
| --- | --- | --- |
| Normalized initial candidate | About 5 minutes; four stationary GE solves in the observed case | An initial 21-case sensitivity panel is about 1.7 core-hours, potentially one parallel wave. Not a full calibration. |
| One 100-date policy/value mapping | 42-56 minutes | The entire backward/forward path at supplied prices. |
| One warm-started equilibrium candidate | Roughly 2.2-2.9 hours if three mappings suffice | Optimistic: three calls allow only one price update and a replay. |
| 23-case joint parameter panel | Roughly 51-68 core-hours under the same optimistic root assumption | Several hours with sufficient independent workers; one derivative round, not finished estimation. |

**Next numerical decision:** after the measurement contract is complete, run the small initial stationary panel and fit first; show every target and diagnostic, inspect parity tails and boundary policies, and assess local rank. Only then freeze the common parameters and launch the scalar-shock PF search. Smoke-test each loop, retain per-case checkpoints and stop on failed gates. This research review launched no new numerical jobs.

# Sources and reproducibility

[1] Sommer, K., P. Sullivan and R. Verbrugge (2013). The equilibrium effect of fundamentals on house prices and rents. Journal of Monetary Economics 60, 854-870. Sections 3.5 and 6; printed pp. 860-861 and 867-868. Primary source.

[2] Greenwood, J., A. Seshadri and G. Vandenbroucke (2005). The Baby Boom and Baby Bust. American Economic Review 95(1), 183-207. Section III.B, printed pp. 189-190. Main transition exercise, distinct from Section IV's illustrative steady states. Primary source.

[3] De Nardi, M., E. French and J. B. Jones (2010). Why Do the Elderly Save? The Role of Medical Expenses. Journal of Political Economy 118(1), 39-75. Section IV, printed pp. 46-48. Primary source.

[4] Borella, M., M. De Nardi, M. Pak, N. Russo and F. Yang (2023). FBBVA Lecture 2023. The Importance of Modeling Income Taxes over Time: U.S. Reforms and Outcomes. JEEA 21(6), 2237-2286. Section 6 and Appendix C.9. Primary source.

[5] Kaplan, G., K. Mitman and G. L. Violante (2020). The Housing Boom and Bust: Model Meets Evidence. Journal of Political Economy 128(9), 3285-3345. Section III, especially pp. 3302-3310. Primary source.

[6] Kohler, H.-P. and J. A. Ortega (2002). Tempo-Adjusted Period Parity Progression Measures, Fertility Postponement and Completed Cohort Fertility. Demographic Research 6(6), 91-144. Sections 1 and 3.1, pp. 92-93 and 102-103. Primary source.

[7] IPUMS USA. MET2013 documentation: geographic assignment, mismatch threshold and comparability across PUMA vintages. Accessed 10 September 2026. Primary source.

[8] IPUMS USA. ROOMS documentation, Comparability: national topcode of nine removed in 2008. Accessed 10 September 2026. Primary source.

[9] NCHS. Births: Final Data for 2007, Table 4, printed p. 25. Source of 2003-2006 period rates; 2.0605 is our equal-year average. Primary source.

[10] NCHS. Births: Final Data for 2023, Table 2, printed p. 14. The 2020-2023 average is 1.64575, versus 1.621 in 2023 alone; window averaging is our calculation. Primary source.

## Local evidence

All local evidence is indexed in **output/model/e5f_matched_pf_20260909a/design_research/README.md**. Housing components, bootstrap draws and source receipts are under housing/; wealth equivalents are under wealth/. The fertility extraction and initial coherence review remain under the adjacent parameter_target_audit/ folder. The active isolated model snapshot is 96a41873. The final numerical receipts are under computation/final_replay/.

The lead checked primary papers and source definitions, reviewed empirical masks, independently recalculated saved components and uncertainty, and checked numerical receipt hashes and fit arithmetic. Claude supplied a public-literature review using two research agents in safe mode; its unverified leads were not promoted. No private project material was exported to Claude. Research statements are separated from unimplemented proposals throughout.

# Appendix: inherited parameters at the converged price path

Objective 94.47557608. These are the unchanged inherited inputs and current 12-row objective, **not estimates under the proposed early target set**. The finite price path converged; horizon certification remains open. Shares and gaps use fraction units.

| Moment | Target | Model | Gap | Weight | Loss |
| --- | --- | --- | --- | --- | --- |
| Completed fertility | 1.918 | 1.7771 | -0.140897 | 1425.74 | 28.3037 |
| Childless share | 0.188 | 0.233597 | 0.0455966 | 17180.7 | 35.7197 |
| Mean age at first birth | 26.0446 | 26.2333 | 0.188626 | 44.4444 | 1.58133 |
| First births at age 30+ (share) | 0.260327 | 0.244093 | -0.0162343 | 10000 | 2.63551 |
| First-birth housing response (rooms) | 0.720246 | 0.4114 | -0.308846 | 137.565 | 13.1218 |
| Rooms gap: 3+ versus 1–2 children, ages 30–55 | 0.3677 | 0.347411 | -0.0202885 | 2958.51 | 1.21779 |
| Parent ownership gap (share units) | 0.167662 | 0.146588 | -0.0210739 | 14229.6 | 6.31949 |
| Ownership share | 0.575472 | 0.587508 | 0.0120357 | 1207.85 | 0.174967 |
| Mean occupied rooms, ages 18–85 | 5.77997 | 6.32621 | 0.546237 | 11.9732 | 3.5725 |
| Wealth / annual gross labor earnings | 6.8731 | 7.05729 | 0.184193 | 6.28767 | 0.213322 |
| Annual bequests / wealth | 0.0088 | 0.00835438 | -0.000445616 | 5.16529e+06 | 1.02569 |
| Wealth / income dispersion, ages 76–84 (p90/p50) | 3.44811 | 3.34635 | -0.101757 | 56.9598 | 0.589791 |

All 11 free coordinates are inherited rather than re-estimated. ψ₀ is normalized to old completed fertility 2.1, which is distinct from the proposed female period-rate anchor. Bounds and near-bound flags reproduce the existing receipt.

| Parameter | Value | Lower | Upper | Role | Near bound |
| --- | --- | --- | --- | --- | --- |
| β | 0.995277 | 0.94 | 0.9995 | Free input | No |
| κE | 2.16817 | 0.02 | 50 | Free input | No |
| κC | 1.77077 | 0.02 | 50 | Free input | No |
| χ | 1.05373 | 0.1 | 5 | Free input | No |
| H₀ | 13.716 | 0.2 | 80 | Free input | No |
| θ₀ | 0.57035 | 0 | 8 | Free input | No |
| θ₁ | 0.103724 | 0.02 | 16 | Free input | Yes |
| hC | 0.247085 | 0.1 | 1.8 | Free input | No |
| F | 4.61973 | 0 | 8 | Free input | No |
| hJ | 0.469765 | 0 | 0.5 | Free input | No |
| Δψ | -0.321388 | -1.5 | 0.2 | Free input | No |
| ψ₀ | 0.290052 | - | - | Normalized | No |
| ψ2023 | -0.0313363 | - | - | Derived | No |
| Tenure dispersion | 0.005 | - | - | Fixed | No |
| Supply elasticity | 0.63 | - | - | Fixed | No |

The full parameter names, unrounded values, every target weight and loss contribution are preserved in the source CSVs. No inference about fit improvement across target systems is made. Existing ACS date/group problems still apply to these inherited rows; see the main report before interpreting them.
