# Post-presentation issues ledger

## This week's review — September 15, 2026

**Start here.** This is the weekly to-do and to-clarify index. M01–M32 and
the empirical evidence below are preserved; M33–M39 bring in missing work
from the decision ledger, numerical status, and earlier code-review agenda.
The detailed [theory decision map](ACTIVE_DECISION_LEDGER.md#simplified-theory-discussion-map)
remains the record of theory choices. [Calibration status](../../CALIBRATION_STATUS.md)
remains authoritative for the implemented model and numerical results.

**Proposed work order, not new specification decisions.** Author priorities
are retained: permanent earnings heterogeneity is highest urgency (M31),
geographic consistency very urgent (M33), and policy attribution very important
(M34). The order also reflects dependencies. Reviewing an alternative does not
adopt it. This consolidation launches no new runs and changes no slides.

| Order | To clarify with Tommaso | To do / concrete deliverable |
|---|---|---|
| 1. Earnings — **highest author priority** | Retain three permanent types alongside persistent risk, or compare a persistent-only specification? Why was the original addition adopted? (M08, M31) | Recover the adoption history and reconcile variance/autocovariance sources; design a comparison separating fixed-parameter effects from recalibration. |
| 2. Target contract — **very urgent** | Which geography, first-birth reference/horizon, sample, controls and weighting? Which initial-period cohort moments? (empirical table; M18, M21, M33, M36) | Full target inventory: builder, population, dates, estimator, uncertainty and model counterpart. Review source timing and sample exclusions before selecting revised estimates. |
| 3. Demography and fiscal closure | How do household units, dependent children, parental death, adult entry and replacement fit together? Confirm the population, pension and rebate rules. (M02–M04, M07, M13–M14, M26, M39) | One accounting specification with checks for births, dependents, entrants, deaths, estates, pensions and rebates; distinguish current rules from proposed alternatives. |
| 4. Other structural assumptions | Sequential or simultaneous nested shocks? Retain the first-birth cost, utility form and rental cap? What disciplines tenure dispersion and supply elasticity? (M05–M06, M09–M11, M17, M19, M24, M27, M29–M32) | Short assumption/source sheet: economic role, identification, evidence and proposed test. Keep the annual beta cap at 0.99 unless explicitly revised; assess why the estimate reaches it. |
| 5. Numerical baseline — parallel work | What terminal and horizon evidence supports the claims? Keep successive surprises distinct from an announced path. (M15, M35, M38) | Continue authorized recovery; require accepted stages and inherited states, and separate market/fiscal, terminal-distance and horizon checks. Reconcile existing code audits before targeted repairs or speed work. |
| 6. Calibration and validation | Which decisions are settled enough to recalibrate? Which 2023 moments remain validation? (M17–M18, M21, M25, M36) | Full fit/parameter tables and identification checks; initial and 2023 lifecycle fertility/housing, completed fertility, young ownership, wealth and intergenerational allocation. Compare old/new utility under the same contract. |
| 7. Policy interpretation | How much of the weak response is experiment design, calibration, numerics, or economics? (M34, M39) | Matched baseline/policy from the same inherited 2023 state; trace tax, capitalization and rebate effects, then compare with Coven et al. on a documented common basis. |
| Parallel theory discussion | Which still-open choices and tentative assumptions need author resolution? (M37 and the theory map) | Reconcile the existing proposal and recorded decisions; preserve settled choices, identify conditional claims, and leave parked extensions parked. No edits to the author-owned manuscript. |

**Suggested sequence:** first discuss 1–3 and select the essential tests in 4.
Data provenance, accounting documentation and the already-authorized numerical
recovery can proceed in parallel. Recalibration follows the relevant decisions;
final policy claims also require baseline and numerical validation. These are
dependencies, not a promise that every long run finishes this week.

**Definition of done:** a recorded decision or answer, evidence/artifact, and
remaining limitations. A wording fix does not close a model-validation issue;
a scheduler completion does not certify equilibrium. Presentation-only cleanup
(M12, M16, M20, M22–M23, M28) follows the substance. M01 and the wording part of
M03 are already closed. The empirical table retains its original priorities;
its target-contract reconciliation is a prerequisite to a new calibration.

**Vintage warning:** M13, M23 and M26 contain historical implementations,
qualified below. The older decision ledger's numbered “Active points” also
contains superseded population rules, loss values and a 1.75 supply elasticity.
These are historical evidence, not current defaults. Its legacy open-entry,
migration, balanced-growth and purchase-grant proposals stay parked unless
explicitly reopened. The presentation coordination notes below are historical
ownership records, not a new assignment of this week's work.

## Author decision — September 12, 2026

For the September presentation, use the historical May empirical plots and
original regression specification. Tommaso has not had time to review the
subsequent measurement and design revisions. This is a presentation-version
decision, not certification that the historical measurement is correct.
Suggested disclosure: **Original May specification; measurement revisions
under review.** Do not describe that specification as timing-corrected.

The May rooms graph is reproducible: all 18 saved points match the recognized
original specification to numerical precision. The rooms survey-year assignment
problem is verified against the assembled source panel, rather than merely
conjectured. Preserve both findings when discussing the historical results.

The new annual and binned comparisons remain diagnostics. Do not replace the
calibration target or its weight, relabel the historical graph, or promote a
new specification without reconciling the empirical and model definitions.
Empirical follow-up is deferred; the completed data jobs and their collection
automation are stopped/paused. No new empirical regressions are launched by
this decision. Resume the items below after the presentation review.

## Numbers that must not be conflated

| Object | Rooms | Interpretation |
|---|---:|---|
| May graph coefficient at +3 | 0.796858555 | Original normalization omits both −2 and −6. |
| May graph +3 minus estimated −1 | 0.740737457 | Four-year contrast calculated from that same historical curve. |
| Current pinned calibration target | 0.720246262 | Later August specification; +3 minus −1, different sample/weighting/control construction. |
| Historical slide scalar | 0.66 / 0.664 | Does not equal the saved May graph's +3 coefficient; source remains to be reconciled. |

The comparable-horizon scalars 0.740737457 and 0.720246262 differ by
0.020491195 rooms (about 2.8% of the May contrast). This numerical proximity
does not establish equivalent estimands or validate either specification.
Tommaso's preferred pre-birth reference remains −2; subtracting −1 answers a
different question. A −2-to-+3 comparison spans five calendar years, whereas
the currently pinned model mapping spans four years.

## Deferred empirical work

All rows remain **open / deferred**, except the reproduction result recorded
above. Close each with a reproducible result and an explicit specification
decision; smoothness or statistical significance alone is not a criterion.

| Priority | Issue and established evidence | Required resolution |
|---|---|---|
| 1 | Rooms answers are stored one interview early in the merged panel. Every assigned value in the 352,250-row prepared common sample matches its survey-year source after alignment. | Reconstruct rooms directly from the source-year crosswalk; recover the 52,317 source-observed values not restored by shifting; validate year-specific response codes. Preserve May reproduction separately. |
| 1 | Annual event-time support alternates after PSID becomes biennial. Some cohorts lack the intended −2 reference; May also omits −6. | Choose exact years versus two-year windows, retain a reference before the final pre-birth year, and document cohort support and aggregation weights. The tested baseline −3/−2 is a changed object, not exactly −2. |
| 1 | The first binned diagnostic is much smaller: 149,402 fitted rows versus 345,751 in the prior annual common fit. | Review the restrictions individually: exclude 2019+ dates, exclude no-recorded-birth-year observations, require all six displayed windows, then estimator exclusions. Do not attribute the full difference to binning or data transfer. |
| 1 | The binned driver excludes 101,848 rows without a recorded first-birth year after its date restriction. These may include childless people and unknown histories. | Reconstruct history status before deciding which observations can be controls. Compare admissible last-treated and confirmed-childless designs separately; report their populations and identifying assumptions. |
| 1 | A last-treated 2019 cohort cannot remain an untreated comparison group in 2019 and later. | Either restrict admissible dates/cohort comparisons or choose and justify another control group. Preserve this change separately from outcome timing. |
| 2 | Ludovica code constructs rounded frequency weights and uses them in csdid blocks; its separate Sun–Abraham command reproducing May has no weight argument. August uses direct IW probability weights. | Compare unweighted and correctly defined survey-weighted estimates on identical samples; state the population represented. Do not claim the original analysis generally forgot weights. |
| 2 | August restricts to women who are reference persons/spouses, selects one per household-year, and excludes multi-family-unit dwellings. | Decide the observational unit and population; quantify each restriction, repeated household outcomes, and appropriate weighting/clustering. These are not interchangeable with mechanical corrections. |
| 2 | August defines first biological birth across 20 child records and all available histories, unlike the original first-child field. | Compare birth dates and exclusions person by person; distinguish measurement corrections from a changed biological-parent population. |
| 2 | HOMEOWN has a recovered 41-wave source-year mapping; no evidence justifies applying the rooms shift to it. Original regressions exclude tenure 'neither owns nor rents'. | Numerically validate the mapping, decide the ownership denominator, then rerun ownership with the chosen event design. For ownership transitions, use the previous observed interview rather than calendar-year L. |
| 2 | Moving variables lack a recovered upstream year crosswalk. In 2019 the question covers moves since January 2017 and dates the most recent move; recall conventions vary by vintage. | Recover/validate source mapping, move dates and recall intervals; distinguish most recent move from any move, first-mentioned reason from all reasons, non-movers from missing responses. Bin only after establishing the outcome clock. |
| 2 | Reason code 3 includes expansion/better housing; code 6 includes neighborhood, schools and proximity to friends/relatives. | Check vintage-specific codes and use accurate outcome names. Revisit move gating, missing-to-zero errors and reason-response denominators. |
| 3 | Presentation numbers, empirical contrasts and calibration mapping differ. | Reconcile 0.66/0.664, 0.797, 0.741 and 0.720; choose the intended horizon and population; regenerate the target estimate, covariance-based uncertainty, weight and provenance together before recalibration. |

## Mock presentation feedback — September 11–12, 2026

Merged from `MOCK_PRESENTATION_FEEDBACK.md` on September 14, 2026 so that all items M01–M32 live in one ledger; numbering is unchanged. Original format: Feedback / question, Slide or topic, Next step and owner, Status, Answer / decision. Status lines are as of the merge and several predate the September 13 deck rewrite.

First-pass review task: **Mock presentation: concise conceptual review** (`01a09818-86af-75a1-9ff4-3a396e394c4d`), max reasoning. Replies are advisory; only individually authorized fixes may be implemented.

**Coordination decision:** route ongoing conceptual model questions to **Model questions for the September slides** (`01a098e5-7fbe-79a1-954c-9d71cc8e2763`), starting with M10. Keep the quantitative task focused on runs and quantitative deliverables; do not send it each slide/model question. This presentation task remains slides-only, and clarification work is read-only unless the author separately authorizes implementation.

### M01 — Earnings terminology

Replace “heterogeneous income” with “idiosyncratic earnings risk.”

**Status:** resolved; Tommaso authorized this edit with “do it.”

**Resolution:** On “This Paper,” replaced “Income heterogeneity, earnings risk, and bequests” with “Idiosyncratic earnings risk and bequests.” The model section still identifies permanent income groups separately.

### M02 — Adult aging and death

Do people die randomly and age randomly at any age? Clarify the actual timing.

**Status:** open; concise first-pass review requested. No changes authorized.

### M03 — Child versus adult aging

Children mature stochastically; do adults age deterministically?

**Status:** resolved for author-facing wording. The concern is preserved for review: the stochastic event is the exit from dependency, while adult household age advances deterministically. No implementation verification is claimed.

### M04 — Number of children: m versus n

Do we still need the distinction between children at home (m) and children ever born (n), now that bequests do not depend on n?

**Status:** notation consolidation discussed and deferred; keep children ever born $n$ distinct from dependents at home $m$.

### M05 — Utility specification and literature

The utility function is broadly unjustified and hard to parse. Is it common in the literature and quantitative work? Explain the role of every component: why divide by e(m) outside; why retain the Stone–Geary housing floor; why combine Cobb–Douglas and CRRA; why is the preference for children linear and outside the consumption/housing aggregator? Saverio suggests showing only a generic u and its shape properties in the model section, then specifying functional forms in Quantification.

**Status:** exposition edit authorized and completed. The generic model utility and exact quantitative functional form are now separated; justification and literature review remain open.

### M06 — Budget-constraint slide

The slide is messy; the constraints need to be made clear.

**Status:** timing checked; slide simplification explicitly deferred. No budget-slide or model edits made.

**Timing finding:** households trade housing out of beginning-of-period liquid wealth, then earn/pay the gross bond return on the resulting liquid position. Current earnings/pensions and rebates enter afterward; consumption, rent, maintenance, and property tax are period flows. The ordinary initial down-payment test uses wealth available before current earnings. These conventions support the budget slide's gross-return factors on housing purchases and sale proceeds.

The transition rental-pricing identity is consistent with this convention: current rent plus next-period asset value, less current maintenance and property tax, equals the gross bond return on the current purchase price. This verifies timing, not the economic justification of frictionless rental pricing for constrained owners with transaction costs (M10).

**Evidence:** read-only inspection and source-extracted checks of the optimized solver, tenure kernel, and perfect-foresight rent function; nine interior tenure alternatives, the pre-income down-payment restriction, the dated rental identity, and constant-price nesting passed. The quantitative task independently confirmed that all four source hashes match the fitted-patch terminal contract and stationary-policy preparation. No model solve was run. Grid-clipping incidence and the estate price-date convention were not audited here. Scratch check: `tmp/september_slides_review/budget_timing_check.py`.

### M07 — Fertility and child-aging slide

The slide is messy and confusing. Some material concerns calibration and should not be in the model exposition.

**Status:** event-order clarification checked; slide cleanup and any model change remain unapproved.

**Timing finding:** the household chooses whether to attempt a birth; success raises both children ever born and current dependents. Housing, consumption, and current utility use this post-birth child count. Maturation then determines next period's dependents. Under independent-count maturation, conditional on the parent's survival, $m_{t+1}\sim\operatorname{Binomial}(m_t+d_t,1-\mu)$. The newborn is included in this first maturation draw: there is no minimum childhood duration. With a four-year period, a child can therefore exit dependency by the next model date after birth. This is a substantive approximation to flag, not evidence of a discrepancy between the household problem and forward transition. No change implemented; parental-death accounting remains M13.

**Evidence:** pinned optimized `solver.py` uses post-birth values at lines 2770–2821, current child-count utility adjustments at 2241–2271, and child aging in continuation values at 2521 and 7056–7079. `parameters.py:896` builds the independent binomial transition. The perfect-foresight population path calls `run_e5f_open_population_transition.py:828`, which advances the post-fertility survivor distribution using that same child transition. Inspection only; no numerical solve.

### M08 — Age earnings profile

Is e_a exogenously fixed, and is it really needed?

**Status:** open; concise first-pass review requested. No changes authorized.

### M09 — Housing-supply timing

Is housing supply instantaneous? Confirm.

**Status:** open; concise first-pass review requested. No changes authorized.

### M10 — User-cost equation

Clarify the user-cost equation very deeply, for Tommaso as well as the audience.

**Status:** conceptual explanation received in the dedicated clarification task; author will review exposition tomorrow. The code imposes rental-asset pricing, rationalized by an unconstrained competitive investor earning the bond return after maintenance, tax and anticipated appreciation. No explicit investor portfolio/balance-sheet problem is solved. Household frictions can coexist only under that investor interpretation. No equation or model change authorized.

### M11 — Gumbel shocks, expectations, and constraints

**Weekly clarification, September 15:** the source check below establishes what
the active sequential implementation does. It does not settle the author's
earlier preference for simultaneous nested choice. Reconcile the experimental
comparison and adoption history before deciding whether to reopen that branch;
do not silently switch the running model or restart an already-completed test.

If a household is constrained, it still receives the taste shock: is that appropriate? Clarify the Gumbel shocks and expectations.

**Status:** equation/implementation sequence checked by the dedicated clarification task; no substantive mismatch found and no shock-related pre-run model change indicated. Source inspection only, not a new occupied-state numerical audit.

The active architecture is sequential: observe fertility tastes, choose an attempt, realize birth success, observe housing tastes, then choose feasible housing/saving. The housing kernel applies down-payment/borrowing feasibility before softmax and the log-sum expected maximum. Fertility values use the expected housing values separately for successful and failed conception, matching the displayed expectation of the housing maximum. Historical joint/nested experiments are not the active specification.

For tomorrow's exposition, make the feasible menu explicit and describe the shocks as mean-zero Type-I extreme value. The implemented expected maximum is the unadjusted scale times log-sum-exp; the option value of a larger feasible menu is retained. Subtracting log menu size or silently switching to zero-location Gumbels would change values/incentives. Sources: pinned optimized `kernels.py:715–763`, `solver.py:2668–2820`; slides household equations at `latex/september_14_presentation.tex:323` and `:353`. No edits to those equations or numerical routines.

### M12 — Household-problem slides

The slides are messy, the notation is not transparent, and it is unclear whether m or n is needed. The preceding introductory slide may be unnecessary.

**Status:** open; concise first-pass review requested. No changes authorized.

### M13 — Parental death and children

**Current-branch qualification, September 15:** the retained original-birth-queue
recovery uses endogenous household propagation and a separate recorded-birth
entry queue, without historical age-mass rescaling or immigration. The
head-age bridge/person-headship passages below describe earlier branches.
Current queue accounting does not resolve the missing care/housing assignment
for dependents of exiting parents. Neither the joint-death alternative nor a
care pool should be inferred to be active from the historical proposals below.
Use the live calibration status to identify the implemented branch.

Very important: how do survival s_a and the number of children interact? If parents die, do their children mature? If so, is market clearing/population accounting consistent?

**Status:** substantive accounting gap quantified. Subsequent author decision in the quantitative task is to plan both the retained demographic branch and a joint-death/surviving-maturation branch, preferring the latter if verified. Both baselines must rebate property taxes. See the authoritative `docs/model/e5f_two_closure_overnight_plan.md`; this slides task does not implement or launch model work. The new branch's unit conversion/formation rule and reproduction normalization remain explicit pre-search decisions. Earlier proposal discussions below are preserved as history.

**Measured size, initial stationary state:** using checkpoint `120ffc45...` (the recovered no-rebate history's initial equilibrium), parent-household exit removes 0.005972242 literal dependents per unit household mass every four years: **1.19965% of current dependents**, or **5.18089% of births**. All loss is at parent ages 66+, and 48.9851% comes from the forced final-age exit. These are flows, not a measured unassigned-child stock; the 2023 measurement is reported below. Full age table, source hash, reproduction script and checks: `output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/child_accounting/README.md`.

**Author decision:** retain stochastic maturation to avoid tracking child birth cohorts. Fixed or bounded dependency durations are set aside. The quantitative task has been asked to recommend an accounting repair that preserves stochastic maturation, explains care/housing needs after parental death, and avoids person/entrant loss or double counting. Recommendation only; no implementation authorized. This task remains slides-only.

**New author alternative under consultation:** assume dependents die jointly with their parent household, and derive maturation/adult renewal only from surviving dependents. The raw household transition already removes dependents at parental death and reports survivor-only maturation. The quantitative task is assessing how to reconcile that interpretation with person deaths, the independent birth-entry pipeline, head formation and migration, and the existing child-dependency clock. This is a proposed alternative, not an approved implementation or a change to the slides.

**Quantitative assessment received:** joint death is coherent as a simplified aggregate dependent-renewal closure. Domestic adult renewal must use only surviving maturation, with explicit person-to-household formation/retention and separately counted migration; do not automatically reuse the existing births/2.1 conversion. Literal dependent births, deaths and maturation must use consistent units with any 3+ adjustment. This changes the demographic closure and requires re-solving initial/terminal equilibria, fiscal balance and the historical fit. It is not implemented. Children-ever-born remains a history of births after child deaths; estate rules stay unchanged unless separately chosen.

**Numerical feasibility for tonight/tomorrow, quantitative task assessment:** plausible without new household state dimensions; the operator is inexpensive, while existence/convergence of the new equilibrium and refitting are uncertain. Proposed bounded sequence: saved-policy multistep accounting check with all parallel entry streams and overriding age/head rescaling disabled; reproduction check with an explicit formation/retention/migration rule; one initial and one terminal equilibrium with housing/PAYGO checks; short transition; expand/refit only after those pass. The old fertility-2.1 normalization does not establish reproduction under this closure. Final conversion and closure are not yet chosen; no run was launched by this consultation. Monday wording can describe an implemented simplified renewal process only after checks; current displayed fits still use the previous closure.

**Current-path clarification from the quantitative task:** the successive-surprise implementation uses the observed historical head-age bridge before 2023, a fixed 2023 person anchor, and the annual age/sex person operator thereafter. The raw maturation flow is computed but does not determine the person-law head totals. The person block has child-age/sex marginals and the household block has parent-age/dependent-count marginals, but no joint links. Applying joint death to detailed person ages therefore needs an explicit allocation/matching approximation or additional aggregate cohort bookkeeping. The minimal aggregate-renewal version avoids child-age states but gives up a literal age-18 entry interpretation: memoryless maturation permits first-four-year entry. Reported loss/birth ratios are contemporaneous flow comparisons, not estimated cohort mortality probabilities or policy effects.

**Quantitative task recommendation, awaiting author choice:** dependents of exiting households enter an aggregate publicly financed care pool and retain stochastic maturation. Specify their consumption and housing resources and balance funding, potentially before rebating property-tax revenue. Surviving households' children-ever-born state remains unchanged. Estates follow the existing allocation unless an explicit care-fund claim is chosen; do not allocate assets twice. Parent death neither destroys child persons nor itself creates adult entrants. Pool maturation must be reconciled with the separate adult-entry register without adding a second entrant stream. A scalar pool addresses missing care demand but does not by itself close the existing dependent/person-stock bridge. Care costs, funding, and the entry reconciliation remain design decisions; no implementation approved.

**2023 measurement from the quantitative task:** recovered-history current literal dependents 0.525204862; births 0.093827175; parent-exit dependent loss 0.005310226 per four-year transition. Loss is **1.011077% of dependents** and **5.659582% of births**; 36.787057% comes from forced terminal-age exit, and all losses occur at parent ages 66+. Source: `recovered_sequence/source/readout_2023/dependent_loss_summary.json` within the same overnight result folder. This is a flow, not the untracked pool's stock. Arrays were obtained during the task's separately authorized fixed-coordinate observer replay; no extra solve for M13.

**Alternative considered and set aside:** bounded child dependency would address the measured source of the gap. Last birth is at model age 42; a maximum 20-year dependency duration ends before parental mortality begins in this specification. A 50/50 departure at 16 or 20 years illustrates how to preserve the current 18-year mean, while removing immediate maturation and the unbounded tail. This would need child-age/cohort information and renewed fit assessment. If younger-age parent mortality were introduced, explicit guardian/care assignment would still be necessary.

$s_a$ is the parent household's survival probability to the next age cell. The household transition multiplies the entire post-birth distribution by $s_a$ before maturation. Dependents attached to dead households therefore leave the household distribution; they are neither reassigned nor included in that transition's reported mature-child flow. At the terminal parent age, the entire household cohort exits.

This does **not** imply that the fitted transition deletes their future adult entry: the separate recorded-birth queue continues independently of parental survival. The newer person-demography operator instead advances child persons using their own age/sex survival and determines household-head totals from headship rates. Neither mechanism, in the inspected implementation, explicitly places surviving children of dead parents into another household or reconciles those children's ongoing housing needs with household dependent states. Thus aggregate demographic/household-head identities do not certify that dependent-child accounting link. Equilibrium and welfare effects of a repair have not been measured.

**Evidence:** pinned `run_e5f_open_population_transition.py:847–895` (survivor-only dependent transition and maturation; terminal exit); `run_e5f_perfect_foresight_transition.py:650–694` and `advance_birth_vintage_queue` at `run_e5f_open_population_transition.py:351` (separate future-entry queue); `run_e5f_perfect_foresight_person_demography.py:807–838` and `demographic_transition/household_person_coupling.py:33–99` (separate person law and age-head-mass reconciliation). The fitted-patch receipt retains the historical head-age bridge and frozen 2023 person anchor; do not represent it as the newer person-cohort operator. The quantitative task independently confirmed the gap: separate person and head accounting gates pass, but no invariant reconciles household-dependent children with demographic children, and there is no explicit orphan reassignment. No solve or implementation change.

### M14 — Equilibrium and population accounting

The within-period equilibrium detail is unclear and probably unnecessary. Give a general equilibrium definition, presumably including the population law of motion. The preceding population slide may be unnecessary. Very important: clarify for Tommaso and the audience the distinction between households, agents, children, and the other population objects.

**Status:** open; concise first-pass review requested. No changes authorized.

### M15 — Speed of convergence

Why does convergence to a new equilibrium take so long? Why not one full generation, or two?

**Status:** open; concise first-pass review requested. No changes authorized.

### M16 — Equilibrium exposition as a whole

The sequence of slides around the equilibrium concept is confusing and needs substantial clarification.

**Status:** open; concise first-pass review requested. No changes authorized.

### M17 — Interest rate and calibration

Is i too low? More generally, the whole calibration needs a much clearer explanation and must be clearer in Tommaso’s own understanding.

**Status:** rate reference checked at the author's request; broader calibration explanation remains open. No parameter change authorized.

The cited source is Greaney, Parkhomenko and Van Nieuwerburgh, *Dynamic Urban Economics*, February 16, 2025 version (local Zotero PDF `MAXJ699L/Dynamic_Urban_Economics.pdf`, physical/printed page 28, paragraph “Interest rate and discount factor”). It sets annual q=0.02 using the average ten-year real interest rate over 1962–2024 in a small open economy. Our July calibration note cites this choice explicitly at `latex/calibration_strategy_eqscale_provisional.tex:103`. The four-year gross factor is 1.02^4=1.08243216. The reference supports the parameter provenance; its underlying real-yield time series was not independently reconstructed here. Our single rate for positive and negative liquid wealth remains a separate simplifying assumption from the numerical value 2%. The author paper has a later December 2025 listing, but this verification concerns the actual February version in the project's reference library.

### M18 — Fertility measurement

Are the slides showing TFR or completed fertility?

**Status:** open; concise first-pass review requested. No changes authorized.

## Model and slide issues — September 13, 2026

Items raised while editing the September 14 deck. Numbering continues M01–M18 above. None is authorized for implementation by the slide task.

### M19 — Tenure smoothing $\kappa_H$

`tenure_choice_kappa = 0.005` is the lower search bound set on June 28, not an interior estimate. Searches returned the bound; the June sweep showed that 0.05 changes the economics materially (old-age ownership near 0.90). The only dedicated moment ever proposed was the PSID four-year-ahead ownership Brier score (0.117113, SE 0.002102); the matching simulated-history exercise was never implemented. Removed from all September 14 slides.

**Decision needed:** treat as external numerical smoothing with a stated value, or implement the auxiliary prediction target and estimate it.

### M20 — Parameter-table vintage

Deck tables report $\psi_0 = 0.160$ and the earlier parameter vintage (e.g. $\kappa_1 = 0.279$, $\xi = 0.127$). The September 13 corrected initial packet gives $\psi_0 = 0.149$, $\kappa_1 = 0.338$, $\xi = 0.266$. Refresh all initial tables together or not at all.

### M21 — Recent-parent ownership gap

Removed from the Identification slide and both target tables on September 13. It remains a weighted target in the live run (13 rows; deck shows 12). If it is dropped from the calibration, name the replacement discipline for the housing block; if kept, the paper table must show it.

### M22 — Owner housing grid in exposition

The draft budget-constraint slide states continuous sizes $h \in [0,\bar h]$; the household-problem and state slides still write $h_{t+1} \in \{0,H_1,\dots,H_K\}$. Choose one exposition (discreteness as computational detail) and harmonize.

### M23 — Equilibrium definition loose ends

**Current-branch qualification, September 15:** the retained baseline rebates
property-tax revenue equally. The zero-rebate clause in the historical record
below is inapplicable to that baseline. Reconcile the written fiscal equations
with M39 before any future exposition update.

Sequence definition adopted September 13. Open: the set $\{V_t,g_t,G_t,N_t,P_t,r_t,T_t,\varpi_t\}$ is not named element by element; $N_t$ is first defined inside the definition after the Population and Households slide was cut; the rebate condition $T_t\int dG_t=\tau_t^p P_t\int h_{t+1}dG_t$ should read $T_t=0$ in the no-rebate baseline; $\mathcal A_t,\mathcal F_t$ defined in words only; four appendix person-accounting frames are now unlinked. Supersedes the exposition part of M14/M16.

### M24 — Preferences cross-partial

Slide states $\partial^2 u_t/\partial s_t\partial m_t>0$ on the utility, not on the aggregator, because $\mathcal C_{sm}$ has ambiguous sign in the implemented form (housing requirement raises it, equivalence scale lowers it). Confirm the sign at $\sigma=2$ or drop the bullet.

### M25 — Identification mapping unverified

Slide assigns $\xi \to$ childlessness, $\kappa_1 \to$ first-birth timing, $\kappa_C \to$ one-child families, $h_P \to$ both room moments. This is a reading of the parameter table, not a Jacobian check. Verify against the local weighted Jacobian before the paper.

### M26 — Population law in slides versus code

**Historical entry, superseded for the current recovery on September 15:** the
description below concerns the earlier imposed-age-mass/person-headship branch.
The retained original-queue exercise carries households forward endogenously,
uses the recorded-birth queue with births divided by 2.1 for future household
entry, and has no immigration or age-mass rescaling. The 2.1 divisor is an
external replacement normalization, not implemented child mortality. A national
headship conversion used to illustrate results is not the model's population
law. Reconcile this distinction with M13, M18 and M33.

The fitted 2007–2023 history imposes observed household age masses (births/2.1 entry queue rescaled to data); the post-2023 forecast converts annual persons to heads with fixed 2023 ACS headship. The deck now says neither. Decide how much to state on the calibration slide.

### M27 — Underwater-debt rollover

No unsecured credit line ($\lambda_d=0$). Debt below the collateral floor arises only after a price fall or a sale with shortfall and rolls over at share $\lambda_{a+1}$ (1 before age 42, linear to 0 at 62). Removed from the deck as niche; belongs in the paper appendix.

### M28 — Schematic frames overfull

Initial Steady State, Impact, and Demographic Adjustment each overflow by 18pt (minipage heights). Cosmetic.

### M29 — Supply elasticity source

The deck cites Baum-Snow and Han (2024) for $\eta = 0.63$. Their headline is an average urban floor-space supply elasticity near 0.5; the independent quantitative audit found no primary receipt deriving 0.63. Establish the derivation or change the value/citation before the paper.

### M30 — First-birth fixed cost $\xi$

The author was not aware the model carries a one-time utility cost at the first birth (`first_birth_fixed_cost`, default zero, estimated at 0.127 in the deck vintage and 0.266 in the September 13 corrected initial). It exists because the first-child housing jump $h_P$ is pinned by rooms moments and childlessness needed its own lever. Decide whether to keep it as a fixed cost of parenthood, replace it with a per-period time/goods cost of children, or test whether childlessness can be matched by $h_P$ alone. Requires recalibration; not for the September 14 deck.

### M31 — Earnings process documentation

**Author priority: highest urgency; substantive retain/remove decision open.**
The author prefers investigating removal of the permanent component, not just
omitting it from exposition. Reconstruct why it was introduced, check whether
the empirical components are compatible, and separate a fixed-parameter
comparison from recalibration. No removal is adopted here. The full rationale
and historical evidence are in the decision ledger's September 14 earnings entry.

Correction (September 14): the live calibration layer (`code/model/intergen_eqscale_seq_optimized/local_panel.py`) builds the earnings state as a five-point Rouwenhorst discretization of an AR(1) crossed with three permanent income types; the three-point grid with persistence 0.85 in `parameters.py` is a module default that the calibration overrides. The slides say only "discretized AR(1)" and omit the permanent types. Open: (i) document the annual persistence/innovation parameters and their four-year conversion in one place with a source; (ii) verify that the PSID-based permanent-type variance and the literature-based persistent process do not double count dispersion; (iii) decide whether the three permanent types stay in the paper exposition.

### M32 — Sources for the selling cost and the rental size cap

The deck cites Greaney, Parkhomenko and Van Nieuwerburgh (2025) for both the 6% selling cost and the 6-room rental cap. The selling cost is a standard transaction-cost value in that literature; the rental cap is a maintained restriction whose level is contentious (earlier sensitivity work found it non-monotone and load-bearing for the ownership fit). Establish a proper empirical basis for the cap, for example the share of 6+ room units that are renter-occupied in the AHS/ACS, or reframe it as a calibrated object.

## Additional weekly work — reconciled September 15, 2026

### M33 — Geographic consistency of targets and quantification

**Author priority: very urgent. Status: open.** Housing targets and lifecycle
comparisons use 42 selected metropolitan areas; fertility and PSID moments
have different, national scopes. This is one pooled housing market, not a
42-location model. Choose the intended population, reconcile targets and
validation series, housing-supply normalization and demographic closure, then
remeasure mismatches. A full-U.S. population equivalent remains an illustrative
display conversion. Source and detailed checklist: the decision ledger's
September 14 geographic-scope entry. Do not relabel existing estimates national.

### M34 — Weak policy effects: calibration versus economic mechanisms

**Author priority: very important. Status: open.** First establish comparable
tax changes, rebates, starting states, geography, demographic closure and
horizons, and separate numerical uncertainty. Then compare credible calibrations
under a common target contract and perform controlled mechanism comparisons at
fixed parameters. Trace holding costs, capitalization and transfers into young
housing, fertility and population; report interactions. Deliver an attribution
table, including the comparison with Coven et al. A weak response is not itself
evidence of a new mechanism. Detailed scope: the decision ledger's first entry.

### M35 — Transition validity, terminal conditions and solver convergence

**Status: open; numerical recovery already authorized.** Distinguish economic
existence, uniqueness/selection, numerical convergence at a finite horizon,
terminal-state distance and stability when the horizon is extended. The first
long recovery candidate now passes finite-horizon market/fiscal checks; the
complete fitted history and terminal/horizon validation remain outstanding.
Old statements that no positive closed endpoint exists concern other branches;
the retained branch has verified stationary endpoints. Neither fact proves
existence or uniqueness of the desired infinite-horizon transition.

Deliver separate acceptance results for each fitted shock and its inherited
state, a terminal-distance and horizon comparison, and sensitivity to starting
guesses. Document the terminal household continuation, population scale and
supply normalization. Reconcile the existing SSJ/derivative-assisted tests and
their measured timings before proposing another solver. Keep the maintained
successive-surprise fit separate from announced and permanent-shock diagnostics;
four-year periods are not calendar years. Read current results in calibration
status rather than treating old scheduler completion as success.

### M36 — Full calibration contract, utility comparison and validation

**Status: open; consolidates M05, M17–M18, M20–M21 and M25.** Produce one complete
target-and-parameter register, including bounds/external restrictions, provenance,
weights, loss contributions and informative-moment/parameter counts. Explain
the annual beta cap (estimated up to 0.99, not fixed there), bequest parameters,
first-birth cost and housing intercept; assess local identification and boundary
behavior. Do not compare losses across changed objective definitions.

Reconcile which old/new-utility tests already exist; the requested static 2023
comparison must use the same targets, fiscal and demographic rules to isolate
the utility change. It is a separate diagnostic from the maintained 2007 initial
calibration and subsequent transition. Keep all 13 active target rows visible,
including the recent-parent ownership gap omitted from the slides. Assess
completed fertility and the initial cohort distribution with correctly dated
data; period TFR and completed fertility need not coincide during transition.
Report 2007 and 2023 data/model comparisons, young ownership, age-housing
allocation, wealth and the stable policy-function diagnostics. Any extra cohort
targets or revised empirical moments require an explicit new target contract.

### M37 — Simplified theory decisions: retain the existing map

**Status: separate author discussion track, not 18 new unresolved tasks.**
The decision ledger has 18 stable theory IDs. Its register records I1 (baseline
institutions) and U0 (unit conventions) OPEN, W3 (reallocation costs) TENTATIVE,
and W4 (welfare/redistribution) REOPENED. W2, Q0 and Q1 are PARKED; V0 records
completed implementation/reconciliation. Preserve the other recorded decisions
and distinguish an adopted direction from still-conditional theorem scope.
Reconcile the existing proposal's feasible financing/compensation, transition
claims and quantitative interpretation with the author; do not duplicate the
theory register or infer new decisions from this weekly index. Suggestions stay
outside the read-only author manuscript.

### M38 — Code correctness, bloat and computation costs

**Status: existing audit to reconcile, not a new broad audit.** Use the
[completed correctness/efficiency review](e5f_full_code_correctness_efficiency_review_20260905.md)
to classify each finding as fixed, still applicable or superseded in the live
code. Verify model-critical policy reuse, retry paths, budgets, feasibility and
population propagation before trusting reuse or optimizing them. Record
measured times for one household solve, stationary equilibrium, transition
mapping and accepted fit stage. Prioritize demonstrated correctness defects,
then targeted speed improvements with reproduction checks. Broad cleanup belongs
in the isolated code-review branch and must not disrupt production recovery.

### M39 — Pension financing and property-tax rebates

**Status: stationary balance verified in the retained branch; complete fitted
transition/policy validation remains open.** PAYGO pensions use current payroll
revenue rather than an accumulated fund. Document the fixed payroll-tax rule
and the benefit adjustment that balances it; verify the budget at every date
and at the terminal equilibrium. The retained 1% baseline and 2% counterfactual
both use equal property-tax rebates. Keep their fiscal and demographic rules
identical and start policy from the baseline's inherited 2023 households.
Report receipts, expenditure, transfers and residuals alongside outcomes.
A repaired stationary budget does not certify an unconverged transition.

### M40 — Transaction volumes (Guido, September 14 talk)

Question: how does the model do on housing transaction volumes; does it give a
good sense of sales volume? The model produces a sales flow directly (owners
who change size or tenure each period, plus estate sales at death) but no
volume moment is targeted or reported. Deliverable: measure the model's annual
sales-to-stock ratio and its age profile in the 2007 stationary state and along
the transition, and compare with the U.S. existing-home turnover rate (sales
divided by owner-occupied stock, roughly 4–6 percent a year in the 2000s;
source to be pinned from NAR/Census before use). With one-period debt, a
6 percent sale cost and no moving shocks, the model's turnover is likely far
below the data; if so, this is evidence for the missing moving shock (P5/X3 of
the structural review). Diagnostic only until a target contract exists.

**Status:** open; raised September 14, recorded September 16.

### M41 — Housing consumption per person in historical perspective (Guido)

Question, separate from M40: put the model's housing consumption per person in
historical perspective. Deliverable: a short data note on rooms (or square
feet) per person in the U.S. over time (Census/AHS: rooms per household and
household size, 1960s onward), the trend of rising space per person alongside
falling household size, and where the model's 2007 and 2023 rooms per person
sit against it. Flag that the model has one price per room and no quality
dimension, so square-foot and quality growth are outside it. Ties to the
geography question (M33): the historical series is national.

**Status:** open; raised September 14, recorded September 16.

## Evidence and reproduction

- Recognized original code: `/Users/tommasodesanto/Desktop/Projects/Fertility/Codes/code_per tommi_addingcontrolsandfixingthings.do`.
- May rooms graph: `/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/Graphs/rooms_f_c_y_all.png`.
- May saved coefficient table: `/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/Tables/rooms_f_c_y_all_estimates.dta`.
- Consolidated audit: [first-birth correction review](../../code/data/psid_followup_mar2026/output/first_birth_correction_review/README.md).
- Binned diagnostic: [results and sample flow](../../code/data/psid_followup_mar2026/output/first_birth_correction_review/binned_rooms/summary.csv), [sample exclusions](../../code/data/psid_followup_mar2026/output/first_birth_correction_review/binned_rooms/sample_flow.json).
- Source-year ownership construction: `/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/Construction_Files/Code/01 Collect housing variables.do`, lines 95–118.
- Moving definitions: [PSID 2019 family codebook](https://psidonline.isr.umich.edu/documents/psid/codebook/FAM2019ER_codebook.pdf), pp. 54–56; [1984 family codebook](https://psidonline.isr.umich.edu/documents/psid/codebook/FAM1984_codebook.pdf), p. 161.

Presentation asset changes are owned by the existing September presentation
task. This data task communicated the author's choice and the distinction
between the historical graph and the pinned calibration scalar. The ledger
does not certify restoration of other historical figures whose source assets
have not yet been verified.
