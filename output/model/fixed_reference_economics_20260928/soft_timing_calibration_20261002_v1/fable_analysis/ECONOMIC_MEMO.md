# Soft financing, interest timing, and early fertility: an independent read

**Revised after lead review: October 3, 2026, 07:06–07:13 New York (revision 1 and audited arithmetic correction).**
The pre-revision memo and analysis files are preserved under `revision1/` with a
SHA manifest; corrected claims are listed in `revision1/REVISION_NOTE.md` and
`revision1/POST_REVIEW_CORRECTIONS.md`. No model solve was run; the saved-array
analysis was rerun after correcting the closing-screen formula, and its figures
were refreshed.

Fable 5.1 analysis, October 3, 2026, 06:06–06:24 New York. Zero model
solves, no production, calibration, Slurm, Git, or manuscript changes. The
supervisor's `ERRATA.md` (age-25 definition as completed interview age in
\([25,26)\); third arm is a new wealth-target arm with entrant wealth unchanged)
was read during the run and is applied throughout. Everything
below is computed from the pinned snapshots in `input_snapshot/` (SHA-256 in
`input_snapshot/manifest.json`), the two winners' exact-repeat saved arrays, the
named packet READMEs, and narrow source reads. Cohort status at launch: complete,
46 of 48 chains verified, chain 23 in each arm without an admissible point.

Reference identities. **Original timing** winner: chain 15, verified loss
18.4453051941, equilibrium price 0.718172. **Alternative (post-interest) timing**
winner: chain 13, verified loss 13.7711314635, price 0.776057. Both use the soft
purchase rule, financed share \(\phi=0.8\), no child-dependent consumption share,
no hard or quarter rule. Neither is adopted. Complete tables:
[original 14-row fit](input_snapshot/original_fit.csv),
[original 31 parameters](input_snapshot/original_parameters.csv),
[alternative 14-row fit](input_snapshot/alternative_fit.csv),
[alternative 31 parameters](input_snapshot/alternative_parameters.csv), and the
[collection readout with both full tables](input_snapshot/collection_readout.md).
Reproduction: `analysis/native_grid_analysis.py` writes
`analysis/out/native_grid_analysis.json` and four supplemental figures; the
standard 17-plot packets of each arm are untouched in their native roots.

## 1. Central claim

The statement "borrowing constraints don't do shit" is **not a tested result of
the current soft-rule calibration**, and nothing in the available evidence
shows that the early-fertility miss is caused by, or curable by, credit. Three
separate things are true at once:

1. In this model the implemented closing screen passes 95–99 percent of young
   childless renters, while 14–24 percent could not buy the four-room rung and
   honor the ending debt floor, and 4–9 percent of young realized owners sit at
   that floor (Section 3). First-birth attempt probabilities vary sharply with
   income and only locally with liquid wealth at the states examined. These are
   descriptive conditioning facts at selected states; they do not prove that
   credit cannot matter. Child utility is tenure-neutral, but credit can still
   reach births through budgets, housing size, transaction costs, wealth
   accumulation and equilibrium prices. A causal test requires a
   current-soft-rule credit-policy counterfactual at the selected point, with a
   solved price and the standard packet; none has been run.
2. Every piece of financing evidence available is either a fixed-price or
   held-coordinate exercise, or a historical hard/quarter-rule run at older
   points under a different, binding closing test (Section 2). There is no
   accepted soft-rule policy result to explain.
3. The age-25 children-ever-born target (0.810) misses by 0.274–0.285 in the
   original arm, the alternative arm, the held-coordinate timing swap, the
   new wealth-target arm, and every historical 80-to-100 financing arm. The
   decomposition in Section 4 shows where: the model is close to the share of
   25-year-old mothers (0.442 versus 0.457) and misses the number of children
   per young mother (1.21 versus 1.77). Conditional on the model's current
   18–21 first-birth flow and 22–25 first-birth hazard, the observer cannot
   exceed 0.676; this is a conditional bound, not a proof that the target is
   unreachable (Section 4).

## 2. What policy result is actually being explained

Precision matters here because the credit claim has been repeated across
contexts with different objects.

| Evidence | Nature | Fertility object and size | Source |
|---|---|---|---|
| Hard/quarter purchase rule, financing 80%→100% | Fixed-price, fixed-parameter, older points (loss 31.28 and the 97.0/51.6 points). Hard rule: closing resources \(A=b+y\ge(1-\phi)Q\); quarter rule: \(A+0.25S_{\mathrm{hist}}\ge(1-\phi)Q\), with \(S_{\mathrm{hist}}\) net saving in that experiment. Renter buyers failing the test are rejected in the historical kernel. | early fertility 0.5312→0.5291; matched 18–21 cell ownership 5.0%→24.7% with first-birth probability 26.04%→26.20% | `purchase_rules_overnight_v1/local_runtime/local_plan.json` lines 148–150 and 3075–3077; `engines/quarter/refactor_lab/engine/kernels.py` lines 314–328 (buyer rejection `bg_b < dpn`) and 912–918 (quarter floor add-on); `independent_assessment_20261001`; Oct 2 audits |
| Hard/quarter permanent 100% steady states | GE, fixed \(H_0\), \(N\) endogenous, older points | hard: first-birth **hazard** among all childless −0.10%, first-birth flow −2.43%, adult population −2.32%, owner share +8.6 pp, price −0.50%; quarter: hazard −0.18%, flow −2.22%, population −2.11%, owner share +9.1 pp; completed fertility pinned at 2.1 by the renewal root (verified in the receipt this session) | `purchase_rules_overnight_v1/mechanism_deployment/permanent_steady_state_comparison.json` |
| Hard/quarter temporary 100% financing paths (one-date \(\phi\) change, 48/64 dates) | Dated transitions at the 97.0/51.6 points | first-period birth flow −0.54% / −0.43% | `purchase_rules_overnight_v1/mechanism_deployment/README.md`, `accepted_*_temporary_h48/h64_*.json`; CALIBRATION_STATUS.md |
| Permanent transitions | All four failed root or terminal gates | none accepted | CALIBRATION_STATUS.md |
| Soft rule, original versus alternative timing at the **same** parameter vector | GE price root solved in each, parameters held | early fertility 0.5276→0.5249; childlessness 0.2022→0.2029; mean age 26.03→26.07; **ownership 30–55 0.665→0.777**; rooms 5.96→6.10 | `smoke_collection/readout.json` |
| Soft rule, independently recalibrated arms | Two SMM searches, no convergence certificate | early fertility 0.5331 versus 0.5338 | this packet |

So: there is **no accepted policy result under the soft rule**. The accepted
permanent steady-state financing comparisons are the historical hard/quarter
runs; accepted temporary dated paths at those older points are also available.
Neither is a soft-rule policy result. In the permanent steady states, the
first-birth hazard changes by less than 0.2 percent; the 2.4
percent fall in the first-birth flow is the 2.3 percent fall in the adult
population under the fixed-supply closure, a population-scale object, not a
fertility response. The strongest "null" is the held-coordinate timing swap:
a change that raises prime-age ownership by 11.2 percentage points changes
age-25 children ever born from 0.5276 to 0.5249 (−0.0027), while mean age at
first birth rises by about 0.038 years. That is a held-parameter comparison with a
re-solved price, not a recalibrated comparison and not a policy transition.

The timing change is economically large. Under original timing a buyer pays
\(R(b+S-Q)\), so mortgage interest accrues on the whole purchase in the purchase
period. Under alternative timing the purchase cost sits outside the interest
factor, a one-period interest waiver worth \((R-1)Q\). For the four-room rung at
the original price, \(Q=4\times0.718=2.87\) and \((R-1)Q=0.237\) in mean annual
earnings, roughly nine percent of a four-year entry-age income of 2.65. A
financing change of that size that leaves the fertility rows within 0.003 at
held parameters is a strong local insensitivity result for that object. It is
not a credit-policy counterfactual under the soft rule, and it says nothing
about borrowing constraints in the data.

## 3. Credit mechanism on the native grid

All numbers below are from the saved arrays of the two winners, inherited-renter
childless households, weighted by the reconstructed pre-fertility distribution.
The reconstruction divides the saved post-fertility childless mass by
\(1-\pi_j\,a(b,z)\), where \(\pi_j\) is the engine's fecundity and \(a\) the
first-birth attempt probability; it reproduces the observer's pre-fertility
childless mass at each age to \(10^{-15}\) (`pre_childless_identity_check`).
Note for the lead: **no saved distribution is pre-fertility.** `g`,
`g_beginning_distribution`, and `g_cross_sectional_wealth_distribution` all
reproduce the post-fertility stock to \(10^{-15}\). Any earlier policy-slicing
that treated one of them as pre-fertility understated young attempt rates by
roughly half (0.13 versus 0.24 at ages 22–25 here), because the households
whose birth occurred this period are missing from the post-fertility childless
stock.

Feasibility tests for the four-room rung, the smallest rung above the parent room
floor \(h_P=2.6\): the implemented closing screen is
\(Rb+y\ge(1-\phi)Q\), equivalently \(b+y/R\ge(1-\phi)Q/R\); the ending floor plus positive consumption implies
\(b+y/R\ge(1-\phi/R+\kappa/R)Q+c_{\min}/R\) with \(\kappa\) the period holding
cost; a cash-in-hand test (the Greaney, Parkhomenko and Van Nieuwerburgh or
Sommer, Sullivan and Verbrugge convention) is \(b\ge(1-\phi)Q\).

| Original timing, childless inherited renters | 18–21 | 22–25 | 26–29 |
|---|---:|---:|---:|
| Pre-fertility mass | 0.0617 | 0.0391 | 0.0283 |
| Share with \(b\le0\) | 0.667 | 0.462 | 0.622 |
| Pass implemented closing screen | 0.971 | 0.951 | 0.992 |
| Can buy and honor ending floor | 0.856 | 0.788 | 0.758 |
| Hold 20% down in cash | 0.128 | 0.154 | 0.088 |
| First-birth attempt probability, mean | 0.273 | 0.239 | 0.268 |
| …if cash-down feasible / infeasible | 0.516 / 0.238 | 0.593 / 0.175 | 0.646 / 0.231 |
| …if ending-floor feasible / infeasible | 0.319 / 0.000 | 0.303 / 0.000 | 0.353 / 0.000 |
| Probability of choosing ownership now | 0.231 | 0.171 | 0.194 |
| Realized owners at the debt floor \(b'=-\phi Q\) | 0.045 | 0.092 | 0.043 |

The alternative arm's closing-screen shares are 0.971, 0.952 and 0.933;
other rows are in the JSON. Attempt
probability by income state at 22–25, original arm, states 1–9:
0, 0, 0, 0.006, 0.305, 0.720, 0.840, 0.879, 0.898, with mass shares
0.006, 0.048, 0.163, 0.294, 0.291, 0.149, 0.042, 0.006, 0.000.

Reading, with populations stated. Denominator for the first three feasibility
rows: pre-fertility childless inherited renters in the age cell. The implemented
closing screen passes 95–99 percent of them in the original arm and 93–97
percent in the alternative arm; 14–24 percent could not buy the
four-room rung and honor the ending floor at positive consumption; 85–91
percent hold less than 20 percent of \(Q\) in cash. Denominator for the floor
row: realized owners of all family states in the cell, stayers and movers;
4–9 percent of young owners end the period at \(b'=-\phi Q\), under 2
percent at ages 38–78. At age 82, the terminal cell, the at-floor share rises
to 30.5 percent original and 70.9 percent alternative; this boundary behavior
needs separate diagnosis. Equality at the floor alone does not show how choices
would change if it were relaxed. The floor is reached in some realized states,
and many would fail a strict cash test; pass rates on one screen do not
measure effects. The households failing the ending-floor test have attempt
probability zero, but they also have low income. In the four lowest income
states, mass-weighted attempt rates are near zero, while wealthy occupied
nodes can have high attempt probabilities. The feasibility split mixes income
and wealth and cannot separate their effects. Within an income state the wealth gradient at these ages
is concentrated between zero and positive liquid wealth: state 5 at 22–25 goes
0.213 to 0.347 from the second to the third wealth tercile; state 6 is flat at
0.68–0.75. The descriptive pattern is **strong income conditioning and local
wealth insensitivity at the examined states**. Whether credit matters for
births is a different question: tenure does not enter child utility, but the
budget, housing menu, transaction costs, wealth path and equilibrium price are
all open channels, and the Oct 2 selling-cost probe shows at least one of them
can flip the sign of a financing effect. A mechanism test needs a
current-soft-rule credit counterfactual, not these cross-sections.

## 4. Early fertility: what the 0.533 is made of

Definitions (per `ERRATA.md`). The target is mean children ever born, capped at
three, among CPS June 2004/2006 women of completed interview age 25, that is
ages \([25,26)\). The model observer takes the \([22,26)\) cell and averages
pre- and post-fertility stocks with post-birth weight 0.875 (uniform birth
timing within the cell). Reproduced from the arrays: 0.5330788 (original) and
0.5338047 (alternative), matching `observers.json`.

| Age 25 | CPS 2004/06 | Original | Alternative |
|---|---:|---:|---:|
| Share with any birth | 0.457 | 0.442 | 0.442 |
| Children per mother (capped) | 1.770 | 1.206 | 1.208 |
| Share with two or more | 0.176–0.323 (bounds) | 0.091 | 0.092 |
| Share with more than three | 0.029 | 0 | 0 |
| Mean, capped at three | 0.810 | 0.533 | 0.534 |

The data bounds on "two or more" follow from the published mean, the any-birth
share and the above-three share; the exact split is not in the target file.

**The extensive margin is close; most of the gap is among mothers.** A
woman enters at 18 and can have at most one birth per four-year cell. Under the
observer's projection she can hold two children at \([25,26)\) only with a first
birth at 18–21 and a second at 22–25; the observer cannot record three. The
model's first-birth probability at 18–21 is 0.268 of all women, and 0.389 of
those mothers have a second birth at 22–25. **Conditional** on that 18–21
first-birth flow and on the 22–25 first-birth hazard of the remaining childless
(as encoded in `analysis/native_grid_analysis.py`, lines 123–131), the observer
is at most 0.676, reached only if every 18–21 mother has a second birth. This is
a conditional bound at the current solution, not an architectural impossibility.
Under the same conditioning, reaching 0.810 would require about 0.349 of all
women to have a first birth at 18–21 (0.268 now) with all of them having a
second birth within four years. Shifting first births earlier is not ruled out
by the mean-age row: the midpoint-weighted mean can be preserved by offsetting
later births, and the 30-plus validation share (0.232 versus 0.249) already
leaves room on the late side. Whether other scored rows and feasibility admit
such a shift has not been established; no impossibility was proven.
That conditional 0.349 scenario would also raise the age-25 share with any
birth to about 0.504, above the CPS 0.457. If the empirical any-birth share
were held at 0.457 and no third child were possible by this age, the 0.810
mean would require a two-child share of about 0.352; even with certain second
births, the early first-birth share would need to be at least
\(0.352/0.875\simeq0.403\). These are accounting scenarios, not a proof that
the target is unreachable.

**The projection convention is not the gap.** Using the post-fertility stock of
the cell (age 26 exactly) gives 0.571; CPS at 26 is 0.923. Using the end of the
next cell (age 30) gives 0.900. The model's children-ever-born profile lags the
data by about one four-year cell through ages 22–30 while matching the mean age
at first birth. The hazard shape is one way to reconcile the two: the model's
first-birth hazard is flat at 0.27–0.29 per cell from 18 to 33 and then falls,
so a midpoint-weighted mean can match while the stock at 25 does not. I did not
verify the empirical age-specific first-birth hazard shape in this session.

**Composition differs by income.** In the model, early births belong to the
upper half of the income distribution (attempt probability 0 in states 1–4,
0.72–0.90 in states 6–9); completed fertility by income type rises from 1.57 to
2.10. Children here carry costs proportional to consumption through the
equivalence scale and the room floor, while the benefit \(\psi m^{1-\gamma}\)
is additive and income-free. Whether US early fertility is concentrated among
low-income women is a literature claim I did not verify here (pointers, not
evidence: Kearney and Levine 2012; Bailey, Guldi and Hershbein 2014). If it
holds, it bears on the child-cost structure and would need its own measurement
design before entering the calibration. Relaxing financing in the executed arms
raised ownership mostly in income states that already had high attempt
probabilities; that is a description of those arms, not a proof that no credit
design could reach the zero-probability states.

**Insensitivity across the executed financial comparisons.** Early fertility is
0.5331 (original), 0.5338 (alternative), 0.5249/0.5276 (held-coordinate swap),
0.5352 (new wealth-target arm, a 35.6 percent lower wealth target with the
entrant distribution unchanged), 0.5312/0.5291 (hard 80/100), 0.528 (soft
reference). The moment moved by less than 0.01 across an eleven-point ownership
change and an eight percent price change in these specific comparisons. That is
evidence of local insensitivity to the financial objects that were varied; it
does not establish that the moment is unidentified by the financial block, and
the comparisons are not a credit-policy counterfactual under the soft rule.

## 5. Counterevidence to the "constraints do nothing" claim

- **The current soft rule has no beginning-wealth-only cash test at closing, but
  the historical rules did impose closing eligibility.** A strict cash-in-hand
  20 percent test fails for 85–91 percent of young
  childless renters (table above), and the Oct 1 tabulation found 72 percent
  of actual young purchases would fail it. The historical hard rule instead
  tested closing resources \(A=b+y\ge(1-\phi)Q\), including current income, and
  the quarter rule tested \(A+0.25S_{\mathrm{hist}}\ge(1-\phi)Q\), with
  \(S_{\mathrm{hist}}\) the historical net-saving object. The historical kernel
  rejects renter buyers below the threshold (sources in the Section 2 table,
  verified this session). Those
  historical closing-eligibility calibrations are the ones in which 80→100 financing
  moved matched 18–21 ownership from 5 to 25 percent with first births
  +0.16 points at fixed price, and the permanent steady states gave hazard
  changes under 0.2 percent. A closing-eligibility constraint **has** been
  tested at older calibrations with their own losses (31.28; 97.0/51.6;
  fresh 88.59/48.32), and showed small fertility responses. What has not been
  run is a credit-policy counterfactual at the **current soft** selected point.
  The pre-revision memo's statement that the hard and quarter rules only
  tightened the ending floor was wrong.
- **Where credit did move in GE, the hazard still did not.** The permanent 100
  percent steady states raise ownership by 8.6–9.1 points and move the
  first-birth hazard by −0.10 to −0.18 percent; the flow and population fall
  together by 2.1–2.4 percent through the closure. This is consistent with the
  null, but it is an older calibration under the hard and quarter rules, and
  the population-scale response is itself an economic result the author may
  not want to dismiss.
- **Channels remain open even with tenure-neutral child utility.** Ownership
  can reach births through the housing menu, the budget, transaction costs,
  the wealth path and equilibrium prices. The Oct 2 probes found the six-room
  rental cap binds for 60 percent of renter parents at a flow cost of 2.3
  percent of consumption, and that setting the selling cost to zero flips the
  sign of the financing effect on births at that older point. Those are
  evidence that the active frictions for parents in those runs were the cap and
  the sale-cost lock-in; they do not show that the down payment is inert in
  general.
- **Reduced-form evidence exists for the mechanism.** Hacamo (2021), per the
  project's reference notes, documents fertility responses to mortgage-credit
  expansion; I did not re-open the paper in this session. A small structural
  response here would speak to the model's channel set, not to that evidence.

## 6. The new wealth-target arm is a different contract

Per the supervisor's errata (seen during the run and applied here): this arm
changes **only the empirical wealth/earnings target**; the entrant wealth
distribution is unchanged. Chain 2, loss 48.170707378 under the **new** wealth
target 4.45838713455674 (business/farm equity, other real estate and vehicles
excluded) with the old weight 7.595098. Model wealth/earnings is 6.5617, gap
+2.1034, contribution
33.60 of 48.17. Every other row sits where the original contract puts it:
ownership 30–55 0.6627, early fertility 0.5352, mean rooms 5.877, first-birth
rooms 1.330, recent-parent gap 0.1271 (within 0.001 of target). Full tables:
[winner fit](input_snapshot/new_wealth/winner_target_fit.csv),
[winner parameters](input_snapshot/new_wealth/winner_parameters.csv),
[results note](input_snapshot/new_wealth/RESULTS.md). This loss is not
comparable to 18.45 or 13.77. The search did not find a point near the new
target; at \(\beta_{\mathrm{annual}}=0.9645\) the model holds too much wealth for
a 4.46 ratio, and no design-certified standard error supports the target or the
weight. It changes nothing about interest timing or early fertility.

## 7. Numerical and measurement artifacts checked

- Occupied region: 99 percent of pooled pre-tenure mass in \([-4.58, 25.73]\)
  on the \([-12, 3000]\) grid; zero endpoint mass in both arms. The young
  childless occupy \([0, 2.8]\) with 46–67 percent at exactly zero. Zero
  endpoint mass is not a grid-convergence certificate. The asset-grid
  diagnosis (older chain 16 / case 0046 point) found occupied-region
  refinement moves ownership by 0.3 points and mean assets by 1 percent; it
  does not certify the newly selected chain 15 and chain 13 winners, whose
  numerical sensitivity to grid resolution is an unresolved check.
- The physical room floor \(h_P\) is at its upper bound 2.6 in the original arm
  and at 2.594 in the alternative. The housing-need lever the optimizer uses to
  chase first-birth rooms is at its search limit; the Oct 1 assessment reports
  no empirical basis for that bound.
- The fertility noise scales (\(\kappa_{\mathrm{fert}}\) 0.108, continuation
  0.356) are flagged "near bound" by the saved screen, which tests proximity
  within 1 percent of a wide search interval, not endpoint contact. That flag
  does not by itself establish that choices are deterministic. The observed
  attempt-probability profile (0 to 0.90 across income states at one age) is
  the direct evidence of sharp income conditioning.
- Saved arrays carry no pre-fertility distribution (Section 3). The observer's
  own accounting is the only pre-fertility record; slicing code should read it.
- Policy probabilities are in \([0,1]\) with no nonfinite entries and no
  occupied negative value steps (`policy_array_summary.json`).

## 8. Limitations

No solve was run, so no fixed-parameter cross-timing policy comparison at the
chain winners exists; the held-coordinate comparison is at the older 23.08 point.
The CPS target file does not report the exact age-25 distribution by number of
children, so the two-or-more share is bounded, not pointed. The income-gradient
claim about US early fertility is from the literature, not re-estimated here.
Attempt probabilities are policies on the native grid; the weights are exact
reconstructions, but "feasible" and "infeasible" groups differ in income, so the
feasibility splits are descriptive, not causal. Entry-queue timing (half at 16,
half at 20) and the 18–21 first cell are taken as the observer takes them.

## 9. Falsifiable next checks

1. **Measure how much of the gap is the one-birth-per-cell projection.** As a
   diagnostic only, not a replacement target: from a birth-history source
   (NSFG, or PSID where early births can be dated), compute the mean at
   \([25,26)\) counting at most one birth per four-year window from 18, beside
   the unrestricted mean. If the restricted statistic falls toward 0.65–0.70,
   the projection accounts for much of the gap; if it stays near 0.80, the
   model's early birth timing is behaviorally too late. Any change to the
   target itself would need a source-consistent measurement design and an
   identification argument; aligning a target to model support is not a
   justification.
2. **Early fertility by income.** Tabulate CPS June children ever born at
   \([25,26)\) by family-income tercile and compare with the model's by-state
   hazards. The model's gradient is strongly positive; the data gradient is
   unverified here. A negative or flat data gradient would bear on the
   child-cost structure.
3. **Make the constraint bind, once.** One bounded fixed-price solve at the chain
   15 point (original timing, price 0.718172, \(H_0=6.7939\)) with the
   cash-in-hand test \(b\ge(1-\phi)Q\) replacing the income screen. Not run
   here, and it cannot be run cleanly without a code change: the income term
   enters the screen through `native_purchase_income`
   (`refactor_lab/engine/household.py:894-901`), and that flag also gates the
   native transaction grid, kernel, and legacy stayer/credit accounting
   (`household.py:682, 703-716`; `shared.py:236`). Switching it off bundles
   other accounting changes. The clean test is a one-line authorized change
   that leaves `dp_choice = ctx.dp_arr` (no income subtraction) with the flag
   on, followed by a solved price and the standard 17-plot packet. This is
   the minimal current-soft credit-policy counterfactual; without it, no
   causal statement about credit and births at the selected point is
   available. Conjecture from the Section 3 cross-sections: young ownership
   falls by more than half and first-birth attempt probabilities change by
   less than 0.01 in every income state. If attempt probabilities move
   materially, the local-insensitivity reading is wrong.

## Files and policy/distribution views

- `analysis/native_grid_analysis.py`: all computations; `analysis/out/native_grid_analysis.json`: full results for both arms.
- Supplemental figures: `analysis/out/F1_age25_decomposition.png`, `F2_ceb_by_age.png`, `F3_first_birth_by_income.png`, `F4_debt_floor_by_age.png`, and `F5_native_asset_policy_mass.png`. F1's dotted lines are the CPS mean plus/minus twice its **working person-bootstrap** standard error, not a design-certified interval. F3 plots raw first-birth attempt policies by income state for childless inherited renters, weighted by reconstructed pre-fertility mass; F4 counts realized owners (stayers via `bp_pol_stay`, movers via `bp_pol`) at the ending floor. F5 uses raw occupied native asset nodes at ages 22–25 for childless inherited renters in income states 4 and 6. Fertility is the pre-tenure first-birth attempt policy; consumption and next-period financial wealth are **conditional on choosing the renter branch**. Its bottom-right panel is reconstructed pre-fertility childless inherited-renter wealth-node mass, normalized separately within each arm and income state and shown on a log scale. The panel keeps nodes with mass above \(10^{-14}\); high attempt probability at a tiny-mass node is not a large population effect. `analysis/plot_native_asset_states.py` and `analysis/out/native_asset_policy_mass.json` reproduce the panel without a solve.
- Standard 17-plot packets, untouched, one per arm. Original timing (chain 15):
  `../collection/production_original_chain_15/run/native_postcheck/selected_postcheck/phase_b_ge/selected_root/standard_diagnostics/`.
  Alternative timing (chain 13):
  `../collection/production_alternative_chain_13/run/native_postcheck/selected_postcheck/phase_b_ge/selected_root/standard_diagnostics/`.
  In each: `policy_childless_renter_age30.png` and `policy_childless_renter_age42.png` (consumption, saving and housing policies by asset node for the childless renter branch; raw saved policies at grid nodes, including initialized infeasible nodes that carry no mass); `wealth_dist_childless_renter_age30.png` and `..._age42.png` (asset-node mass for that cell); `liquid_wealth_by_age_income_state.png`, `housing_by_age_income_state.png`, `ownership_by_age_income_state.png`, `fertility_policy_by_age_income_state.png` (age profiles by income state, realized distribution weights); `housing_market.png`, `market_clearing_by_market.png`, `market_clearing_residuals.png`, `housing_prices.png`, `owner_rungs.png`, `tenure_services.png` (prices, quantities, residuals); `fertility_by_age.png`, `ownership_by_age.png`, `income_state_outcomes.png`.
- No additional native-grid panels were added in the revision; the packets above provide the requested consumption, saving, housing and distribution views for both arms.
