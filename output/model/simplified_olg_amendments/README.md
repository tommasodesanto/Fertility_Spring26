# Simplified OLG amendment checks

## Overnight continuation — September 8–9

The author submitted the consolidated review in
[this Pro chat](https://chatgpt.com/c/6aa0cfc4-2494-83e9-a1a7-99f5317cab64)
and authorized focused overnight iteration. The live
[work record](overnight_review.md) identifies the active response, independent
proof checks, next action, and morning deliverables. The existing heartbeat
is active through the morning assessment. The packet-preparation status below
describes the earlier handoff, not the current running status.

## Consolidated review packet, prepared for manual submission

The [complete packet](oracle_consolidated_theory_bundle.md) combines
[the new prompt](../../../docs/prompts/oracle_simplified_olg_consolidated_theory.md),
[current context and model excerpts](oracle_consolidated_theory_context.md),
the corrected general-preferences memo, two prior Pro responses, and the
author's writing conventions. The general branch now starts from gross
utility, with linear child needs treated only as a later specialization.

The packet asks for the full dated planner comparison, fertility effects,
and a funded equilibrium transition without adopting a standalone patience
restriction. It records existing partial proofs and unresolved cap, estate,
existence and convergence issues. The author requested the packet for manual
submission; this review has not been launched by the agent. No main model,
paper, or slide specification was changed to prepare it.

## September 8: housing restrictions and fertility follow-up completed

The [third Pro packet](oracle_housing_fertility_bundle.md), from
[this prompt](../../../docs/prompts/oracle_simplified_olg_housing_fertility.md),
was submitted in the existing chat at 7:53 PM Eastern. The browser displayed
the new attachment `Pasted text(20260908-235319).txt` and Stop answering.
It asks for three connected results: utility and primitive conditions for
housing toward the young; privately chosen fertility after a specified
reallocation; and the same static planner choosing fertility with welfare
valued through current parents only. A short transition bridge is secondary.

The existing dated planner's financial permissions remain fixed. Competitive
old households still cannot borrow; the suggested change to the planner's
old-age restriction was not adopted. The author expects old households to
occupy larger homes, but that competitive ordering remains a result to prove.
A bounded Astra/max preflight checked the fertility differential, child-cost
resource accounting and incumbent-commitment scope. It did not solve the new
tasks. Pro completed the third response in 22 minutes 19 seconds. The
[full response](oracle_housing_fertility_response.md) was captured, preserving
all 130 mathematical expressions, and the monitor is now paused.

Its simpler sufficient housing restriction is alpha>=gamma and
beta Gamma/q >= alpha+vartheta, with Gamma covering both old-estate regimes.
The lead checked the main housing comparison, conditional fertility
differential, joint-planner optimality conditions, and the uncapped joint
fertility argument. This is an initial assessment, not adoption or an
independent audit of every auxiliary claim. The strongest average-fertility
result additionally requires the joint planner solution to be uncapped;
the core housing theorem allows binding caps. The work record contains the
exact conditions and verification limits. No model, note or slide revision
was made.

## September 8: dated allocation is now the active question

The author rejects the stationary-cohort/permanent-transfer result as an answer
to the housing question. The [new follow-up](../../../docs/prompts/oracle_simplified_olg_dated_allocation.md)
asks for a full one-date consumption/housing comparison with fixed fertility,
tenure and future real opportunities. The [paste-ready packet](oracle_dated_allocation_bundle.md)
is intended for the existing Pro chat, which already contains the exact model
source. Equal cohort tenure shares follow from a stationary competitive
reference; they do not impose equal housing quantities. The main remaining
task is to derive the competitive housing ordering or a sharper planner
comparison from primitives.

The author submitted this follow-up in the existing chat. Pro completed it
in 23 minutes 41 seconds; the [complete response](oracle_dated_allocation_response.md)
was captured from the browser on September 8, preserving all 125 mathematical
expressions. It gives separate individual and aggregate sufficient conditions
and an analytical counterexample to an unconditional aggregate claim. The
lead's assessment and independent checks are recorded in
[the work record](../../../docs/model/simplified_olg_utilitarian_work.md).
In particular, its claim of unrestricted beta needs qualification: the price
bound used in its sufficient conditions implies an upper bound when other
primitives are fixed. These are discussion results; the model, note and slides
have not been revised. The monitor was paused after this response; the new
fertility follow-up described above is now the active run.

## September 8: full stationary planner under discussion

The author paused further note/slide changes to settle the full planner with
consumption free and fertility fixed. The [new Pro handoff](../../../docs/prompts/oracle_simplified_olg_stationary_planner.md)
asks for a literature-grounded choice of welfare weights, complete estate and
external-resource accounting, and a stationary comparison at the same cohort
masses. The [paste-ready packet](oracle_stationary_planner_bundle.md) includes
the exact current note as background. Its existing dated welfare criterion is
explicitly subject to review. The author will submit the packet manually.
This handoff is not an adopted model revision or a new theorem.

Pro completed the review in the new chat on September 8. The
[complete response](oracle_stationary_planner_response.md) was captured from
the browser, preserving all 31 numbered equations. It proposes a full planner
with outside estate recipients and entrant remittances, domestic intermediaries,
and fixed external wealth measured before estate payouts. It distinguishes
general failure of stationary utilitarian optimality from a conditional
increase in young housing, and includes a proposed analytical reversal.
These are advisory results; the accounting choices are not adopted and the
full argument has not yet received independent verification. Monitoring is
paused after completion. The note and slides are unchanged.

## September 8: theory slides for the September 10 deadline

Read the [seven-slide theory extract](../../pdf/simplified_olg_theory_slides.pdf)
or [full seminar deck](../../pdf/september_14_presentation.pdf). Both use
[the same source](../../../latex/september_14_presentation.tex). The seven
main theory slides and nine theory appendix slides are the only revised
deck sections; quantitative material is unchanged.

The allocation diagram shows raw marginal housing utility and the local
utilitarian transfer. The two demographic panels show an assumed fertility
decline followed by an assumed positive policy fertility effect, with
population computed from the exact cohort law. They do not establish an
equilibrium policy path, convergence, or positive stationary equilibria for
the revised household model. The earlier model's figures remain historical
artifacts and are not inputs to this deck revision.

Rebuild the diagrams, full deck, and extract without running a model solver:

~~~sh
python3 code/model/tools/build_simplified_olg_theory_slides.py
~~~

The [slide verification receipt](theory_slides_utilitarian_checks.json)
records the illustration checks, source preservation, builds, and visual
inspection. The note and theory slides use bars for cohort averages and
integrals for aggregation. Their household equations and analytical results
are unchanged by this notation revision.

## September 8: separate utilitarian note

Read the [utilitarian discussion PDF](../../pdf/simplified_olg_utilitarian.pdf)
([LaTeX](../../../latex/JMP_DS_suggestions/simplified_olg_utilitarian.tex)):
seven pages of main text and seven pages of proofs/extensions. The environment,
household problems and notation are preserved from the conventional-finance
note below. That note and the author manuscript remain unchanged.

The main direct-allocation result uses \(\beta\ge q\), equal weights on
living households' remaining utility, and explicit primitive sufficient bounds
for constrained owners. It reallocates housing from old to young and raises
utilitarian welfare. The proof does not require compensation. The conditional
aggregate housing optimum requires additional cap/floor restrictions.

The note distinguishes three transfer results: simple welfare-improving
redistribution within old age under the original private choices; a
young-directed program conditional on committed fertility and tenure; and an
Appendix E construction allowing both choices to remain private at
\(\phi=q\). The last uses individual, taste-informed targeting and offsets
births within the young cohort, so aggregate fertility remains unchanged.

The original young-directed grant pair raises a treated owner's fertility
after reoptimization at fixed prices and rebates. The all-date equilibrium
fertility effect and higher limiting population remain open. Dated fertility
tests and finite-date population accounting do not require convergence;
endpoint comparisons require positive stationary limits. Welfare weights and
the authority's information/timing remain proposals for author discussion.

The [work record](../../../docs/model/simplified_olg_utilitarian_work.md) and
[verification receipt](utilitarian_checks.json) separate these conclusions.
Independent derivations and hostile reviews cover:

- [Direct allocation](utilitarian_direct_review.md).
- [Transfers](utilitarian_transfers_review.md),
  [hostile check](utilitarian_transfers_hostile_review.md), and
  [assembled statement](utilitarian_assembled_transfer_review.md).
- [Fertility and paths](utilitarian_fertility_path_review.md),
  [hostile check](utilitarian_fertility_hostile_review.md), and
  [assembled statement](utilitarian_assembled_fertility_review.md).
- [Old-only redistribution](utilitarian_old_redistribution_review.md).
- [Full-choice four-group construction](utilitarian_free_choice_transfer_review.md)
  and its [independent hostile review](utilitarian_free_choice_hostile_review.md).
- [Completed Pro response](oracle_utilitarian_response_capture.txt) and
  [assessment of what was adopted](oracle_utilitarian_review.md). Pro's cap
  simplification and exact-grant fertility proof were independently checked;
  Pro did not review the later four-group construction.

Reproduce the 47 exact symbolic and source-preservation checks with:

~~~sh
python3 code/model/tools/verify_simplified_olg_utilitarian.py
~~~

The checker regenerates the algebra receipt; the archived receipt additionally
records final compilation and visual inspection. No numerical equilibrium or
parameter-neighborhood certificate is used for this note.

## September 8: conventional-finance results

Read the new
[discussion PDF](../../pdf/simplified_olg_conventional_finance.pdf)
([LaTeX](../../../latex/JMP_DS_suggestions/simplified_olg_conventional_finance.tex))
for the complete proposed household model, equilibrium, lifetime housing
reallocation, fertility conditions and transition limitations. This is a
separate proposal; the existing reading note and author manuscript are unchanged.

The stationary welfare condition
\(\beta(1+\gamma+\omega_B)(1-\phi+q\tau^p)\ge\phi-q\)
allows conventional mortgage finance and high patience. Explicit primitive
income and capacity inequalities deliver eligible owners in equilibrium.
The Pareto comparison adjusts the selected young households' own future
consumption, housing and estates while fixing individual fertility and tenure.
It does not establish higher current consumption-equivalent housing marginal
values for the young. Fixed-price fertility results, including tenure choice,
are separate. A convergent policy transition and a higher final population
remain unproved for this revised model.

The six bounded analytical reports are:

- [Lifetime housing and welfare](conventional_efficiency_review.md).
- [Hostile welfare review and required repairs](conventional_welfare_hostile_review.md).
  Its heterogeneity and individual-finance repairs are incorporated in the
  new LaTeX appendix. Its warning that the household algebra alone does not
  prove primitive equilibrium coverage is addressed by the separate next report.
- [Primitive stationary existence and owner coverage](conventional_stationary_primitives_review.md),
  including an explicit nonempty heterogeneous parameter family.
- [Fertility, physical caps and tenure](conventional_fertility_review.md).
- [Mortgage repayment timing](conventional_mortgage_timing_review.md), including
  alternatives that have not been adopted.
- [Transitions and remaining proof requirements](conventional_transition_review.md).

The [verification receipt](conventional_finance_checks.json) records independent
reviews, exact lead algebra checks, source hashes and the ten-page PDF inspection.
The [work record](../../../docs/model/simplified_olg_overnight_work.md) separates
proved results from specification choices and unfinished transition work.

## Conventional finance proposal for discussion

The [household specification and analytical financing conditions](../../../docs/model/simplified_olg_conventional_finance_proposal.md)
develop the author's September 7 request to let current income fund the down
payment and add old-age income. They propose allowing old owners to resize,
while retaining lifetime tenure, both physical size limits, and the existing
fertility and estate preferences. These are proposed changes, not an adopted
replacement for the main note. The document gives the mortgage budget,
an explicit conditional income test, the rental alternative, sufficient
housing conditions, and an exact counterexample to an unconditional housing
response. It has received an independent analytical review; its counterexample
budgets, first-order conditions, and derivative were also checked with exact
rational arithmetic. The September 8 results above supply subsequent
stationary welfare and fixed-price fertility analysis; the original proposal
is retained as the preceding discussion step.

## ChatGPT Pro mathematical review

The [review prompt](../../../docs/prompts/oracle_simplified_olg_transition_math.md)
asks whether the finite two-stage transition can be proved with substantial
renting and positive child goods costs. The
[complete Oracle bundle](oracle_transition_math_bundle.md) contains that prompt
and five named mathematical sources. No quantitative data or session history
is included. Source hashes and the delivery state are recorded in
[oracle_transition_math_status.json](oracle_transition_math_status.json).
After the user signed in, the packet was submitted in the in-app browser
with the visible model set to 6 Pro. Its
[completed answer](oracle_transition_math_response.md) proposes a finite
mixed-tenure certificate using an infinite boundary-value operator.
The retrieved exact verifier passes locally, with output identical to its
archived certificate. The proof, code and original household problems are
checked separately. Read the
[short assessment](oracle_transition_math_assessment.md) first.

The result covers every preference decline and later credit rise up to
\(1/20000\), at any actual baseline intervention date. It retains positive
child goods costs, the rebated property tax and more than 47.5% renting.
Both infinite paths converge; initial fertility and terminal population
have the stated signs. This remains a small reform: at most 80% to 80.005%
financing. It does not imply a welfare gain or higher fertility at every date.

One command replays the exact bounds and checks the original household
problems and stationary derivatives:

~~~sh
PYTHONDONTWRITEBYTECODE=1 python3 code/model/tools/verify_simplified_olg_pro_transition.py
~~~

The [receipt](oracle_transition_math_checks.json) matches Pro's
[archived certificate](oracle_transition_math_expected.json), apart from the
source hash changed by packaging. The
[original source](oracle_transition_math_source.md),
[fixed coefficients](oracle_transition_math_coefficients.json),
[reviews](oracle_transition_math_reviews.json), and
[delivery and verification record](oracle_transition_math_status.json)
are preserved. The reviews explain the boundary rows and the summed tail
matrix used for stationary derivatives. The original source also replayed
under Python 3.11.9; the packaged driver passed under Python 3.9.6.

## September 6–7: supporting extensions

[Housing allocation and demographic adjustment](transition_extensions.md)
contains a wider mixed-tenure stationary theorem, a simple primitive
sufficient condition for local convergence, an explicit finite comparison
in the all-owner limit, and an exact finite household fertility condition.
It also states why the dated housing-misallocation comparison applies along
the paths. Two new supporting proofs extend the positive-cost and finite
mixed-tenure results:

- [Positive child goods costs](positive_child_costs.md) gives explicit
  stationary sufficient conditions and a feasible counterexample in which
  more credit lowers terminal population. The conditions use household
  ratios, with no borrowing multiplier; they are not purely primitive.
- [Finite mixed-tenure transition](mixed_finite_transition_proof.md) proves
  a preference decline and a credit reform at any later baseline date with
  positive costs, tax and about 47.6% renting. Its certified shock bounds,
  \(10^{-11}\) and \(10^{-8}\), are very small. This is an explicit existence
  result, not an economically large reform.

The main note, both figures and earlier proof files are unchanged.
These results do not select a planner institution or a paper specification.

One command checks the new results without an equilibrium simulation:

~~~sh
PYTHONDONTWRITEBYTECODE=1 python3 code/model/tools/verify_simplified_olg_transition_extensions.py
~~~

The receipt is [transition_extension_checks.json](transition_extension_checks.json).
It includes symbolic original-equation checks, exact finite interval bounds,
original household checks and source hashes. The consolidated driver also
runs [the positive-cost checks](../../../code/model/tools/verify_simplified_olg_positive_costs.py)
and [the finite mixed-tenure certificate](../../../code/model/tools/verify_simplified_olg_mixed_finite.py).
The research reports, independent reviews and lead checks are preserved in
[transition_extension_reviews.json](transition_extension_reviews.json).
Broad conditions and economically large mixed-model reforms remain open.
The stronger Pro bounds above have their own verification receipt.

## September 6: integrated note and two-stage illustration

Start with `output/pdf/simplified_olg_amendment_proposal.pdf`: seven main-text
pages and six appendix pages. Source:
`latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex`.
It combines the original model and compensation argument with a preference
decline and a later credit reform from the inherited baseline state. The
original separate utilities, all four W/V problems and tenure display are
preserved. The detailed local sequence argument is in Appendix B.

The existing allocation figure is reused unchanged. The new two-panel
`combined_transition_figure.pdf/png` shows the old steady state, the baseline
decline, the two policy-date equilibria and their different endpoints. It uses
analytical first-order responses at the original mixed-tenure economy. Dashed
lines are tangents to constant-price schedules, and arrows compare selected
states; no finite reform, nonlinear axis transformation or monotone path is
claimed. The shock direction is delta = epsilon/2 with intervention at date one.

One command builds this note and its figure without rebuilding slides or proofs:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 code/model/tools/build_simplified_olg_theory_note.py
```

`integrated_note_figure_checks.json` records the original-equation derivatives,
common inherited-state boundary, pre-announcement forecasts, stationary
intersections, full first-order arrays and source hashes.
`integrated_note_verification.json` records the lead's exact sign checks,
equation-preservation check, six bounded reviews, final two-pass compilation,
page-by-page visual inspection and final hashes. The lead replayed the existing
preference-shock checker; the stationary checker and original household
optimizations were also checked in the independent first pass. See
`docs/model/simplified_olg_overnight_work.md` for the concise synthesis.

The direct proof compensates the old in current consumption and offsets the
changed estate title with a bond, preserving all future real allocations.
The conditional fertility inequality applies at a date; the fully primitive
specialization additionally needs zero taxes and stationary or nonincreasing
expected house prices. The combined equilibrium theorem remains local, with
no stated shock-size radius, all-date fertility ordering or policy-welfare
claim. All author decisions remain recorded; no new instrument or convention
is adopted. Earlier work below is retained as supporting or historical evidence.

## September 6: first combined baseline-and-policy assessment

The latest work record establishes a local route in the original mixed-tenure
model: an initial preference decline and a permanent credit reform introduced
later from the same inherited state. Read the newest section of
`docs/model/simplified_olg_overnight_work.md` for assumptions, proof and limits.
This supersedes the earlier preliminary baseline-sign status. The integrated
note above now presents this result; the earlier slides and other PDFs below
remain as they were at this assessment.

The new stationary and baseline certificates can be replayed with:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 code/model/tools/verify_simplified_olg_stationary_endpoints.py
PYTHONDONTWRITEBYTECODE=1 python3 code/model/tools/verify_simplified_olg_fertility_decline.py
```

Receipts: `stationary_endpoint_checks.json` and `two_stage_transition_checks.json`.
The latter includes source hashes and lead check output. The policy-date extension
was checked analytically; no finite paired-policy simulation was run. No broad
primitive conditions, numerical reform radius, all-date fertility ordering or
policy-welfare conclusion is established. Proposed shocks are not author-adopted.

**Author clarification preceding the integrated note:** the desired illustration compares policy
against a baseline already adjusting after a fertility decline, with today
as an inherited point on that path. Explaining the initial decline is outside
the exercise; the suggested historical date remains tentative. The current
credit-only figure starts from a stationary baseline and does not yet show
that comparison. The checks below validate the existing artifact, not its
alignment with the author's clarified request. See the live decision ledger
and `docs/model/simplified_olg_overnight_work.md` for the two-path question,
unchanged steady-state notion, and preliminary local baseline-sign check.

## September 6: seven slides with equilibrium curves and the complete model

Tommaso asked to follow the earlier Raquel-style diagram more closely and
show every equilibrium element. The review PDF is
`output/pdf/simplified_olg_theory_slides.pdf`: seven main frames, identical to
pages 6–12 of the 75-page `output/pdf/september_14_presentation.pdf`.
The latter is also saved beside `latex/september_14_presentation.tex`.
Supporting theory is on pages 44–56; page 55 defines the plotted curves.

The sequence is environment/preferences, housing/finance, household choices,
competitive equilibrium, housing misallocation and compensation, conditional
fertility, and equilibrium transition. All four W/V problems are visible in
the main section, with their original notation and dated budgets. Equilibrium
includes tenure, rental entry, housing clearing, fiscal balance, cohort growth,
the old-state distribution, initial conditions, and external goods/bond trade.
The longer original household-problem frames remain in the appendix.

The two-panel transition again has curves, symbolic price/population/fertility
axes, replacement fertility, and initial A, impact I, and terminal A-prime.
Its visual reference is the August 17–21 advisor deck's demographic-adjustment
figure; retrieved records do not establish the exact Raquel meeting attribution.
The older quasi-linear preference-shock model is not substituted for the
accepted model. Instead, the original household equations give conditional
local housing-clearing curves in the same certified mixed-tenure credit example.
Initial and impact curves keep inherited old assets and the respective future
prices/rebate fixed. The long-run curve uses constant prices and repeated
lifetime choices. The right panel evaluates fertility along the same housing
curve. Only the marked states satisfy all equilibrium conditions.

The dotted impact curve conditions on inherited assets and expectations that
differ from their long-run values. This is a presentation choice, not a general
requirement for a third curve. The arrows compare the three
states; they do not assert monotone dynamics or a static one-state transition.
Axes are explicitly schematic. A monotone arcsinh vertical transformation
makes the small fertility impact visible; neither numerical distances nor
curvature are claimed. All three marked states and model parameters are
unchanged. The allocation figure still uses exact old-utility compensation;
the separate fixed-price fertility condition and all planner decisions remain.

Rebuild the figures, main deck, and extract from the single main TeX source:

```sh
python3 code/model/tools/build_simplified_olg_theory_slides.py
```

The driver uses `tmp/pdfs/september_14_theory_compact/` for build products;
`--figures-only` omits PDF compilation. `theory_slides_figure_checks.json`
records the curve definitions, derivatives, original-equation central checks,
matching equilibrium points, full analytical response, and compensation checks.
`theory_slides_verification.json` records final hashes, the page map, two-pass
builds, visual review, and preservation of both quantitative source blocks
and all pre-existing appendix frames. No overflow or undefined references;
the pre-existing main-deck appendix-bookmark warning remains.

No calibration, finite-horizon transition solve, new planner permission,
protected-manuscript edit, or full-note integration. All 18 decisions remain.
`september_slides_review.md` records this review and retains the earlier audit;
`september_slides_verification.json` concerns only the superseded long deck.

## September 6: transition with substantial renting and positive child costs

Start with the revised **four-page** `output/pdf/simplified_olg_paper_core.pdf`.
Its allocation and conditional-fertility pages are unchanged. Figure 2 now
uses the analytical mixed-tenure response. The **eight-page** supporting proof
is `output/pdf/simplified_olg_mixed_transition_proof.pdf`; both LaTeX sources
are under `latex/JMP_DS_suggestions/`. The earlier limiting figure, ten-page
assessment and all-owner proof are preserved below.

The new exact six-variable representation retains original dated budgets,
positive child costs, a positive rebated property tax, finite logistic tastes,
endogenous tenure and the actual initial old's assets and housing titles.
An explicit economy has `11/21` owners, child goods cost `3/20`, and property
tax `467/9250`. The same stationary bundles support a family of taste scales.
Within each credit reform both taste parameters remain fixed.

- For every positive taste scale: four stable/two unstable roots, a positive
  stationary population derivative, and negative stationary conditional
  lifetime-value derivatives in both tenures. Exact polynomial coefficient
  signs prove these statements throughout the positive half-line.
- For every scale in `[1,4]`: the actual initial-old boundary is transverse
  and initial fertility rises. Rational intervals cover the entire interval,
  with 913 adjacent cells; this is not just a list of parameter-point solves.
  The sequence implicit-function argument gives locally unique nonlinear
  converging paths after small permanent credit relaxations.
- Near scale `1`: a nonzero dominant complex mode and a uniform nonlinear
  remainder bound prove that small finite reforms eventually approach
  replacement fertility, and final population, from both sides. The initial
  fertility increase and final population gain do not imply monotonicity.

Small primitive and compact-support entrant-distribution changes preserve the
strict conditions. The reform radius and a multidimensional parameter box are
not quantified. This remains an illustrative local result, with fixed entrant
endowments and world bond price. Population growth is distinct from the
compensated allocation welfare gain and from a social population ranking.
No planner permission, author convention, or quantitative-model change is made.

Evidence:

- `mixed_transition_proof.md`: full derivation, general local theorem,
  positive-tax inversion lemma, exact construction, finite-reform oscillation
  proof, and the wider family result.
- `mixed_transition_map_review.md`, `mixed_transition_certificate_review.md`,
  `mixed_transition_family_review.md`: three sequential bounded read-only
  reviews of distinct new steps. The second requested an explicit uniform
  tail lemma; the third confirms that the added lemma closes that issue.
  The last review audited source and receipt but could not replay them in its
  model environment because `sympy` was absent. The lead executed the complete
  checks in the working Python environment.
- `mixed_transition_certificate.json`: exact rational matrices, root bounds,
  boundary determinant, response signs and complex-mode signal; also records
  the analytical Figure 2 values.
- `mixed_transition_family_certificate.json`: all-positive-scale Routh,
  discriminant and stationary-value proofs; complete `[1,4]` interval coverage;
  original-equation comparisons at four declared scales. Completed in about
  91 seconds with a 600-second/12,000-subdivision fail limit.
- `mixed_transition_smoke.json`, `mixed_transition_checks.json`: exact original
  branch checks, four short original-equation paths, 24 household optimizations,
  central derivatives and 24/40-date horizon comparison. Maximum equilibrium
  residual is below `2.9e-15`, budget error below `1.6e-15`, and initial/final
  derivative discrepancies below `3e-9`. Finite paths check arithmetic; they
  do not supply the convergence proof or the plotted curves.

One new driver reproduces the work; the old limiting driver is unchanged:

```sh
python3 code/model/tools/verify_simplified_olg_mixed_transition.py --smoke
python3 code/model/tools/verify_simplified_olg_mixed_transition.py
python3 code/model/tools/verify_simplified_olg_mixed_transition.py --certificate-only --figure
python3 code/model/tools/verify_simplified_olg_mixed_transition.py --family-only
```

The verified interpreter is `/Library/Developer/CommandLineTools/usr/bin/python3`,
with NumPy 2.0.2, SciPy 1.13.1, SymPy 1.14.0 and Matplotlib 3.9.4. Use that
interpreter for these commands if the active model virtual environment lacks
SymPy. No package or environment changes were made for this proof.

All four receipts match the final new source hash. The supplemental figure
files are `mixed_transition_analytical_figure.pdf/png`; the compact note uses
this response as Figure 2, with the allocation figure unchanged.

## September 5 late-evening scope and prose pass

Start with `output/pdf/simplified_olg_paper_core.pdf`: **four pages**, comprising
three pages of proposed results after the existing household setup and one page
of economic assessment. Both existing figures are retained. Source:
`latex/JMP_DS_suggestions/simplified_olg_paper_core.tex`. The ten-page assessment
below remains the fuller reading note; neither document is an adopted paper edit.

The new section 11 of `local_transition_proof.md` explains the limiting
stationary credit result for heterogeneous entrants. Price rises in proportion
to common borrowing capacity, each young type's housing and fertility remain
unchanged, and every old type uses less housing. Stationary population rises,
but conditional owner lifetime utility falls for every same-type entrant.
This holds at the all-owner, zero-child-cost/tax limit with the original
uniformly strict branches and fixed income, wealth distribution and world bond
price. It is separate from transitional-cohort welfare, variable-population
social welfare, and the compensated allocation result.

- `local_transition_scope_review.md`: one bounded read-only second opinion,
  independently deriving the general welfare formula and assessing the economic
  assumptions. The lead qualifies its estate wording: the original estate
  inequality is not a minimum-retention rule and does not prohibit a sale.
- `local_transition_welfare_checks.json`: exact symbolic utility and housing
  identities, original budgets and constraints for five fixed types at three
  financed shares, ten original household optimizations, and utility derivatives
  checked to `8.3e-9`. Regenerate with the existing driver:

```sh
python3 code/model/tools/verify_simplified_olg_local_transition.py --welfare-only
```

Recommendation, pending author discussion: retain the allocation and conditional
fertility arguments as the core; treat the second figure as a separate local
demographic extension, with the price and welfare interpretation explicit.
Near-all ownership, small child goods cost, fixed entrant endowments and
unchanged household constraints limit its economic coverage. The simple
stability inequality does not make all these restrictions innocuous. No new
planner power, ownership preference, population convention, or model is chosen.

## September 5 continuation: analytical transition and primitive conditions

The longer reading note is **ten pages**, at
`output/pdf/simplified_olg_simple_assessment.pdf`; its LaTeX source remains
`latex/JMP_DS_suggestions/simplified_olg_simple_assessment.tex`. Pages 1–4 keep
the allocation and conditional fertility arguments. Page 5 replaces the
prescribed demographic curve with the analytical equilibrium response.
Appendix B gives primitive stability and initial-fertility conditions;
Appendix C gives a stationary population condition with positive child costs.
No preference or planner-instrument decision is adopted.

The all-owner limit with zero child goods cost and property tax has an exact
aggregate recurrence, including with heterogeneous entrants on uniformly
strict household branches. The condition `ell + (3+q) C > 4D`, with each
coefficient defined from primitives in the proof, gives exactly two stable
roots and one unstable root. The predetermined initial cohort and actual old
assets select a local converging equilibrium. The exact initial-fertility test
uses a primitive polynomial evaluation; a simpler sufficient inequality is
also available. These extend locally to positive renter mass, child costs and
property tax through the full dated equilibrium operator.

The original model's initial-old accounting, endogenous tenure, extra price
leads, normalized type distribution and all private constraints are retained.
The positive-parameter neighborhood and reform must be small; their numerical
size and a global transition theorem remain unproved. Initial fertility and
terminal population rise under the stated conditions. All-date finite-reform
fertility and monotone population are proved only for the plotted limiting
example. An admissible second economy instead has lower initial fertility and
a larger final population. The uncompensated credit reform lowers stationary
household welfare in the plotted limiting example; it is not the compensated
Pareto allocation proof.

Supporting files:

- `local_transition_proof.md`: full derivation, general primitive conditions,
  exact aggregation, uniform infinite-tail proof, welfare distinction, and
  counterexamples. Sections 10.4 and 10.2 contain the broad stability condition
  and exact initial-fertility test.
- `local_transition_independent_review.md`: first review of the explicit
  recurrence, household choices and original initial-old boundary; identifies
  the subsequently closed mixed-extension and uniform-tail gaps.
- `local_transition_closure_review.md`: independent check of those two newly
  completed proof steps and the negative welfare derivative.
- `local_transition_heterogeneity_review.md`: independent aggregation and
  positive-child-cost stationary-condition review; required probability,
  inherited-state and uniform-operator qualifications are incorporated.
- `local_transition_general_review.md`: independent general convergence review
  under `w>d`, the exact/sufficient fertility tests and negative initial-sign
  example. This report predates the sharper condition in section 10.4.
- `local_transition_sharp_stability_review.md`: final narrow independent check
  of that sharper root condition, including zero/repeated stable roots and
  carryover of the original initial-old boundary and primitive fertility tests.
  Its requested standalone operator/decay-weight qualification is incorporated.
- `local_transition_checks.json`: eight original-equation cases, original
  household optimizations, a central difference versus the analytical path,
  and 24/40-date horizon comparison. Maximum equilibrium residual is below
  `2.3e-14`; derivative error is `1.48e-9`. Finite paths are arithmetic support,
  not the proof of infinite convergence or an admissible-neighborhood bound.
- `local_transition_stationary_checks.json`: exact symbolic derivative and
  three original stationary cases with both population signs at positive cost.
- `local_transition_heterogeneity_checks.json`: zero Hessians of eight owner
  quantity/asset maps, four-corner constraints and 24 original household
  optimizations, plus exact aggregate dated-demand and inherited-old checks.
  Its dated prices are diagnostic inputs, not equilibrium paths.
- `local_transition_general_checks.json`: six symbolic polynomial/transform
  identities, three initial-boundary identities, 270 declared coefficient
  cases (including unstable cases), and an original-equation counterexample
  with initial fertility derivative `-0.0597632` and final population derivative
  `0.0952381` per proportional credit-cap increase. This grid is an arithmetic
  check, not a calibration or an economically admissible-parameter map.
- `local_transition_smoke.json`: the exact zero/tiny-reform loop smoke check.
- `local_transition_analytical_figure.pdf` / `.png`: Figure 2, obtained directly
  from the analytical formula, never from the finite-path solver.

All verification modes take seconds locally and record the current driver hash.
From the repository root:

```sh
python3 code/model/tools/verify_simplified_olg_local_transition.py --smoke
python3 code/model/tools/verify_simplified_olg_local_transition.py
python3 code/model/tools/verify_simplified_olg_local_transition.py --stationary-only
python3 code/model/tools/verify_simplified_olg_local_transition.py --heterogeneity-only
python3 code/model/tools/verify_simplified_olg_local_transition.py --general-only
python3 code/model/tools/verify_simplified_olg_local_transition.py --figure
```

The first unscaled optimizer check at the deliberately tiny rental cap was
ill-conditioned and did not report success. Scaling the optimizer's choice
units resolved it; the original objective and constraints were retained. All
final modes pass with independent household optimization and strict private
conditions. No failed check was accepted as proof.

## September 5: earlier simple allocation assessment

The initial reading note was the eight-page
`output/pdf/simplified_olg_simple_assessment.pdf`, with source
`latex/JMP_DS_suggestions/simplified_olg_simple_assessment.tex`. It restores
the author's priority: first establish housing misallocation simply; then
assess fertility, population and constrained-finance extensions. It does not
adopt new preferences or public collection powers.

- `simple_assessment_checks.json`: four explicit mixed-tenure equilibria,
  nine independent original household optimizations, 36 finite compensated
  reallocations, owner-to-renter replication, conditional fertility signs,
  demographic accounting and stationary stock scaling. Driver/specification
  hashes are recorded. The family varies income and ownership-taste location
  with the tax rate; it is not a tax comparative static.
- `simple_assessment_direct_review.md`: independent direct-proof and
  constraint-role review. Two transcription corrections are disclosed at the
  top. This review preceded the new equilibrium constructions.
- `simple_assessment_new_claims_review.md`: separate independent derivation
  of replication, the equilibrium family, the slack-rental equilibrium, and
  stationary scaling. Its scope and incorporated clarifications are explicit.
- `simple_allocation_figure.pdf` / `.png`: exact compensated reallocation for
  the eligible pair in the constructed equilibrium; the crossing is the best
  allocation on this comparison, not a full economy-wide first best.
- `simple_population_figure.pdf` / `.png`: a prescribed fertility sequence
  and the resulting demographic identity. This is not a solved equilibrium
  transition or a prediction of adjustment speed. This earlier figure is retained
  but the latest note uses the analytical Figure 2 described above.

Regenerate this bounded packet in a few seconds, without calibration or
equilibrium search:

```sh
python3 code/model/tools/verify_simplified_olg_simple_assessment.py --plot
```

The full derivations and assessment are also recorded near the top of
`docs/model/simplified_olg_overnight_work.md`. A complete equilibrium example
now establishes that a positive group of eligible owner pairs can exist.
A second example establishes that the binding rental cap is not necessary
with the existing ownership taste. The continuation above supplies local
analytical credit-reform results; global transitions and anonymous public-loan
implementation remain unproved.

## September 5: external discussion assessment

`external_discussion_review.md` preserves the independent agent's response
returned by the author. It recommends public down-payment finance and revising
the interpretation of the existing conditional proof. It is a discussion input,
not acceptance of planner powers or verification of its proposed anonymous
pure-loan implementation. The lead's qualifications are in the live decision
ledger; the existing mathematical sources and receipts remain unchanged.

## September 5: focused constrained-efficiency result (latest)

The author asked to settle constrained inefficiency before moving to other
issues. The two existing Astra reviewers and the lead completed a new pass.
The current four-page synthesis is
`output/pdf/simplified_olg_constrained_efficiency.pdf`, from
`latex/JMP_DS_suggestions/simplified_olg_constrained_efficiency.tex`.
This result supersedes the earlier statement below that the complete-path
constrained comparison is entirely unproved.

- `constrained_full_path_review.md` proves a complete infinite-path Pareto
  improvement with committed household-specific young and old transfers,
  enforceable future taxes, and an explicit passive initial title owner with
  outside bond access and participation in the fiscal account. Young housing
  rises at every date. Initial young owners can gain strictly. Every market
  clears and all initial claims are included. It also proves a local obstruction
  to one-time gifts around the verified stationary regime, not a global no-go.
- `constrained_instruments_review.md` gives a separate finite-support fallback
  using a current-income advance and two targeted housing incentives. It is an
  expanded-instrument result, not a claim of globally minimal instruments.
- `planner_benchmark_checks.json` now also records three finite committed
  reforms, a fourth case with a strict initial-young gain, original-budget and
  full-resource identities, the infinite contraction bound, 12 independent
  household optimizations, full tenure-deviation checks at fixed individual
  fertility, and the one-time-gift price roots.
- `committed_cash_path.csv` contains dated receipts for the analytical path.
  `committed_cash_diagnostics.png` is a supplemental six-panel internal check;
  the agreed theory and earlier numerical figures are unchanged.

Regenerate the complete receipts and supplemental diagnostic with SciPy and
Matplotlib installed (about five seconds locally):

```sh
python3 code/model/tools/verify_simplified_olg_planner_benchmarks.py --original-optimizers --plot
```

The positive theorem is conditional on the stated transfer and ownership
permissions, stationary group structure, strict private branches and tenure
margin, and a contraction condition. Author acceptance of those permissions
remains OPEN in the decision ledger. In particular, public commitment provides
financing across ages; an unchanged mortgage share does not remove this extra
economic power. No new fertility externality is used. The author manuscript
and quantitative model were not edited.

## September 5 planner benchmarks

The author requested two Astra agents at maximum reasoning, one for direct
allocation and one for cash transfers followed by markets. The lead reviewed
both arguments and performed independent algebra and finite-allocation checks.
The five-page synthesis is `output/pdf/simplified_olg_planner_benchmarks.pdf`;
its source is `latex/JMP_DS_suggestions/simplified_olg_planner_benchmarks.tex`.

- `first_best_review.md`: exact current-goods compensation, original title and
  estate accounting, conditional marginal-value diagram, and global choices.
- `constrained_efficiency_review.md`: exact cash-transfer envelopes, a local
  obstruction and four-household-group improvement. The latter moves housing
  toward the old and omits future market clearing and intermediary-owner welfare.
  Neither example establishes constrained inefficiency of the complete OLG model.
- `planner_benchmark_checks.json`: independent receipts, including source hashes,
  27 direct perturbations, three exact conditional cash cases, and the obstruction.

Reproduce these small checks without equilibrium or calibration solves:

```sh
python3 code/model/tools/verify_simplified_olg_planner_benchmarks.py
```

The live decisions and remaining ownership/transfer/estate questions are in
`docs/model/ACTIVE_DECISION_LEDGER.md`. The earlier amendment and its evidence
remain below as a separate completed pass.

## Earlier amendment checks

This is the supporting output folder for the September 4–5 theory development.
The live phase/checkpoint record is `docs/model/simplified_olg_overnight_work.md`.
The proposal is `latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex`;
its readable version is `output/pdf/simplified_olg_amendment_proposal.pdf`.

Regenerate the bounded analytical checks from the repository root:

```sh
python3 code/model/tools/verify_simplified_olg_amendments.py
```

- `verification_summary.json`: check counts, errors, source/code hashes,
  constructed original-owner bundles, and common-cap correction.
- `fertility_checks.csv`: direct scalar household solutions and finite differences.
- `welfare_checks.csv`: exact compensation, uniform premium, and ordinary
  market-price comparisons; dated financing and utility changes.
- `analytical_checks.png` / `.pdf`: the stable two-panel diagnostic for these
  analytical examples. No equilibrium or calibrated result is shown.
- `independent_review.md`: completed bounded read-only reviewer report. The
  lead adopted its tax-reserve and seller-receipt clarity corrections. Its
  original line references refer to the reviewed checkpoint, not later edits.
- `independent_review.log`: transient wrapper log; not a canonical claim source.

The primitive sufficient condition and allocation proof retain their explicit
stationary/fixed-price/branch qualifications. These checks do not establish a
GE policy sign or broad transition theorem. No production calibration is run.

Regenerate the small illustrative equilibrium checks and the full six-panel
packet with:

```sh
python3 code/model/tools/verify_simplified_olg_amendments.py --transitions
```

The predeclared experiment is four financed shares (80%, 81%, 85%, 90%), with
identical initial cohorts, fixed fertility preferences, zero gains tax, and
property tax 5%. A horizon-40 comparison, second solution method, original
household optimization, and endpoint derivative checks are included. Expected
local runtime is about three minutes, with a ten-minute stop. This is a theory
example, not a calibration or production quantitative policy run.

To rebuild the compact report figure and check population accounting from
existing paths, without solving another equilibrium:

```sh
python3 code/model/tools/verify_simplified_olg_amendments.py --saved-transitions
```

- `transition_verification.json`: all policy summaries, parameters, original
  optimizer discrepancy, longer-horizon and second-solver differences.
- `transition_phi80.csv`, `transition_phi81.csv`, `transition_phi85.csv`,
  `transition_phi90.csv`: full path quantities and original-equation checks.
- `transition_phi85_h40.csv`: longer-horizon comparison.
- `transition_stability_checks.json`: initial and final endpoint linearizations,
  root counts, inherited-state projection, and derivative step comparisons.
- `original_household_optimizer_checks.csv`: 24 direct original-budget
  constrained optimizations, separately checking the reduced choices.
- `population_and_composition_checks.json`: cohort-product identity and exact
  initial-tenure-share decomposition of the fertility change.
- `credit_policy_transitions.pdf` / `.png`: full six-panel diagnostic packet.
- `credit_policy_summary.pdf` / `.png`: compact fertility/population figure for
  the consolidated proposal, with the full packet retained above.
- `transition_independent_review.md`: bounded independent mathematical review
  of the continuation theorem and finite-history argument.
- `final_independent_review.md`: final independent synthesis review; the lead's
  resolution is in the work record. Line numbers refer to the reviewed version.
- `delivery_manifest.json`: source/output hashes and the final verification
  record. Numerical solver-function hashes are distinct from presentation edits.

The 90% reform has a slightly negative initial fertility response despite a
larger final household population. The owner borrowing limit is slack in that
case. These findings are retained; the example does not establish a positive
fertility response at every date or verify interval-wide continuation hypotheses.

Build the note from the repository root, directing all intermediate products
outside the protected author manuscript. Run this twice:

```sh
/Library/TeX/texbin/pdflatex -interaction=nonstopmode -halt-on-error -output-directory=tmp/pdfs/simplified_olg_amendments latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex
```

The final PDF is copied to `output/pdf/simplified_olg_amendment_proposal.pdf`
after rendering and visual inspection. Wrapper logs, progress JSON, runtime
records, and temporary build products are not canonical evidence.
