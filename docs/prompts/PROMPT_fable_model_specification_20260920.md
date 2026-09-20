# Fable review: choose a provisional model before serious recalibration

Prepared for Tommaso's Claude Max session, September 20, 2026.
Project: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.

## The decision Tommaso wants help making

We have been considering several changes to a lifecycle housing and fertility
model while also trying to improve its calibration. Tommaso now wants to settle
the model specification first, then undertake a serious joint calibration of a
provisionally frozen specification. The immediate task is to recommend that
specification and identify the few author decisions needed to settle it.

Treat this as an independent economics review. Evaluate the case for retaining
each current feature as well as the case for changing it. Earlier Fable and
ChatGPT recommendations are inputs to reassessment, not author decisions or
established conclusions. Recommend a coherent, parsimonious baseline rather
than a collection of every available extension.

Tommaso prefers a conventional macro earnings specification and has asked to
try persistent plus transitory earnings without a separate permanent type.
That is a preference to assess explicitly, not an instruction to select whichever
process produces the desired mechanism. Identify the exact referenced
Boar/Gorea/Midrigan paper and its earnings equation before claiming equivalence.

## Scope and source discipline

One bounded review pass, approximately 20 minutes if the session supports a
time limit. Return a complete first assessment and clearly identify unresolved
items. Use existing evidence and targeted source inspection. Do not launch
calibrations, simulations, cluster jobs, nested review campaigns, or code edits.
Return the assessment in this conversation. Do not edit author LaTeX, the active
decision ledger, model defaults, targets, or weights.

If this session has repository access, follow its startup instructions, then
read the following in order. If it lacks access, say which sources need to be
supplied; do not pretend to have inspected local paths.

1. `CALIBRATION_STATUS.md`, especially the September 20 reviewed completion.
2. `docs/prompts/HANDOFF_chatgpt_restart_20260919.md`, particularly the
   sandbox distinction, switches, and pending author decisions.
3. `docs/model/structural_model_review_consolidated_20260915.md`, including
   the follow-up log. Read individual topic sections as needed.
4. `output/model/native_financing_diagnostic_20260919/README.md`, the
   September 20 reviewed result at the top.
5. `output/model/native_financing_diagnostic_20260919/overnight/final_search/readout.md`
   and `overnight/final_review.json` under the same experiment folder.
6. `output/model/native_financing_diagnostic_20260919/overnight/final_mechanisms/comparison.md`,
   its `comparison.csv`, and `cohort_accounting_check.json`.
7. `docs/model/ACTIVE_DECISION_LEDGER.md`, only the relevant open entries on
   earnings, geographic scope, and calibration versus mechanism. It contains
   older and unrelated theory entries; it does not override current status.

For exact specification claims, distinguish active development code, the
September-14 frozen native source, and Claude's sandbox. The frozen source and
checkpoint are pinned in the experiment's receipts. The old sandbox is not an
interchangeable numerical benchmark. Default-off implementation and passed
tests do not establish author adoption. Some reviewed source files have
unrelated uncommitted edits; reading them does not authorize modifying them.

Use the primary economics paper for a literature-dependent recommendation,
with its relevant equation or passage identified. Existing literature memos
are pointers whose claims may need checking. State uncertainty if the source
cannot be verified. Do not add a broad literature survey.

## What the overnight work establishes, and its limits

The new-income candidate combines a deterministic age profile, five persistent
AR(1) states and three iid transitory states, without a permanent type.
The finite joint search evaluated 96 proposals, with 89 valid and seven
rejected. The selected loss is 353.6588729140903, compared with 502.74561411262744
for the prior local candidate and 179.2984242480252 for the original benchmark.
The 13 target definitions, values and weights are identical in that comparison.
Two selected native repetitions reproduce fit, loss and price exactly.

This is not a converged or adopted calibration. Selected ownership at ages
30–55 is 0.48427 against 0.64833, and capped mean rooms 6.67698 against 5.56110.
Annual discounting and the first-child housing requirement are at the actual
upper bounds, 0.99 and 2.3. Read the complete 13-target and 17-parameter table
before judging the fit. The raw scorer's generic beta upper bound 0.9995 is not
the search bound. The normalized child-preference level is 0.1922702943;
0.2394995040 in the search summary is its initial normalization seed.

All three 48-cell financing/rental grids completed: original earnings,
new earnings before structural refitting, and the selected new-income refit.
Within each family, preferences, prices and the initial population are fixed.
Across families these objects, including native entrant wealth, can differ.
Financed mortgage shares are 0.8/0.9/0.95/1, unsecured limits are
0/0.25/1/5 times four-year after-tax earnings (zero in retirement), and rental
room caps are 6/8/10. Thus 0.25 times period earnings is one annual earnings
amount; the largest dose is twenty annual earnings amounts, not five years.

Mortgage relaxation changes both deposit and collateral access. At the
six-room rental cap, total birth flows rise about 0.690%, 0.174% and 0.230%
in the three families; explicit lifetime cohort births rise 0.234%, 0.085%
and 0.081%. Larger rental support substantially reduces these increments.
In the refit, unsecured credit equal to one annual earnings amount raises
snapshot birth flow 2.456% but lowers lifetime cohort births 0.769%. At four
annual earnings the changes are +2.380% and -2.415%. These cohort arms have
exactly identical entry populations and their age-by-age totals reproduce
the reported lifetime sums: early births rise, later births and lifetime first
births decline. The economic or numerical cause remains unresolved. The equal
outcomes at credit multipliers 1 and 5 also need explanation.

The separate 2.1 stationary normalization is not the explicit lifetime cohort
birth total, a population forecast, or a female period TFR. First-birth age in
the cohort and the calibration observer also use different age conventions;
do not compare their raw levels as if they were the same target.

These are fixed-price diagnostics, not a literal frictionless benchmark,
identified mediation, or general-equilibrium policy results. Poor rooms,
wealth, ownership and fertility lifecycle fit can put prospective parents in
the wrong states for mechanism evaluation. Therefore the latest discussion
explicitly withdrew a strong verdict against the mortgage mechanism based on
these experiments alone. Conversely, poor fit does not establish that a future
calibration will restore a large mortgage effect. Assess what is actually known.

## Specification choices to resolve

Use this as a grounded starting inventory, grouping closely related choices.
Correct it where the source establishes a different implementation or status.

| Topic | Open choice and source |
|---|---|
| Earnings | Permanent heterogeneity versus persistent plus iid transitory risk; a coherent income definition, age profile, time aggregation and tax treatment. Review E1; latest native candidate and full fit above. |
| Child benefits and costs | Linear versus concave child benefit; fixed first-birth utility cost versus separately disciplined goods/time/earnings costs. Review P1/P3 and the child-cost literature memo. Distinguish a primitive change from a calibration adjustment. |
| Dependent children | Memoryless departure, parent-age-dependent departure with newborn exemption, or an explicit child-stage state. Address surviving dependents when parents die consistently. Review F3/F4. |
| Mortgage and unsecured borrowing | Present rolling collateral/deposit contract versus amortization and origination-only borrowing; economically defensible unsecured credit. Review H1/H2. Separate internal accounting requirements from optional realism. |
| Housing services and tenure access | Family space requirement, ownership premium, hard rental-size cap, and a size-dependent rental wedge. Review H3/H6 and the housing-floor discussion. Which objects need independent data and which are substitutes? |
| Tenure choice shocks | Retain and identify the small tenure shock, or use deterministic tenure choice. Review F2. Distinguish economic heterogeneity from numerical smoothing and explain the equilibrium implications. |
| Bequests and estates | Net-of-selling-cost estate valuation, recipients, and any proposed late-life downsizing mechanism. Review P4/P5. Keep valuation consistency distinct from adopting a richer recipient or health model. |
| Household and population accounting | Reproductive-household units, persons, dependents and entrants; outside inflow versus contraction. Review D1/D2/F4. Identify what must be settled for initial calibration and what specifically blocks transitions or welfare. |
| Geographic scope | National versus selected-metro empirical population, with compatible housing, earnings, fertility and demographic evidence. Review D4 and the geographic entry in the active ledger. Existing mixed samples cannot be relabeled. |
| Housing supply over time | Reversible static supply versus an inherited stock with construction/depreciation. Review H4. Distinguish a stationary baseline choice from the requirements of a transition experiment. |
| Additional fertility or household margins | Unintended births, household formation, widowhood or health states appear in the earlier review. Explicitly justify inclusion now or defer them; do not let an optional extension silently become a prerequisite. |

Existing switches are summarized in handoff section 4 and
`code/model/sandbox/specs/`. Their numerical constants are diagnostic choices,
not automatically defensible baseline parameters. The pending recommendations
in handoff section 8 remain recommendations unless Tommaso adopted them.

## Requested response

Keep the main response around 1,500–2,000 words plus one compact decision table.

1. Recommend one coherent provisional baseline, stating its economic logic and
   its main tradeoffs. Provide one conditional fallback only if a genuinely
   decisive unresolved choice prevents a single recommendation.
2. Give a decision table: topic; retain/change/defer recommendation; reason;
   supporting evidence and remaining uncertainty; state-space/computational
   cost; identifying data or external restriction; blocks specification freeze
   or can wait. Identify precisely which recommendations require author choice.
3. Name the three highest-value author decisions to discuss first and the
   smallest additional piece of evidence needed for each, if any. A ranked
   shortlist is more useful than reopening every past issue equally.
4. Describe the subsequent calibration contract: which moments/profiles must
   be fitted or validated, what disciplines each changed parameter block, and
   which measurement definitions must be reconciled. Do not add free parameters
   without identifying variation or external restrictions, or silently drop or
   reweight difficult targets. An aggregate loss alone is not sufficient.
5. Explain which earlier mechanism conclusions survive this evidence and which
   should remain conditional on a credible baseline. Distinguish source-code
   facts, mathematical implications, outcomes at a particular calibration and
   economic conjectures. In particular, reassess the old claims about children
   being net costs, unsecured credit increasing fertility, and a frictionless
   benchmark rather than inheriting them unchanged.
6. End with the proposed freeze point and the ordered work after it: only the
   necessary implementation and measurement changes, then joint calibration,
   lifecycle validation, and finally comparable mechanism/policy tests.

Use plain economics prose and LaTeX for equations. Use “number of children,”
“children ever born,” or “fertility” rather than “parity.” Cite relevant local
sections and verified primary sources. No model modification is being adopted
by requesting this assessment.
