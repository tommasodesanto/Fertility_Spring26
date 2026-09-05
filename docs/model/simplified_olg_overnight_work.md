# Simplified OLG: development plan and overnight record

Authorized by Tommaso on 2026-09-04: “yeah, let's do this.. so you can work on
this overnight while i sleep hopefulyl.” This begins development after the
discussion pass. It does not convert provisional economic choices into final
author decisions. The issue tree remains in `ACTIVE_DECISION_LEDGER.md`.

## Scope, deadline, and deliverables

**Latest author instruction:** work continuously until finished; do not wait for
hourly runs. Tommaso also requested three checks of every substantive result
and a simpler writing style, guided by the recent conversations. The hourly
automation is PAUSED. Continue in this task, saving progress as each part is
completed. The earlier 08:00 checkpoint is a reporting aim, not permission to
stop with work that can still be completed. Stop when the agreed development,
verification, and revised note are finished; state precisely any mathematical
claim that remains unproved or requires an author choice.

Use three distinct checks: (1) derive the result, (2) compare it with the
original household and equilibrium equations, and (3) test it independently,
using an adversarial proof review or an independently constructed calculation
as appropriate. Repeating a solver three times is not a substitute for these
checks. Review the finished prose and PDF as well.

Produce one consolidated, readable LaTeX/PDF amendment proposal, supporting
reproducible checks, and a short morning readout in this file. For each claim,
distinguish **proved under stated assumptions**, **numerically illustrated**,
**candidate**, **not established**, and **requires author decision**. A failed
conjecture is a valid research outcome; do not replace a missing theorem with
assertions. No particular theorem is guaranteed by the deadline.

Proposed manuscript: `latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex`.
Readable output: `output/pdf/simplified_olg_amendment_proposal.pdf`.
Supporting evidence: `output/model/simplified_olg_amendments/`.
Use existing theory drivers where suitable; one focused verification driver
may own the new analytical checks. Preserve the independent review as historical
evidence. No files in `latex/JMP_DS_draft/` may be written, including build files.

## Accepted directions and retained decisions

- W0/W1: establish an allocation gain at fixed fertility and cohort masses;
  borrowing limits may be relaxed. Compare housing willingness to pay in
  consumption units, and respect real resources, creditors, and estates.
- W3/W4 are **provisional**: household-specific compensation and a nonnegative
  real reallocation cost. Show the actual transactions and a simpler-transfer
  comparison. Return these concrete choices to the author before presenting
  the welfare benchmark as final. Keep marginal costs and fixed moving costs
  distinct. Seller compensation is a transfer, not resource destruction.
- I0/Q0/W2: capital-gains taxation is deferred. Develop the core at zero gains
  tax. Property tax is separate and is not silently removed. A zero-property-tax
  specialization can be informative if explicitly labeled as such.
- H0/H1: one young completed-fertility choice and one old period. Seek an
  explicit sufficient condition, preferably with no prices or multipliers;
  distinguish household partial effects from equilibrium policy effects.
- D0/D1/D2: compare policy and baseline transitions from the same initial
  state, then total population at their endpoints. Seek the broadest valid
  transition result. Do not assume small reforms in advance or equate finite
  terminal closure with an infinite-horizon theorem. Fixed housing stock is
  the benchmark; assess a simple stationary supply extension.
- U0: repair inconsistent units by an explicit change of variables; keep
  adult households, children, and resident persons distinct. Do not change
  the quantitative model's demographic contract.
- P0/P1: presentation is secondary. Preserve all 18 issue IDs and 14 review
  repairs, including parked alternatives and the conditions for revisiting them.

## Work sequence and acceptance checks

| Phase | Economic task | Required evidence | Status |
|---|---|---|---|
| 1 | Construct the compensated old-to-young reallocation at fixed fertility. | Dated goods, title, debt, sale, and estate ledger; admissibility/slackness conditions; costs and welfare coverage; simpler-transfer alternative. | CHECKED local proposition; compensation/cost remain provisional |
| 2 | Derive and audit the fertility threshold. | Exact sign calculation; a condition without multipliers; examine a genuinely primitive condition; negative-response example; capped-old-renter and common-cap correction. | CHECKED local threshold and primitive stationary specialization; aggregate policy sign open |
| 3 | Establish the transition result and its limits. | Separate existence, convergence, uniqueness, and policy signs; attempt broader continuation under explicit bounds; identify failures or missing hypotheses; check the zero-shock system when using a local theorem. | CHECKED conditional finite-reform theorem; illustrated credit transitions with property tax retained. Interval-wide applicability remains unproved. |
| 4 | Reconcile population units and stationary supply. | Exact rescaling, transition population identity, common initial state, terminal household/person definitions; supplier-income/financing qualifications. | CHECKED identities, original laws, and independent rescaling calculations; literal-person convention still requires author review |
| 5 | Consolidate amendments and deliver the PDF. | Reconcile proposal against section, appendix, builder, and claim ledger; relevant small checks; compile twice, inspect rendered pages; preserve frozen evidence and unrelated edits. | COMPLETE: consolidated 12-page proposal, all 14 repairs reconciled, three checks recorded, final review passed, all PDF pages inspected. |

Prioritize phases 1–2, then 3–4, then presentation. If a broad transition
theorem cannot be justified, deliver the strongest valid result and the precise
remaining obstacle. Do not silently weaken the requested objective.

## Execution and file ownership

The lead owns economic judgment, derivations, synthesis, and final verification.
Use at most one independent worker at a time for a bounded adversarial check
through the repository wrapper. Default reviewer limit: 20 minutes. Inspect its
actual output; do not accept its conclusion without checking the argument.
No duplicate worker if an earlier one is still active. No silent retries after
timeout. Record any worker PID/session, output, scope, and completion below.

Only small illustrative theory computations are in scope. No calibration,
production policy runs, target changes, or quantitative sweeps are authorized
by this theory request. Before any numerical experiment, record solve count,
estimated time, and stopping rule; smoke-test any loop before expansion.
Keep individual local experiments under 10 minutes and total local numerical
experiments under 30 minutes. Anything requiring a longer run needs a new
explicit plan under the project's cluster procedure, not an unattended launch.

The checkout contains substantial unrelated work, including another overnight
quantitative task. Never stage all changes. Prefer the new suggestion/check
files and targeted discussion-map edits; inspect diffs before changing existing
theory sources. Source reconciliation can document a proposed replacement
without promoting provisional assumptions into the active manuscript. Keep
`CALIBRATION_STATUS.md` and quantitative code unchanged. Commit and push only
coherent, verified changes owned by this task. Do not force push.

## Checkpoints and morning readout

### Initial checkpoint — 2026-09-04, approximately 22:10 EDT

- Startup context and decision tree loaded; active theory section inspected.
- Work plan recorded. Hourly same-task continuation is active as `simplified-olg-overnight-theory`;
  it must pause at completion or the stated morning deadline.
- New fertility lead: on the strictly binding down-payment branch, the
  stationary housing expenditure at the cap can be independent of the house
  price. This may yield a primitive sufficient condition when lifetime
  transfers are exogenous (in particular in a zero-property-tax specialization).
  It has not yet passed an independent check and is not a transition theorem.
- Next: write and verify the candidate, and construct the dated compensation
  ledger while one bounded reviewer checks the fertility argument.

Append substantive checkpoints here; update the phase table rather than
creating a separate status file each hour. The morning readout should lead with
what is established, then what failed or remains conditional, then the small
set of concrete choices for Tommaso. Include source/PDF/check links.

### First analytical verification budget

- Worker: `reviewer_strong`, read-only, wrapper session 39561, 20-minute limit;
  output `output/model/simplified_olg_amendments/independent_review.md`.
  Scope: the two new welfare/fertility arguments only; no transition audit.
- Local verification: 6,075 candidate reduced-household configurations, with
  only binding cases solved and three scalar roots per retained case; smoke
  the same loop on one preference triple first. Add 108 direct welfare-ledger
  evaluations and one common-cap finite-difference check. Expected under one
  minute; stop at ten minutes or any failed assertion. No equilibrium solve.
- Source proposal created. It gives an exact seller compensation function,
  an alternative small uniform premium, explicit debt repayment, the primitive
  stationary fertility restriction, and exact population accounting. Claims
  remain subject to checks and the independent review below.

### Verified checkpoint — 2026-09-04 22:13 EDT

- **Phase 1:** explicit exact compensating transfer, bond reserve for property
  tax, seller's estate reserve, buyer repayment out of old consumption, and
  unchanged government receipts and next-date housing supply. A small uniform
  seller premium is a simpler local alternative to exact household-specific
  compensation. These are checked conditional propositions, not final W3/W4
  instrument choices. Fixed moving costs require the separate finite-trade test.
- **Phase 2:** a sufficient threshold with no multiplier, plus a stationary
  exogenous-resource specialization with no price either. The primitive
  inequality is `(1-q)b/[(1-phi)(y+b)] >= alpha/(rho+alpha)` when property tax
  and its rebate are zero and old choices are interior. It is conditional on
  strict down-payment binding and maintained branches; it is not a GE sign
  theorem. The report includes a complete conditional-owner example satisfying
  saving, retention, and estate inequalities.
- **Verification:** 4,280 binding-cap configurations;
  796 satisfy the primitive sufficient condition
  and all have positive local responses. 548 other
  configurations have negative responses. Maximum derivative/finite-difference
  discrepancy is 6.69e-08. The 108 constructed
  welfare cases satisfy the budget/repayment checks; exact compensation keeps
  the seller indifferent and a uniform premium makes both parties better off.
  The common-cap derivative correction also agrees with finite differences.
  These are constructed household checks, not calibrated or GE results.
- **Reviewer:** session 39561 completed successfully, within its 20-minute
  limit. Both arguments passed the local, fixed-price/branch scope. Its two
  accounting clarity corrections (explicit tax-reserve asset and explicit
  seller title receipt) were adopted and checked by the lead. Review output:
  `output/model/simplified_olg_amendments/independent_review.md`.
- **PDF:** seven pages, compiled twice; no undefined references or overfull
  boxes. All pages visually inspected, with the transaction table inspected
  at full size. The analytical diagnostic figure was also inspected. An early
  summary-output attempt encountered a NumPy integer JSON serialization error;
  it was corrected, and the complete verification driver now passes in under
  two seconds. No model result was accepted from the failed output attempt.
- **Phase 4:** exact transition-product and terminal-scale identities, complete
  child-unit rescaling, and stationary supply accounting drafted. A resident-
  person counting convention remains an author choice. No quantitative units
  or targets were changed.
- **Next wakeup:** phase 3. Read the existing transition residual and local
  descriptor code; assess a broad continuation result, then run only the
  smallest justified tax-free diagnostic with explicit overrides and outputs
  under this workflow. Existing driver defaults use positive gains tax and
  overwrite older figures: do not invoke its main unchanged. Preserve those
  historical outputs. No worker is still active for this theory task.
- **Still unfinished:** broad transition existence/convergence/uniqueness;
  aggregate policy fertility and terminal-population signs; reconciliation of
  all 14 original repairs against the final source proposal; final morning
  summary. Do not pause the automation merely because the first PDF exists.

### Continuous work and writing revision — September 4 evening

The author rejected idle gaps between scheduled runs and asked for continuous
work until completion, with three checks of all substantive results. The
schedule is paused. Recent conversations reviewed: `Review transition closure
deck` and `Reconstruct theory and draft history`. The latter contains the
explicit instruction “make simple things simple” and says that simple theory
proofs should not announce an elaborate proof strategy. Retain the economics
and mathematical qualifications, but replace shorthand such as “maintained
branch” in explanatory prose with what it means: the same constraints bind.
Keep the proof details in the appendix and the work record out of paper prose.

The next numerical check keeps property tax at 0.05 and sets only the gains
tax to zero. Existing historical figures/results will be preserved. First
smoke the exact policy loop at financed shares 0.80 and 0.82 with horizon 8;
then, if correct and fast, examine predeclared shares 0.80, 0.81, 0.85, and 0.90
with horizon 28 and a longer-horizon check of one nonzero reform. These are
constructed theory examples, not calibration or a search for a positive effect.
Budget: at most eight initial transition solves and under ten minutes locally;
stop on an unexplained failed condition. Independently check the original
budgets, demographic laws, and housing/fiscal clearing, and use a second
numerical method for one path. Record a new plan before any additional run.

### Transition and accounting checks — September 4, 23:20 EDT

The exact two-case smoke loop passed before expansion. The four predeclared
credit reforms and the longer-horizon/second-method checks completed in 164
seconds. The baseline has gains tax zero and property tax 0.05; all cases use
the same initial old cohort, initial young/old masses, and fertility preference
0.35. Numerical results are illustrative theory examples, not calibration.

| Financed share | Initial fertility (model units) | Final price | Final household change |
|---|---:|---:|---:|
| 0.80 | 0.500000 | 0.342684 | 0% |
| 0.81 | 0.501891 | 0.355879 | +1.4297% |
| 0.85 | 0.503893 | 0.415235 | +5.9656% |
| 0.90 | 0.499167 | 0.445716 | +6.4316% |

All endpoints have replacement fertility 0.5. The 90% reform initially lowers
fertility: conditional fertility changes at initial tenure shares contribute
+0.019564, while the shift toward renters contributes -0.020398. The owner
borrowing limit is slack in that case. This prevents presenting the strictly
binding-cap theorem, or a smooth aggregate-map assumption, as valid throughout
that reform without further work. The decomposition is accounting, not a
causal experiment separating price and credit changes. Smaller policies also
have later dates below replacement. None proves a uniform positive policy sign.

The compensated allocation proof and the credit-policy example answer different
questions. The former specifies payments and holds fertility fixed; the latter
lets prices, tenure, and fertility adjust and does not establish Pareto gains.
The fertility theorem lets the young household reoptimize after a housing change;
it does not hold its future consumption allocation fixed as in the explicit
feasible welfare construction.

A conditional continuation proposition now covers finite reforms, with a true
infinite converging tail, if feasible solutions stay in a bounded region,
remain locally unique, and enter a feasible local stable graph. The proof uses
the implicit-function theorem and compactness. Those assumptions have not been
established over a policy interval from primitive conditions. The example's
endpoint eigenvalues and finite-horizon checks cannot supply that missing proof.

The old cohort also has a finite description for model-generated initial
conditions: with fixed entrant wealth and two-period lives, its distribution
is determined by the prices and policies facing it when young. An unexpected
reform must preserve pre-reform expectations in the inherited cohort. This
qualifies the older arbitrary-distribution warning; it does not invalidate it
for arbitrary initial distributions or models with inheritance.

The bounded transition reviewer (session 24755) completed. Its substantive
corrections were incorporated: aggregate differentiability is explicit;
model-generated initial distributions and pre-reform expectations are required;
the cohort price list includes two future prices; and the continuation function
is defined on an open neighborhood. The reviewer also suggested requiring all
off-solution terminal points to lie in the local stable neighborhood. The lead
instead requires the stable graph itself to be feasible and forward invariant:
at a zero of the matching equation the terminal point is exactly on that graph,
which suffices for the infinite tail. This choice is explicit for final review.

### Three-check record

| Claim | 1. Derivation | 2. Original equations or independent accounting | 3. Separate check |
|---|---|---|---|
| Fixed-fertility allocation improvement, exact compensation, small uniform premium | Finite utility difference and positive derivative; fixed-cost inequality derived separately | Buyer title/tax reserve/loan/resale ledger; seller estate bond; government receipt equality | First independent proof review and 108 constructed finite transactions |
| Fertility sign and primitive restriction | Implicit derivative; unconstrained expenditure-share inequality; full positive-saving owner example | Original two dated household budgets and Euler conditions; separate old-renter reduction | Independent proof review; 4,280 scalar household cases, 796 satisfying the restriction |
| Common rental-limit derivative | Chain rule includes next-period rent deduction | Original old-renter budget gives the $q^2 r_{t+1}$ deduction | Centered difference: -0.085122524329 vs analytical -0.085122524341 |
| Finite-reform continuation | Compactness and local unique continuation; stable graph supplies infinite tail | Finite-date matching equation and inherited-state timing checked against equilibrium definition | Independent transition proof review; local stable theorem checked against McGehee–Sander's original result |
| Heterogeneous cohort history | Distribution is the pushforward of fixed entrant types through two-period choices | Price list includes $P_t,P_{t+1},P_{t+2}$; original old states and pre-reform expectations retained | Independent transition review plus final synthesis review |
| Child units, population, supply | Exact change of variables, cohort-product identity, housing clearing | Original demographic and fiscal equations; household vs resident-person counts | Four direct unit rescalings, independent supply/person arithmetic, and product checks on saved transitions |
| Illustrative policy paths | Original market/household system, endpoint accounting | Every-date original budgets, FOCs, housing/fiscal clearing; 24 direct constrained household optimizations | Independent residual with least-squares solution, longer horizon, and endpoint derivative step check |

Numerical maxima: housing error 3.45e-10; government error 3.85e-12;
terminal-state gap 2.93e-10; 28-vs-40 path difference 2.13e-13;
second-solver log-path difference 1.10e-12; original household choice gap
1.07e-6. Unit-rescaling error is 1.78e-14; cohort-product error is 1.12e-15.
These tests establish the reported checks, not missing interval-wide hypotheses.

### Reconciliation of the original 14 repairs

The new consolidated proposal is the proposed replacement. The earlier theory
section, appendix, numerical example, and frozen independent review remain
historical evidence. No provisional choice has been promoted into the protected
manuscript, and the old builder's defaults have not been changed.

| Review row | Resolution and source location | Remaining author or proof obligation |
|---|---|---|
| 1 | Proposal allocation proposition and Appendix A separate market transactions from the added planner loan. Real resources and repayment are explicit. | Same-private-credit benchmark remains parked under W1. |
| 2 | Core proof charges compensation and real costs to the buyer's funded transaction. It does not reuse the disputed gains-tax/common-rebate sufficiency claim. | Positive-cost gains-tax extension remains parked under W2; no silent certification of the historical claim. |
| 3 | Main transition section and Appendix C report finite numerical paths, initial cohorts, terminal errors, and distinct household/population outcomes. | No quantitative endpoint or uniform aggregate fertility sign is certified. |
| 4 | Capital-gains tax is zero in proposed core and new illustrative calculations; property tax remains positive except the labeled primitive specialization. | Quantitative gains-tax inclusion is still Q0, parked. |
| 5 | Main fertility result is conditional; Appendix A includes the common young/old rental-limit resource effect and its verified derivative. | Other binding constraints require their own reductions. |
| 6 | Population section gives the complete child-unit transformation. Mapping 0.5 to 2.1 requires a factor 4.2, not 2.1. | Author chooses literal-child and resident-person interpretation; arithmetic itself is settled. |
| 7 | New environment gives all dated budgets, external finance, closing sequence, title/occupancy limits, warm-glow estates, government budget, and one equilibrium definition. | These clarify the existing institutions; broader institutional changes require I1. |
| 8 | Zero-reform path and initial steady-state linearization checked separately from finite-reform roots. Appendix B states the actual theorem hypotheses. | Endpoint numerical derivatives do not verify interval-wide smoothness or uniqueness. |
| 9 | No taxable gains or relevant acquisition basis in the new zero-gains-tax examples; the old evidence and its thresholded gains dates remain intact. | Restore and verify acquisition basis and meaningful gains dates if Q0 is reopened. |
| 10 | Main text places the fertility inequality beside its economics; Appendix A contains the proof, nonvacuity example, and adverse case. | No multiplier is needed; the fully primitive form uses the stated stationary specialization. |
| 11 | Main note includes a compact fertility/population figure; the full six-panel diagnostic is retained in supporting output. | The old wedge figure is not selected for the proposed core; final manuscript figure choice remains low priority. |
| 12 | The owner cap becomes slack in the largest example; direct original optimization checks both tenures at three dates in all four cases. Original equations checked at every date. | No claim of exhaustive corner coverage or arbitrary-type solver coverage. |
| 13 | Exact finite utility gains are available for the compensated transaction. No policy consumption-equivalent welfare number is inferred from fertility or population. | A broader policy welfare exercise depends on the W3/W4/I1 instrument and fiscal rules. |
| 14 | Appendix B supplies a conditional exact infinite-horizon continuation result and identifies its missing model-wide assumptions. | Primitive conditions covering all desired reforms, especially constraint crossings, remain unproved; no unconditional global theorem is claimed. |

### Final independent review — September 4, 23:24 EDT

The final `reviewer_strong` pass (session 46179) completed within its default
20-minute limit. It found no substantive mathematical, timing, or numerical
error. It checked the welfare transaction, fertility reoptimization, primitive
condition, original constraints, units and supply, finite-reform proof, cohort
history, and the saved numerical findings. It also confirmed that feasibility
of the local stable graph suffices at a zero of the matching equation.

The lead adopted its one wording correction: equation residuals and constraints
are checked at every date, whereas direct original-household optimization is
checked at dates 0, 2, and 28. No new mathematical claim was added after review.
The report is `output/model/simplified_olg_amendments/final_independent_review.md`.

### Read this before the next discussion

The completed proposal establishes a local allocation gain at fixed fertility
under the stated financing and compensation, and an explicit conditional
fertility prediction, including a primitive stationary specialization. The
illustrative credit reforms show larger terminal household populations even
though fertility can fall at some transition dates. The finite-reform theorem
is conditional: its model-wide assumptions, especially at constraint crossings,
have not been proved. Neither population growth nor a household partial effect
is presented as a general policy welfare or fertility theorem.

Start the next discussion with **W4**, now a concrete choice: exact compensation
or a small per-unit seller premium. Then settle **W3**'s real cost and **I1**'s
institutional interpretation if needed. **U0**'s arithmetic correction is done;
the literal-child and resident-person conventions still need interpretation.
No quantitative calibration, target, model unit, or protected manuscript was
changed. All 18 decision IDs and all 14 repair rows are retained.

The hourly automation remains paused. This development pass ends with the
consolidated proposal and explicit unresolved choices; it does not silently
adopt provisional assumptions or claim the missing global results.

### Delivery verification — September 4, 23:30 EDT

The final 12-page PDF compiled twice with no LaTeX warnings, undefined
references, or overfull boxes. All pages were visually inspected; after the
last formatting changes the unchanged page images were byte-compared and
every changed page was inspected again. The transaction table and the enlarged
figure legend are readable. The completed final independent review found no
substantive error. Its one wording correction was adopted.

The analytical driver and saved-path accounting passed after the final plotting
change. The only subsequent source change clarified the dates covered by direct
optimization. The original transition solve, second method, and longer horizon
were not rerun merely to refresh prose hashes; the delivery manifest records
the numerical function hashes separately. No model equations changed after the
completed checks. All 18 map IDs and the ordered 14 repair rows were validated.
The hourly automation remains paused and no reviewer remains active.

### Additional prose pass — September 4, 23:47 EDT

Tommaso authorized any safe remaining work before going to sleep. This pass
simplified the opening, the explanation of payments in the compensated sale,
and several appendix passages. It also kept a subsection heading with its
explanation across a page break. The specific Menzio and Fernández–Rogerson
passages read for exposition are now recorded in Section 9 of the writing guide.

Three checks covered these edits: a reading of the prose diff for changes in
meaning; a comparison against commit `509e773` confirming unchanged inline
mathematics, displayed equations, propositions, definitions, tables, and labels;
and two clean LaTeX compilations followed by visual inspection of all 12 final
pages. All cross-references resolve. The existing numerical files and code
match their delivery hashes. The inspected PDF replaces the previous delivery,
and the manifest records this prose pass separately from the numerical checks.

No economic choice or formal result changed. There was no new model run,
theory review, protected-manuscript edit, or scheduled work. This completes the
bounded prose pass; the next discussion still starts with W4's compensation
choice, and the automation remains paused.

### Author conventions restored — September 5

Tommaso objected that the proposal had changed presentation and notation he
had chosen himself. The earlier section and appendix,
`latex/simplified_olg_paper_theory_section.tex` and
`latex/simplified_olg_paper_theory_appendix.tex`, are the comparison sources.
The protected manuscript contains a placeholder and was not edited.

| Object | Unnecessary change in the proposal | Restored convention |
|---|---|---|
| Preferences | One expanded lifetime-utility expression | Separate flow utilities $u_t^Y(c,h,n)$ and $u^2(c^2,h^2,e)$ |
| Household problems | Budgets followed by unnamed optimization and renamed young values | Explicit $W_t^R(i),W_t^O(i),V_t^R(a),V_t^O(a,H)$ problems |
| Old-age choices | $c_2,h_2$ and $c_{2i}$ | $c^2,h^2$ and $c_i^2$ |
| Tenure | $\pi_t$, $\sigma$, inverse-logit expression | $\pi_t^O$, $\sigma_\xi$, and the original exponential-share expression |
| Entrant distribution | $F_I$ | $F$; the new appendix matching function uses $\mathcal F$ to avoid taking this name |
| Housing values and costs | $m_Y,m_O,k_c$ | $MV^Y,MV^O,K_{ijt}$; local household indices are suppressed |
| Lifetime reduction | $\rho$, $h$ as the housing choice | $1+k$, with $k=\beta A$ or $\beta(1+\omega_B)$, and the original $x,s,n$ choice presentation |
| Markets and averages | New user-cost shorthand throughout, changed housing-average subscripts | $R_f=q^{-1}$, the original rental-intermediary equation, $\bar h_t^Y,\bar h_t^O$ |

The environment, household problems, and equilibrium are separate again.
The earlier housing-menu inequality and rule against renting retained owner
housing are stated explicitly, as in the pre-amendment source. The reduced
household comparison suppresses indices locally, with the correspondence stated
before its equations. The agreed removal of gains tax from the core, planner
credit, conditional fertility result, and qualified transition result remain.

The restored household objectives and constraints were read against both the
earlier source and the previous proposal. Substituting $R_f=1/q$ gives the
same dated budgets; the two tenure-share formulas are algebraically identical;
and $s=h-\kappa n$ gives the same reduced budget and feasible housing limit.
An automated comparison against `20d9c0b` verifies all four propositions,
the equilibrium definition, all 16 result displays, and all tables after the
documented notation substitution. All old labels remain and references resolve.
The final PDF compiled successfully, with no warnings or layout errors in the
last pass, and all 13 pages were visually inspected. Supporting numerical code
and results retain their earlier hashes; no numerical model was rerun.

The correction is recorded in the live discussion map, writing guide, and
cross-session memory. Future economic amendments must preserve the author's
existing conventions unless a specific necessary change is proposed separately.
