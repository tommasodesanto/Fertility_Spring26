# Simplified OLG: development plan and overnight record

## September 6 author response — full note still to be updated

Tommaso assigns the new mixed-tenure transition and stationary welfare work
to the appendix as secondary material for now, and returns the priority to
housing misallocation. The four-page `simplified_olg_paper_core.pdf` is a
results extract plus an assessment. It is not the full rewritten note.
The earlier `simplified_olg_amendment_proposal.tex` contains the full setup,
household problems and equilibrium, but has not been updated to incorporate
all subsequent discussion. That integration remains unfinished; this turn
records the distinction and placement decision, without producing another
partial note or changing the PDFs. Preserve the author's conventions, all
other decision branches, and the earlier two-figure preference.

## September 6 ten-hour proof window — proof and verification completed

Tommaso authorized a more ambitious proof while sleeping. Start 03:28 UTC
September 6; hard stop 13:28 UTC (09:28 New York), or earlier when a meaningful
result and its checks are complete. The target is an analytical local
transition at an interior mixed-tenure economy with material renter mass and
positive child goods costs, rather than only perturbations of the all-owner
zero-cost limit. Retain original household equations, endogenous tenure,
physical restrictions, actual initial old claims and the accepted author
notation. No public-finance permission or preference change is inferred.

Stages and stop conditions: (1) at most 90 minutes to derive and verify a
finite-dimensional exact equilibrium map and its initial conditions;
(2) at most three hours for a primitive or exact-certificate stability and
boundary argument at a nondegenerate mixed economy; (3) at most three hours
for fertility/population signs and a carefully delimited extension, only if
the preceding proof is valid; (4) reserve at least 90 minutes for independent
review, original-equation checks and a compact amendment to the existing
reading material. At a failed stage, record the precise obstruction before
changing method; do not repeat the same broad investigation. Long numerical
searches would require a separate smoke-tested Torch design, but none is
planned. This is primarily analytical work with short arithmetic checks.

The lead owns the model reduction, construction, proof and synthesis. One
bounded read-only reviewer at a time will check genuinely new steps, each
with its default profile limit; no broad parallel duplicate proof search.
Keep progress and all proved/conjectured distinctions in this record and the
existing evidence folder. Preserve the short allocation result and all 18
decision entries while this optional extension is explored.

Completed the scientific work and PDF inspection by 04:44 UTC, about 76 minutes
after starting. Stop early because the selected proof and checks are complete;
the remaining paper-placement and institutional choices belong to the author.
No continuing theory job or scheduled run is created. Final backup details
are recorded in the daily memory.

The result is materially broader than the earlier all-owner limiting proof.
An exact six-variable equilibrium representation retains original positive
child costs and rebated property tax, endogenous logistic tenure and the
actual initial-old distribution. The initial old-owner housing coefficient
must revalue inherited title and the surprise rebate; it cannot be fixed as
if it were predetermined financial wealth. A positive-tax local inversion
lemma applies when old owners occupy at least as much as capped old renters.

The exact constructed economy has `11/21` owners, `chi=3/20`, positive tax
`467/9250`, and a finite taste scale. It admits a locally unique converging
transition with higher initial fertility and a larger final population after
a small permanent credit relaxation. Exact rational roots and the initial
boundary, followed by a one-sided sequence theorem, prove the nonlinear
infinite-path claim. The noninvertible stable zero root causes no obstruction.
An explicit uniform tail lemma also proves eventual fertility and population
crossings for every sufficiently small finite reform near this economy.

The family extension holds the same baseline bundles at different taste
scales by specifying the corresponding taste location in primitives; both
taste parameters are held fixed within each credit reform. A rank-one matrix
identity yields four stable/two unstable roots at every positive scale.
Rational intervals certify the full scale interval `[1,4]`, with 913 adjacent
cells, for boundary transversality and positive initial fertility. Stationary
population rises and both conditional lifetime values fall at every positive
scale in this family. Thus the same stationary entrant loses welfare at any
fixed taste draw, despite population growth. Initial-old and full-transition
welfare remain separate. The reform radius and a general multidimensional
parameter box are unquantified; branch changes and global transitions are open.

Three sequential read-only reviews, each using the default 20-minute strong
profile limit, checked distinct work: exact aggregate reduction and initial
claims; rational one-point certificate and finite-reform argument; newly
broadened family and stationary welfare signs. All completed within their
limits. The second requested explicit uniform parameter/remainder bounds;
the added lemma closes that issue, as confirmed by the third. The third
reviewer audited source and saved receipts but could not replay them in its
model environment because `sympy` was absent. The lead ran the entire checker
in the working environment and independently checked the algebra and original
equations.

The initial full-family interval attempt and a compiled-arithmetic version
each stopped at their four-minute bounds: direct repeated interval
dependencies made the all-family oscillation test too conservative. A changed
method used a scale-independent rank-one eigenvector polynomial and retained
the oscillation claim near scale one, concentrating the family certificate
on initial fertility and final population. That completed in 84 seconds in
the prototype and about 91 seconds in the final driver, within its explicit
600-second/12,000-subdivision fail limits. No calibration search was launched.

One new driver, `code/model/tools/verify_simplified_olg_mixed_transition.py`,
records all exact bounds and reproduction metadata. The four final receipts
match source SHA-256
`f22bdc662ca470b763805abe0adec28cb1548c0b921260e57636bd84be08daa4`.
The old limiting driver is unchanged. Four short original-equation paths and
24 independent household optimizations pass; maximum equilibrium residual
is below `2.9e-15`, budget error below `1.6e-15`, and initial/final derivative
discrepancies below `3e-9`. The 24/40-date price paths agree within `3.2e-15`.
Four separate original stationary comparisons confirm the all-scale rational
derivatives to `6.9e-9`. These numerical checks support the algebra; they do
not supply the plotted curves or the infinite-path proof.

The final backup check removed trailing whitespace from the new driver.
Its complete parsed Python syntax tree is unchanged. All four receipt modes
were replayed after this formatting-only change; numerical fields and the
plotted PNG remain identical to the reviewed version. The final source hash
above includes that cleanup.

The revised `output/pdf/simplified_olg_paper_core.pdf` remains four pages;
its first two source pages are unchanged. Figure 2 now uses the analytical
mixed response, retaining the allocation/transition pair of figure roles.
All earlier figures and limiting results are preserved. The new supporting
appendix, `output/pdf/simplified_olg_mixed_transition_proof.pdf`, is eight pages.
Both documents compile without warnings or overflowing boxes and every page
was visually inspected. PDF hashes respectively:
`e6298f04410475bb23eda9d1f34a239e77ef6059a0306697b0bb5d5d295feb62`
and `b8b950597dc32bb0d9c1093e1f1d80741b53fc68cdbc9b04ec46427d812e7c7c`.
All 18 decision rows remain unchanged. No author convention, preference,
planner permission, quantitative-model source, or protected manuscript file
is changed. Proofs, reviews and commands are indexed in the existing evidence
folder. Recommend keeping the short allocation result central and deciding
the transition's paper placement together.

## September 5 late-evening scope and prose pass — completed

Authorized when Tommaso returned tired and asked for safe autonomous work.
Start: approximately 02:23 UTC September 6 (September 5 in New York).
Bounded objective: a three-page proposed results section plus one page of
economic assessment, preserving the existing setup, notation, preferences and
both figures. Assess the substantive transition assumptions and independently
check whether the plotted stationary welfare loss extends across its limiting
family. One read-only `reviewer_strong` pass has the profile's 20-minute limit;
the lead owns the derivation, text, original-equation checks and final judgment.
Stop when the compact text, checks and PDF inspection are complete, with a
90-minute ceiling for this pass. No new calibration, global transition search,
public-finance design, planner permission, or author-manuscript edit is planned.

Completed at 02:49 UTC, about 26 minutes after starting. The deliverable is
`output/pdf/simplified_olg_paper_core.pdf`, with source under
`latex/JMP_DS_suggestions/`: three pages of proposed results after the
unchanged household setup, followed by one page of assessment. Both existing
figures are preserved. The longer ten-page note remains available for its
proofs, counterexamples and institutional discussion.

The new result explains the limiting steady-state welfare loss generally,
including heterogeneous entrant income and wealth. At the all-owner
zero-child-cost/property-tax limit, increasing the common financed share
raises price proportionally, leaves every young type's housing and fertility
unchanged, and lowers every old type's housing. Stationary population rises
while each same-type entrant's owner lifetime utility falls. The full short
derivation is section 11 of `local_transition_proof.md`; no claim is made
about transitional-cohort welfare or a social ranking of different populations.

The economic assessment separates the mild-looking stability inequality from
near-all ownership, small child goods costs, interior old choices, fixed
entrant endowments and a world bond price. The explicit estate restriction
is stated in primitives. The independent review's phrase about estates
requiring retention is qualified: the original estate inequality is not a
minimum-retention rule. No change to it or to ownership tastes is proposed.
Recommend the allocation and conditional fertility arguments as the core,
with the transition presented as a separate illustrative extension. This is
a recommendation, not acceptance of local scope as the author's final goal.

Three checks cover the new welfare claim: direct derivation from the original
dated budgets and separate flow utilities; one completed read-only independent
review; and symbolic identities plus original-budget/inequality checks for
five fixed types at three financed shares, including ten independent household
optimizations. The utility derivative discrepancy is below `8.3e-9` and budget
errors below `4.5e-16`. Every prior economic/numerical/plotting function is
AST-identical. All five earlier numerical receipts reproduce exactly except
for source hashes and elapsed time; all six current receipts match the final
driver. Figures were not regenerated or redesigned.

Revisited Menzio's JPE 2007 model exposition (pp. 751--753) and the actual
Fernández--Rogerson QJE 1996 prose (pp. 137--138). No model assumptions were
borrowed. Final TeX has no warnings, undefined references or overflowing
boxes. All four pages were visually inspected; the final change is confined
to page 4, with pages 1--3 byte-identical to the inspected renders. Final PDF
SHA-256: `882ced5e6bd997a14d6a8cdba4e60bc8f192f17c66840ccc1eec8bd2088e8f5a`.
All 18 decision rows are unchanged. Source and evidence backup is recorded
in the daily memory. No worker or numerical job remains active from this pass,
and no scheduled run was created.

## September 5 continuation before the author's 9pm return — completed

Start: 21:33 UTC, authorized when Tommaso said he would not be back before
9pm America/New_York (01:00 UTC September 6). The analytical work and document verification were completed at 23:12 UTC;
the source/output backup outcome is recorded in the daily memory note. The earlier
eight-page simplicity pass remains a historical checkpoint. Its reading note
is now ten pages, with a derived equilibrium Figure 2 and primitive transition
conditions. All work remains outside the author-controlled manuscript.

### Results and scope

The new proof in `output/model/simplified_olg_amendments/local_transition_proof.md`
starts from the exact original dated budgets in an all-owner, zero-child-cost,
zero-property-tax demand limit. It retains initial old assets and titles rather
than reoptimizing their past choices. Its example admits a local nonlinear
converging transition, initial and all-date above-replacement fertility, and a
larger final population. A uniform tail bound proves the finite-reform all-date
sign in this example; pointwise positive derivatives alone would not suffice.

The general result is stronger than the initial numerical parameter example.
Writing all coefficients from primitives as in proof section 10, local
convergence has the exact root condition
\[
(1-q)+(3+q)C>4D.
\]
With the maintained uniformly strict household branches, the initial cohort
and actual pre-reform old demand select a unique local converging path. A
separate exact initial-fertility test evaluates an explicit cubic at a point
computed from primitives. Under \(q\ge1/2,r>1\), the simpler sufficient
condition is
\[
(1-q)C[(1-q)+(1+q)C]>L[(1-q)-D+C].
\]
These are primitive conditions, not conditions on an unknown price or borrowing
multiplier. The local theorem permits heterogeneous entry: owner quantities
aggregate exactly by income and wealth means at the limit, and integration is
smooth under uniform original private conditions. An explicit nondegenerate
support is verified in section 9. The theorem admits small positive renter
mass, child goods costs and property tax by a weighted infinite-sequence
implicit-function argument that retains endogenous tenure and its added price
leads. Neither the allowed neighborhood nor the reform size is numerically
certified. Arbitrary heterogeneity crossing constraints and global continuation
remain outside the theorem.

Initial fertility and terminal population are separate signs. A second
admissible economy has lower initial fertility and a larger final population.
The all-date monotonic population result is therefore limited to the plotted
example, not the whole family or the positive-renter extension. Separately,
section 8 gives an explicit stationary population derivative with positive
child costs and examples with both signs. The uncompensated credit reform
lowers stationary household welfare in the plotted limiting economy:
\(\partial_d W^O=-289/475\). This is not the direct compensated Pareto
comparison and does not implement it automatically.

### Verification and review record

- Lead derivation against the original dated utility, budgets, demographic law,
  normalized entrant distribution, tenure rule, inherited old state, and rebate.
- Symbolic recurrence and boundary differentiation, rational root/sign bounds,
  two algebraic stability arguments (direct roots and a unit-circle transform),
  exact primitive stationary derivative, and uniform infinite-tail proof.
- Five sequential read-only `reviewer_strong` passes, each with the profile's
  20-minute limit and a distinct new deliverable: initial example; newly closed
  operator/tail gaps; heterogeneity and positive-cost stationary formula;
  general convergence under \(w>d\) and exact/sufficient fertility tests;
  the final sharper stability condition, including zero/repeated roots. Each
  completed within its 20-minute limit. The last report independently confirms
  the section 10.4 algebra, alongside the lead’s second proof and symbolic/root
  receipts. All requested scope and operator qualifications are incorporated.
- One deterministic driver, no long search: eight finite original-equation
  cases, three stationary sign cases, four-corner conditional/integration
  checks, 270 declared coefficient cases including unstable ones, and a
  negative-initial-fertility original-model counterexample. Finite terminal
  closures only check arithmetic; they never substitute for the infinite-path
  proof. The largest equilibrium residual is below \(2.3\times10^{-14}\),
  central derivative error is \(1.48\times10^{-9}\), and 24/40-date paths
  agree within \(4.45\times10^{-16}\).
- The deliberately tiny rental cap caused one unscaled SLSQP check to fail
  its success criterion. Scaling numerical choice units fixed the conditioning
  without changing the original objective or constraints. Every final mode
  passes and its JSON records the final driver hash.
- Figure 2 now plots the analytical first-order formula, including the
  second-cohort fertility hump. Earlier figures are preserved. The core
  allocation figure and original flow/value notation remain unchanged.

The final ten-page PDF was compiled without warnings or overflowing boxes,
rendered, and inspected page by page. Its final SHA-256 is
`d4f9316baead87946dc6af5f069493f710291f6ede40e362fc057028ee327c86`.
Evidence and regeneration commands are indexed in
`output/model/simplified_olg_amendments/README.md`. No planner permission,
ownership preference, population-welfare objective, or quantitative-model
change is adopted. The new theorem is a proposed theory extension for discussion.

## September 5 simplicity pass — completed

Authorized while Tommaso is on a three-hour run. Start: 20:38 UTC; finish by
23:38 UTC or earlier when the assessment, independent checks and compact PDF
are complete. The aim is a problem-set-sized allocation argument, followed by
a ranked discussion of extensions. No further author choice is inferred.

- Core: verify the short compensated reallocation, the roles of both access
  constraints, and a nonempty mixed-tenure equilibrium case.
- Extensions: reconcile the meaning of constrained efficiency, the existing
  public-finance construction, the fertility threshold, and demographic
  identities. Do not turn the later finance branch into the core result.
- Figures: an exact compensated-allocation diagram and a clearly conditional
  demographic-transition diagram; retain all earlier figures unchanged.
- Deliverable: approximately six to eight pages, source
  `latex/JMP_DS_suggestions/simplified_olg_simple_assessment.tex`, PDF
  `output/pdf/simplified_olg_simple_assessment.pdf`. Supporting evidence stays
  in the existing `output/model/simplified_olg_amendments/` folder.
- Route: one read-only `reviewer_strong` worker, profile limit 20 minutes,
  reviewing the direct proof and roles of the constraints. Lead owns the
  economics, independent household/equilibrium arithmetic, extension synthesis
  and final prose. No model search or policy simulation is planned.
- Three checks: analytical derivation, original-equation/feasibility checks,
  and the bounded independent review; then compile and visually inspect the PDF.

### Assessment

The simple allocation result works. Keep the existing household problems and
ownership taste. A positive group of young owners can value extra housing more
than old owners value keeping it; current-goods compensation then produces a
local Pareto improvement with fertility and future allocations fixed. This is
the direct benchmark already allowing relaxed private finance, not a full
first-best characterization or a claim about an unspecified constrained
planner. A binding physical rental cap alone supplies no feasible improving
direction for a renter who must retain that tenure.

New analytical constructions close the nonemptiness question. They include
mixed-tenure steady states with both rental caps binding, and another economy
with both rental caps slack but a positive owner value gap. Ownership tastes
make the latter possible. Thus rental segmentation can limit an alternative
route to space without being necessary for this particular inefficiency.
No preference change is adopted merely to make both restrictions essential.

The fertility derivative and its price-free sufficient specialization remain
separate conditional predictions. A common relaxation of young and old rental
caps includes both dated rent costs. Population accounting and stationary
stock scaling are established; a broad equilibrium transition following a
credit reform is not. The new second figure prescribes fertility and applies
the demographic identity. It is not a solved equilibrium path.

The existing committed-transfer proof and one-time-gift obstruction retain
their distinct, narrow scopes. Public collection of future income supplies
financing across ages. The obstruction does not prove that loans are necessary;
an anonymous loan programme is not proved by renaming identity-based transfers.
These remain later branches, with no new author acceptance of planner powers.

### New derivations retained for reproducibility

**Owner-to-renter replication.** With zero gains tax and the same prices and
rebates, set
\[
a_R'=a_O'+qP_{t+1}h-\phi P_t h.
\]
Then \(a_R'+u_th=a_O'+[(1-\phi)+q\tau^p]P_th\), and
\(a_R'/q=a_O'/q-\phi P_th/q+P_{t+1}h\). These reproduce the young
budget and the owner's old resources, including title value. The rental menu
must contain both housing choices and replicated saving must be nonnegative.
Without an ownership benefit, values agree; if current rental space can rise,
a small change \((\mathrm dx,\mathrm dh)=(-u_t\epsilon,\epsilon)\)
improves utility when \(MV^Y>u_t\). With an ownership taste, the material
replication need not match utility. This is a conditional deviation, not a
general-equilibrium cap-removal experiment.

**An explicit family with binding rental caps.** Use \(P=1,q=.5,\phi=.8\),
\(b=.2,\bar H=2,\nu=2\), \(\alpha=\beta=\omega_B=.4\),
\(\gamma=.3,\chi=.15,\kappa=.5\), rental cap .25, owner cap 2.
Select young-owner \((x,h,n)=(1,1,.75)\); its fertility condition sets
\(\vartheta=.3525\). For each \(0\leq\tau^p<.28\), let
\(u=(1+\tau^p)/2\), \(w=1.7925+u\), and
\[
(c^2,h^2,e)=(.8,.24/u,.64),\qquad MV^Y=.64>MV^O=u.
\]
The old owner has strict retention and estate slack. Both renter caps bind
at the unique root \(n^R\in(0,.5)\) of
\[
\frac{.3525}{n^R}=\frac{.15}{x^R}
 +\frac{.2}{.25-.5n^R},\qquad
x^R=\frac{1.7925+.625u-.15n^R}{1.56}.
\]
The left side strictly decreases from infinity, the right side strictly
increases to infinity, and \(x^R>2.03/1.56>1.30\). Consequently young
renter housing value exceeds 2.08 and old renter housing value exceeds 1.248,
both strictly above \(u<.64\). These choices have the right cap multipliers.

Let \(\pi^O=(.5-n^R)/(.75-n^R)\),
\(d_h=\pi^O(1+.24/u)+(1-\pi^O).5\), and choose
\[
Y=O=2/d_h,\quad T=q\tau^p d_h/2,\quad
y=w-b-(1+q)T,\quad
\bar\xi=\sigma_\xi\log[\pi^O/(1-\pi^O)]-W^O+W^R.
\]
Any \(\sigma_\xi>0\) supports the required owner share. Average fertility
is .5, housing and fiscal budgets clear, and the inherited old distribution
is generated by these same young choices. The entrant distribution is a point
mass in income and wealth, with the original continuous ownership taste.
Since \(d_h<1.48\), \(T<.1036\), \(y>1.9371\), and owner saving
\(a'_O=.98-.5T>.9282\). Renter saving is
\(a'_R=.56x^R+.125u-.5T>.7387\). The original objectives are concave
on convex feasible sets. The fertility, saving, estate, housing and cap KKT
conditions therefore establish global household optimality, not feasibility
alone. Income and taste location vary across this family: it is an existence
construction, not a tax comparative static holding primitives fixed.

**A complete equilibrium with slack rental caps.** Keep the other displayed
primitives, but set \(\tau^p=0\), rental cap 1.5, and young-owner
\((x,h,n)=(1,1,.49)\). Set
\(\vartheta=61397/302000\), \(y=2.0535\), and \(w=2.2535\).
Writing \(\rho=1+\beta(1+\gamma+\omega_B)=1.68\) only in this
supporting derivation, the unconstrained renter has
\[
x^R=\frac{w}{\rho+\alpha+\vartheta},\quad
n^R=\frac{\vartheta x^R}{\chi+\kappa u}
=\frac{276716279}{551645600}>.5>.49,\quad
h^R=\frac{\alpha x^R}{u}+\kappa n^R.
\]
The same share and stock construction gives \(\pi^O=.139389713\),
\(Y=O=1.325073899\), \(h^R=1.040368344<1.5\), and
\(h^{2R}=.473735108<1.5\). Owners still satisfy
\(MV^Y=.529801325>.5=MV^O\). At taste scale .12 the needed location
is \(\bar\xi=-.217848432\). Removing the redundant cap leaves this
second equilibrium unchanged. Comparing it with the first construction is not
a causal rental-cap reform: preferences and income also differ.

**Fertility condition.** The short derivative agrees with the earlier checked
proposal. A direct sufficient-condition proof is also available without an
unconstrained-demand detour. Put \(\rho=1+k\), \(w=\rho x+\chi n+uh\).
Strict purchase rationing implies \(\alpha x>us\). If
\(uh/w\geq\alpha/(\rho+\alpha)\), then
\[
0<\alpha\rho x-\rho us
 =\alpha w-(\rho+\alpha)uh+n(\rho\kappa u-\alpha\chi)
\]
forces \(\rho\kappa u>\alpha\chi\). Hence
\(\rho\kappa x>\chi s\); for any \(p<MV^Y=\alpha x/s\),
\(\alpha\kappa/s^2-\chi p/(\rho x^2)>0\). Stationary zero-tax
purchase spending is \(uh=(1-q)b/(1-\phi)\), giving the displayed
primitive restriction. No aggregate policy sign follows. The original-owner
negative-response example added to the new verifier uses
\(q=.5,\phi=.8,\beta=.4,\gamma=.5,\omega_B=1,\alpha=\chi=1\),
\(\kappa=.1,\vartheta=1.1,b=.22,y=3.33,P=1,\tau^p=0\), and
\((x,h,n,a',c^2,h^2,e)=(1,1.1,1,1.33,.8,.8,1.6)\). It is an
original household optimum with strict credit rationing and a negative
fertility response when \(p=u=.5\), not a policy equilibrium.

**Population.** The cohort ratio is the product of fertility ratios from a
common initial state; total household population includes the lagged old
cohort. At a positive steady state, \(\bar n=1/\nu\) and
\(N_{hh}=2\bar H/(\bar h^Y+\bar h^O)\). Under the actual fixed-entrant,
world-bond, ordinary-rebate closure, scaling \(\bar H\) and all cohort
masses by \(a>0\) leaves prices, rebates and individual problems unchanged.
This is a statement about corresponding stationary equilibria. An additional
supplier-income/finance channel, an endogenous entrant distribution, or
domestic bond clearing could invalidate it. For the new schematic,
\(\log(\bar n_t^1/\bar n_t^0)=(1-\delta)\log(a)\delta^t\), so
\(Y_t^1/Y_t^0=a^{1-\delta^t}\). The choice \(a=1.2,\delta=.6\)
illustrates accounting only. It proves no equilibrium transition or speed.
Stationary person counts require common adults-per-household and fixed rules
for the timing and residence of children; transitional children must be counted.

### Checks and review scope

Reproduce the new receipts and figures with:

```sh
python3 code/model/tools/verify_simplified_olg_simple_assessment.py --plot
```

The driver checks four constructed equilibria (three points of the proved
tax-indexed family and the slack-cap example), nine independent optimizations
using all seven original dated choices, 36 finite reallocations including
real-cost cases, replication, fertility finite differences, and demographic
identities. Maximum original equilibrium residual: \(1.30\times10^{-15}\);
maximum optimizer choice discrepancy: \(1.19\times10^{-6}\). The JSON pins
the driver and original specification hashes. There is no calibration, policy
simulation, numerical equilibrium search, or claim that these examples fit data.

Two sequential read-only `reviewer_strong` passes completed within their
20-minute limits. The first verified the direct proof and role of each
constraint; the second independently derived replication, the equilibrium
constructions and stationary scaling. Their corrected/preserved reports are
`simple_assessment_direct_review.md` and
`simple_assessment_new_claims_review.md` in the evidence folder. The first
report's accidental inclusion of K in the old-value definition and one broken
path were explicitly corrected; no economic conclusion was changed. The
second review's tax-family and person-count qualifications were incorporated.
Its scaling display was corrected to retain the normalized old-state
distribution G; the original specification scales the cohort mass, not that
probability distribution. The draft and driver already use this convention.
Earlier fertility proof review and original receipts remain independently
available; the new negative-response original optimization is a lead check,
not attributed to the second reviewer.

For style, the lead read Menzio's actual *A Theory of Partially Directed
Search* model section and the Fernández--Rogerson NBER abstract. The latter's
full PDF was inaccessible, so no full-text read is claimed. These are prose
references, not authority for this model's welfare or transition claims.
The author-controlled manuscript, original notation, original figures, and
quantitative model remain unchanged.

Final PDF verification: eight pages, two clean compilation passes after the
last text edits, no undefined references or overfull boxes, and all eight
rendered pages visually inspected. Figures, equations, table, citations and
page breaks are readable. The final PDF is in `output/pdf/`; source and
verification artifacts are ready for the author’s discussion.

## Earlier September 4–5 amendment pass


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
