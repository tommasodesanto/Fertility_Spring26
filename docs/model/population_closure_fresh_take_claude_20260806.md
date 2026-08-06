# Population closure: independent assessment (Claude fresh take, 2026-08-06)

Response to
`docs/prompts/PROMPT_population_closure_fresh_take_20260806.md`, written under
that prompt's protocol: the assessment below was frozen before reading the
companion survey memo (`docs/model/population_closure_literature_memo_20260806.md`);
a reconciliation coda follows at the end. Disclosure: while locating files I saw
the memo's first ~60 lines (its Section 1 headline that no published paper makes
this population claim). Nothing below defers to it.

Numbers used throughout: \(E_0=0.062\) (required entrant flow per unit
population per four-year period), \(B_0=0.053\) (matured city-born flow),
\(s^{out}=0.169\) (outside-origin share of entrants), \(q^*\approx0.969\)
(baseline entry probability implied by the identity), and the production policy
readouts recorded in `memory/daily/2026-08-06.md`.

## Verdict

The accounting core of the current construction is sound and, after the
August 6 normalization repair, empirically disciplined: two observable numbers
— the outside-origin entrant share 0.169 and the \(S=1\) scale normalization —
pin everything the baseline actually uses. The choice-theoretic wrapper around
that core (the representative outside value \(\bar V^{out}\), the taste scale
\(\lambda=2\), the logit) is decoration at the baseline and actively harmful in
counterfactuals: \(\lambda\) is unidentified, it generates 80–90 percent of the
current production population responses, and the near-saturated margin
(\(q^*\approx0.97\)) makes those responses mechanically asymmetric. The author's
long-standing suspicion of the disembodied outside value is correct, and the
right response is deletion, not repair.

Recommendation in one line: adopt menu item (iv) — stationary level equilibrium
closed by a policy-invariant outside inflow — as the primary closure; report
population effects as renewal arithmetic times a general-equilibrium housing
damping term the model itself supplies; demote the logit protocol to a labeled
sensitivity row. This changes no baseline moment, deletes code rather than
adding it, and moves the honest headline population effects from roughly +10
percent to roughly +3 percent for the current experiments.

## 1. The menu conflates three questions

The six protocols answer three logically separate questions as if they were
one:

1. **Existence**: what device gives the model a stationary equilibrium at all,
   given endogenous below-replacement fertility? (Some inflow must top up the
   missing 14.5 percent of each entering cohort.)
2. **Invariance**: which closure object is held fixed in counterfactuals — the
   outside value, the entry probability, the inflow level, or the population
   scale?
3. **Statement**: what population-level sentence does the paper write?

Once separated, most of the menu collapses. Items (ii), (iii), and (iv) are
**one model reported in different units**, not rivals. Under any closure that
holds the outside inflow's level fixed, the stationarity identity
\(S E_0 = \bar M + \bar R\, S B_0\) (with \(\bar R\) the native entry/retention
rate and \(\bar M\) the outside inflow) gives

\[
\frac{d\ln S}{d\ln B_0} \;=\; \frac{\bar R B_0}{E_0-\bar R B_0}
\;=\; \frac{1-s^{out}}{s^{out}} \;\approx\; 4.92 .
\]

The same derivative can be quoted three ways, and all three are the same
number: births +0.67 percent (the grant experiment) implies (a) stationary
population +3.3 percent at fixed immigration, (b) required immigration −3.3
percent at fixed population, and (c) 3.3 percent of the fertility gap closed.
Protocol (ii) reports (a); protocol (iii) reports (b); protocol (iv) is the
implementation both run on. The only genuinely different model on the menu is
(i), which adds a second term
\(\frac{1}{s^{out}}d\ln q\) driven by the unidentified \(\lambda\).

Three framing errors follow:

- **"Population must respond" does not mean "entry must be behavioral."**
  Desideratum 2 is fully satisfied by (iv): population responds through the
  margin the model actually disciplines — births — amplified by transparent
  renewal arithmetic. The author rejected (iii) as "population must respond,"
  which is right for the headline; but the fixed-\(S\) unit (b) should survive
  as a reporting row, because "this pronatal policy substitutes native births
  for immigrant headcount roughly one-for-one" is a bulletproof,
  policy-relevant sentence that no referee can attack, and it is closure-
  invariant.
- **The local-versus-national ambiguity is prior to the closure choice**, and
  the menu smuggles it inside the protocols: (i) is an open-city device, (ii)
  is a national device, and choosing between them by taste hides that they
  answer different questions. Settle the interpretation first. Everything
  below assumes the national reading, which the prompt says is intended.
- **The (ii) ≡ (iv) equivalence claim is wrong as a claim about the world**,
  even though the numbers coincide. See next.

### Protocol (i) is unsalvageable as a headline

Three independent reasons, any one sufficient:

1. \(\lambda\) is a free parameter that generates most of the result. From the
   production record: floor-arm tax-only population +9.98 percent, of which the
   arithmetic term \(4.92\times d\ln B_0\) explains about +0.9; the remaining
   ~+9 points are \(\frac{1}{s^{out}}d\ln q\) with \(d\ln q\) set by
   \(\lambda=2\). Roughly 90 percent of the headline is an unidentified
   parameter. The project's own SMM identification rule — never let an
   unidentified object carry a result — applies to counterfactual objects with
   the same force as to estimated ones.
2. Near-saturation makes the response set asymmetric by construction. With
   \(q^*=0.969\), the maximum possible upside is
   \(\frac{1}{s^{out}}\ln(1/q^*)\approx+18\) percent of population, while the
   downside is unbounded. Policy conclusions should not inherit a one-sided
   response window from a normalization.
3. Under the national reading the experiment concept is incoherent: holding
   \(\bar V^{out}\) fixed while a *national* policy moves all inside values is
   the single-treated-city thought experiment. The paper is one market.

Keep (i) only as a labeled upper-bound sensitivity row ("if entry responded to
values with \(\lambda=2\), an assumption we cannot identify, population effects
would triple"). The existing runs are then not wasted.

### (ii) versus (iv): same numbers, different — and differently defensible — stories

They coincide for aggregate national experiments: both end up holding entry
rates and the inflow fixed. They are not equivalent as economics:

- (ii)'s justification ("the outside is more of the same economy; a nationwide
  policy moves \(\bar V^{out}\) in parallel, so the gap and hence \(q\) are
  unchanged by optimality") requires the outside to be *replica US cities*.
  That contradicts the national one-market reading, where the outside is
  abroad: a US property-tax reform does not move foreign values in parallel.
  The author's "moving the goalposts" discomfort is detecting a real
  inconsistency.
- (iv)'s justification is institutional and survives scrutiny: US lawful
  immigration is quantity-rationed and oversubscribed. The marginal entrant is
  determined by the quota, not by indifference; a policy-induced value change
  of well under one percent consumption-equivalent does not move the binding
  constraint, so \(\partial(\text{inflow})/\partial V\approx0\) locally is a
  statement about institutions, not a behavioral assumption. (Caveat for
  honesty: the rationing argument is cleanest for lawful permanent admissions;
  the empirical drivers of other components are wages and enforcement, both
  exogenous and unchanged here.)

So the right primary closure is (iv), *justified by rationing*, with (ii)'s
parallel-shift story retired. Removing the choice-theoretic layer is a feature,
and not merely "honesty about what is identified": the layer was never
identified even at the baseline — \(\bar V^{out}\) and \(\lambda\) were backed
out, never tested, and the data cannot distinguish the logit model from the
quota model at all. They differ only in counterfactuals, which is exactly where
unidentified structure is most dangerous.

### (v) fails by arithmetic, not by taste

The prompt instructs me to challenge the claim that a balanced path sacrifices
the scarcity margin. On examination the claim is right, and stronger than
stated: **a stationary-in-per-capita-terms path with population growth
\(n\neq0\) is inconsistent with any aggregate housing supply curve in levels.**
With \(H^S=H_0(r/\bar r)^\xi\) and constant \(r\), the stock \(H^S\) is
constant while \(S\) declines geometrically, so per-capita housing grows
without bound; rents must then trend, contradicting balancedness. A flow
(construction) version fails the same way. The only escapes are (a) supply
proportional to population (kills desideratum 3), (b) \(\xi\to\infty\) (kills
it too), or (c) a land-microfounded per-capita margin, which reintroduces the
trend under \(n\neq0\). Conclusion: given the level supply curve the author
wants to keep, a stationary population level (\(n=0\)) closed by an inflow is
*forced*, not chosen. That resolves the menu's biggest item on logical grounds
and makes the closure much easier to defend in a seminar: aggregate long-run
scarcity and permanently nonzero population growth are incompatible
desiderata; the paper picks scarcity.

Two further points on stance:

- The immigration-closed stationary state is not a modeling trick; it is the
  current US demographic regime. Natural increase is near zero and recent
  Census vintages attribute the large majority of US population growth to net
  international migration. The model reproduces precisely this: natives below
  replacement, an outside inflow making up the difference.
- Household-level fertility results (prompt question 1a) are best supported by
  the existing fixed-population decomposition rows — household responses at
  unchanged aggregates — with the (iv) general-equilibrium rows as the
  headline. Current practice already does this; keep it.

## 2. Off-menu: the statement the model can own

The genuinely new content is not a new closure; it is that **the model's own
scarcity margin disciplines the population statement**. Under (iv), the
stationary system is

\[
S=\frac{\bar M}{E_0-\bar R\,B_0(r)},\qquad
H_0(r/\bar r)^\xi = S\cdot \bar h(r;\text{policy}),
\]

and a policy that moves births directly by \(d\ln B_0^{pol}\) moves population
by

\[
d\ln S=\frac{m\, d\ln B_0^{pol}}{1+m\,|\varepsilon_{B,r}|\,\kappa},
\qquad m=\frac{1-s^{out}}{s^{out}}\approx4.92,
\]

where \(\kappa=d\ln r/d\ln S=1/(\xi+|\varepsilon_{h,r}|)>0\) comes from the
level supply curve and \(\varepsilon_{B,r}<0\) is the fertility–rent
elasticity — the paper's core object. The chain is: more births → more
households in the stationary state → permanently higher rents (supply is in
levels) → slightly fewer births. The knife-edge renewal amplification
(\(m\approx5\), a hyperbola in the native replacement ratio) is tamed by
endogenous housing scarcity, and the amount of taming is the model's
contribution. This serves desiderata 2 and 3 simultaneously: population
responds, and the scarcity margin is not just preserved but load-bearing. At
the current calibration's small fertility–price responses the damping is
modest (order 10 percent; one extra stationary solve computes it exactly), so
the arithmetic dominates — say so plainly rather than overselling the GE
content. The macro sentence the paper can own: *a reform that raises births
0.7 percent raises long-run population about 3 percent and long-run rents
about \(3.3\kappa\) percent, so housing scarcity makes pronatal policy partly
self-limiting.* The rent leg has an external validation point: the identified
population→rent elasticity in the immigration literature (Saiz 2007: an
immigration inflow of 1 percent of city population raises rents about 1
percent) can be compared to the model's \(\kappa\).

Other off-menu items, in decreasing order of usefulness:

- **Report the trio, always.** (a) \(\Delta S\) at fixed inflow, (b) required
  inflow at fixed \(S\), (c) percent of the fertility gap closed. One number,
  three units; makes clear it is renewal accounting; unit (b) speaks directly
  to the births-versus-immigration policy debate.
- **Report the convergence half-life.** Stationary comparisons hide that the
  renewal gap closes at rate \(1-\bar R B_0/E_0\approx0.17\) per generation:
  half-life on the order of a century. Computing it is pure demography at
  frozen prices (two price bounds — old and new stationary prices — bracket
  the path with no forward-looking solve). Reporting it is what makes
  excluding full transitions (desideratum 1) intellectually respectable rather
  than convenient, and it protects against the referee who multiplies the
  stationary claim by a misread time horizon.
- **Closure-device doctrine (the capital-market analogy the prompt asks
  about).** Schmitt-Grohé and Uribe (JIE 2003) showed that small-open-economy
  models need a stationarity-inducing device (debt-elastic premium, etc.), that
  several devices deliver nearly identical local behavior, and that the honest
  practice is to pick the simplest and show device-robustness. The outside
  value was this model's device pretending to be economics. Frame the closure
  section exactly this way and show the trio is invariant across
  (ii)/(iii)/(iv).
- **If a price-elastic entry margin is ever wanted, the identified one is
  household formation, not immigration.** The model already carries children
  at home (\(m\)) with a maturation hazard (2/9 per period). Making the
  maturation-to-household transition respond to rents is the *headship*
  margin, with a real empirical literature (young adults' living arrangements
  respond to housing costs and labor conditions; e.g. Ermisch, JUE 1999;
  Kaplan, JPE 2012 on moving back home; the ACS doubled-up literature). It is
  the economically right margin for a housing paper — expensive housing delays
  household entry and, through the fecundity clock, births — and it changes
  the stationary count of *households* even with the count of *bodies* closed
  by the quota. This is real surgery (entry-age distribution feeds
  everything), so it is a future-work or robustness item on a two-month
  runway, not the closure. But it converts the author's "entry should respond
  to prices" instinct into an identified object instead of \(\lambda\).
- **Land microfoundation: one paragraph, no code.** \(H_0(r/\bar r)^\xi\) in
  levels *is* the reduced form of fixed land plus elastic structures; say so
  (Saiz-type geography), so that "population permanently moves rents" reads as
  fixed-factor economics rather than an arbitrary functional form.
- Considered and rejected: perpetual-youth/stochastic-aging recasts (destroy
  the lifecycle timing and fecundity architecture the model is about);
  Barro–Becker dynastic closure (requires dynastic altruism; the warm-glow
  architecture and the entire calibration are built against it); endogenous
  mortality (wrong margin); marriage/partnership states (previously rejected
  by the author for this model).

## 3. What the literatures actually do

Branch by branch, with confidence labels. Verified here = I am confident of
venue and content from my own knowledge; survey = taken from the project's
two-scout survey memo, which reports reading the papers directly; WP = working
paper.

- **Quantitative endogenous-fertility policy models are partial equilibrium.**
  Doepke and Kindermann (AER 2019) and Sommer (JME 2016): exogenous wage
  processes, no market clearing; policy holds prices fixed (verified). The
  Handbook of the Economics of the Family chapter (Doepke, Hannusch,
  Kindermann, Tertilt 2023) keeps household fertility and macro-demography in
  separate sections (survey).
- **Demographic GE / pension-aging OLG models take fertility as exogenous
  projections** (Auerbach–Kotlikoff lineage; for below-replacement Japan:
  Kitao JEDC 2015, Braun and Joines JEDC 2015 — verified as to exogeneity;
  the detail that Kitao manufactures a terminal steady state by assuming
  fertility recovers is survey). Immigration, where present, is an exogenous
  policy-set flow: Storesletten (JPE 2000) has the government choose age- and
  skill-specific inflows (verified). Nobody models national entry as a
  value-based choice against an outside option.
- **Urban and spatial models fix the national total.** Open-city versus
  closed-city (Rosen 1979; Roback JPE 1982): the open city takes reservation
  utility as given precisely because the closed national system pins it
  (verified). Quantitative spatial models (Redding and Rossi-Hansberg, ARE
  2017; Monte, Redding, Rossi-Hansberg AER 2018; Hsieh and Moretti AEJ:Macro
  2019) normalize national population and reallocate it (verified).
- **Structural housing–fertility papers fix the count.** Greaney,
  Parkhomenko, Van Nieuwerburgh (NBER WP): fixed unit mass, exogenous
  turnover, entrants with zero wealth, bequests leave the economy (survey;
  WP). Moreno-Maldonado and Santamaría (CEPR DP): demographic composition is
  an exogenous input; timing moves, the count does not (survey; WP).
- **The exception that proves the rule**: Barro and Becker (Econometrica
  1989) close population through fertility itself, but only via dynastic
  altruism (verified). Golosov, Jones, Tertilt (Econometrica 2007) study
  efficiency with endogenous population, not quantitative closure (verified).
- **Reduced-form ancestors of the two legs the model combines**: demographics
  → housing demand → prices (Mankiw and Weil, RSUE 1989, verified);
  immigration → rents (Saiz, JUE 2007, verified); housing costs → fertility
  (the project's own reduced-form anchor literature).

Implication: there is no convention to inherit — every adjacent branch either
fixes the count, fixes fertility, or drops prices. Interpreted through the QSM
lens, closure (iv) is the standard national boundary condition (exogenous net
international migration) with one new endogenous component, natural increase.
That is exactly where the novelty should sit, and it is why the closure should
be maximally boring: the paper's contribution is the births margin and its
housing feedback, and every unit of structure spent on the entry margin is
structure a referee will (correctly) attack.

## 4. Identifying an entry elasticity: mostly, don't

Split by interpretation:

- **National (the intended reading).** The object would be
  \(\partial(\text{net inflow})/\partial(\text{US value change from housing
  policy})\). No credible moment exists: lawful inflows are quota-rationed
  (elasticity locally zero); gravity-style estimates of international
  migration flows (e.g. Ortega and Peri) are identified off cross-country
  income differences orders of magnitude larger than a sub-percent
  consumption-equivalent policy change, and mapping them through gives an
  effect indistinguishable from zero at this scale. The honest answer is the
  one the prompt anticipates: a preregistered sensitivity band, with
  elasticity zero (quota) as the point specification, a gravity-implied small
  positive value as the middle, and the legacy \(\lambda=2\) logit as the
  labeled extreme upper bound. Three rows, existing machinery, no new claims.
- **City (not the current paper).** Here elasticities are identified —
  cross-city migration and preference estimates (Diamond AER 2016; QSM
  gravity parameters) — but importing them into a one-market national model
  is exactly the treated-city error of protocol (i). If the paper ever goes
  spatial, the closure question reopens with real data behind it.
- **The identified housing-side moment is the reverse one.** Saiz (2007)
  identifies population → rents, which validates the damping leg \(\kappa\),
  not the entry leg. Use it that way.
- **Re-anchor \(s^{out}\) under the national reading.** The 0.169 is
  across-CBSA arrivals at ages 18–22 — a *city* object. The national analog
  is the foreign-origin share of new young households (ACS nativity, same age
  window; plausibly in the same 0.14–0.20 range, but it must be measured, not
  assumed). This is a half-day empirical task and it matters doubly, because
  feasibility of the closure is the knife-edge condition
  \(s^{out}\ge 1-B_0/E_0\) (equivalently \(q^*\le1\)), which the production
  driver now correctly enforces per candidate.
- **Do not over-interpret \(q^*=0.97\).** Since
  \(q^*=(1-s^{out})E_0/B_0\), its distance from 1 moves one-for-one with the
  fertility block's fit: a 3 percent shortfall in \(B_0\) saturates the
  margin; hitting the completed-fertility target instead would put \(q^*\)
  near 0.90. Near-saturation is a property of the current calibration's
  fertility fit as much as of the economy. (Another reason the logit, whose
  local slope lives entirely off \(1-q^*\), should carry no weight.)

## 5. Minimal implementation of the recommended closure

**Objects.** Replace \(\{\bar V^{out},\lambda,\text{logit}\}\) with two
calibrated constants: \(\bar R\) (native entry/retention rate) and \(\bar M\)
(outside inflow level). Baseline: \(\bar R=(1-s^{out})E_0/B_0\) — numerically
the current \(q^*\) (floor arm 0.969225, tilt 0.969844) — and
\(\bar M=E_0-\bar R B_0\) so \(S=1\). The baseline equilibrium is bit-identical
to the current one; zero recalibration.

**Counterfactual.** Hold \((\bar R,\bar M)\) fixed. Solve jointly for the
stationary scale and prices: \(S\,E_0=\bar R\,S\,B_0(\theta';r)+\bar M\) and
housing-market clearing with supply in levels (plus the balanced-budget
transfer, as now). This is the existing funded driver minus the logit block —
a code deletion, roughly an afternoon including reruns of the current policy
table.

**New assumptions, stated in the paper.**
(A1) The outside inflow is a policy-invariant level — quota rationing;
defended institutionally and by the QSM national-closure convention.
(A2) The native entry rate is policy-invariant at 0.969 (its complement, 3
percent early out-migration, is small; sensitivity trivial; note it cannot
exceed 1, so no saturation asymmetry arises).
(A3) Comparisons are across stationary states; convergence is generational —
report the renewal half-life.
(A4) The supply curve is invariant in levels (fixed land), so population
changes move rents permanently. This is desideratum 3 doing the damping work.
Also note the maintained simplifications that all entrants share the
native entrant wealth/age profile and that warm-glow estates transfer to no
one, so population changes carry no mechanical wealth flows.

**Reporting contract.** For every policy: the trio (population at fixed
inflow / required inflow at fixed population / share of the fertility gap
closed), the arithmetic-versus-GE-damping decomposition, the permanent rent
effect of the population change, the half-life, and the entry-elasticity
sensitivity band with the \(\lambda=2\) row explicitly labeled unidentified.
Retire the +9 to +10 percent figures as headlines; the honest headline for
the current experiments is roughly +1 percent (tax-only, births +0.19
percent) to +3.3 percent (grant arms, births +0.67 percent), damped slightly
by the housing feedback.

**What this buys.** Desideratum 1: no transitions, and the half-life makes
that respectable. Desideratum 2: population genuinely responds, through the
disciplined margin. Desideratum 3: scarcity is load-bearing, not just
preserved. Desideratum 4: the suspect object is gone; the closure rests on
two measured numbers. Desideratum 5: weeks, mostly deletion plus one ACS
re-measurement.

**If the current construction were kept instead** (the prompt asks me to say
so if I think it is fine): the baseline is fine — it always was, being
observationally equivalent to the quota model. The counterfactual protocol is
not fine, and no amount of sensitivity analysis rehabilitates a headline that
an unidentified parameter multiplies by ten.

---

## Coda: reconciliation against the project survey memo

*Written after the above was frozen, upon reading
`docs/model/population_closure_literature_memo_20260806.md` in full.*

**Convergence.** The survey memo and this take reach the same primary
recommendation independently: closure (iv), exogenous inflow held fixed,
population responding through births with multiplier
\((1-s^{out})/s^{out}\approx4.9\), logit/\(\bar V^{out}\)/\(\lambda\) retired,
fixed-population rows kept as decompositions, the same Storesletten and
US-demographic-regime defenses, the same verdict that the open-city device is
a category error for a one-region national model, and the same conclusion
that no credible identifying moment exists for a national entry elasticity.
Two independent assessments landing on the same closure is itself evidence
the author can use.

**What this take adds on the memo's two open questions.** (1) On
equivalence: (ii), (iii), and (iv) are one model in different reporting
units — the trio — so they are equivalent in every reportable object for
aggregate national experiments; but they are not equivalent as claims about
the world, because (ii)'s parallel-shift justification is internally
inconsistent under the national reading while (iv)'s rationing justification
survives. Only (i) is a different model. (2) Off-menu additions: the GE
damping fixed point (the scarcity margin tames the renewal multiplier; the
model's own contribution to the population number, with Saiz JUE 2007 as
external validation of the rent leg — a different paper from the Saiz QJE
2010 supply elasticities the memo cites for the supply curve); the
convergence half-life (~a century) as a required caveat statistic; the
quantified indictment of protocol (i) (roughly 90 percent of the current +10
percent headlines is the unidentified \(\lambda\) term, with a one-sided
response window); the observation that \(q^*=0.97\) saturation moves
one-for-one with the fertility block's fit and would be ~0.90 at the
completed-fertility target; the recommendation to re-anchor \(s^{out}\) to a
national object (foreign-origin share of new young households, ACS nativity)
rather than the across-CBSA 0.169; and the headship/maturation margin as the
identified price-elastic entry margin for future work.

**One genuine tension.** The memo's Section 3 keeps \(M\) calibrated to the
across-CBSA 0.169 "already adopted in May"; I argue that under the national
reading this is a city-level object and the anchor should be re-measured as
the foreign-origin share of new young households (with the feasibility
condition \(s^{out}\ge1-B_0/E_0\) re-checked, since the closure sits within
a few percent of that knife edge). Half a day of ACS work either confirms
0.169 or corrects it; either way the paper's anchor becomes the right object.

**What the memo adds that this take missed.** The welfare trap: with
different population sizes across counterfactuals, welfare statements must
stay per-household for fixed cohorts — never total-utilitarian across
allocations with different sets of people. This belongs in the reporting
contract of Section 5 above. The memo also surfaces two unread field-journal
leads (Okamoto 2021; Hagiwara 2025) that should be skimmed before the paper
claims novelty for endogenous-fertility-plus-immigration closure, and the
"two children form one household" mapping, which confirms the \(B_0\)
accounting used here.

---

## Coda 2: comparison with the second independent fresh take (2026-08-06)

A second independent take (GPT, supplied by the author) reached the same
destination: retire the logit/\(\bar V^{out}\)/\(\lambda\) layer; close the
national model with a fixed absolute flow of net international immigrants;
keep fixed-population rows as the household-level decomposition; treat the
national entry elasticity as unidentified, with zero (a fixed immigration
regime) as the baseline and a sensitivity band otherwise; read (ii)'s
"parallel movement" as an assumption a one-market model cannot deliver as an
equilibrium result; reject perpetual youth; and distinguish the
household-formation margin from international migration. Three independent
assessments now agree on the closure.

**Adopted from the second take** (it is better on empirical anchoring
discipline):

1. **Sequencing.** Its sharpest point: the 0.169 outside-origin share is a
   cross-CBSA *domestic* migration statistic, and domestic moves net to zero
   nationally; until \(\bar M\) is re-anchored to net international migration
   in entrant-household units, the 4.92 multiplier and the +3.3 percent
   headline are transparent algebraic sensitivities, not nationally
   calibrated estimates. Section 4 above flagged the re-anchoring as a task;
   the second take correctly makes it a *precondition* for quoting the
   numbers.
2. **An explicit \(\kappa\).** Make the child-to-entrant-household mapping an
   explicit external object (survival \(\times\) two-children-per-household
   \(\times\) headship \(\times\) retention/emigration), externally
   restricted and reported — never silently one. Do not back out both
   \(\kappa\) and \(\bar M\) from \(S=1\): keep \(S=1\) as the units
   normalization for \(H_0\), pin \(\kappa\) externally, and use Census net
   international migration as an *overidentification check* on the implied
   \(\bar M\). This converts the closure from calibrated-to-fit into
   data-disciplined-and-testable.
3. **Cleaner identification algebra.** With a common \(q\) across pools,
   \(s^{out}=qM/[q(M+B)]=M/(M+B)\): the observed share identifies pool
   composition only; \(q\) itself is an accounting residual of stationarity,
   and no baseline moment touches its elasticity. Sharper than the version
   in Section 1 and worth quoting in the paper.
4. **Consumption-equivalent units for the band.** Write the entry
   sensitivity as \(I(\Delta c)=I_0\exp(\varepsilon_I\Delta c)\) rather than
   \(\lambda\) in model-value units (its point that \(\lambda(1-q^*)/s^{out}
   \approx0.355\) per "value unit" is empirically meaningless is correct);
   convert the legacy \(\lambda=2\) runs into their implied \(\varepsilon_I\)
   so they remain usable as the labeled extreme row.
5. **Naming.** Call the object stationary *households*, not population,
   unless a persons mapping is added.
6. **A canonical cite.** Espenshade–Bouvier–Arthur (1982): constant
   immigration plus below-replacement fertility converges to a stationary
   population — the demographic theorem under this closure. Verify against
   the paper before citing (project rule); the second take's other
   migration cites (Kennan–Walker 2011; Notowidigdo 2020; Zabek 2024;
   Fehr–Jokisch–Kotlikoff) are standard.

**Where this assessment stands its ground:**

1. **The balanced path.** The second take calls (v) "feasible in principle,"
   claiming fixed land can preserve scarcity on a BGP. Maintained
   disagreement: with any level (non-proportional) supply object — fixed
   land included — \(n\neq0\) forces trending per-capita land and hence
   trending rents; balancedness fails generically (the knife edge
   \(n=-\delta\) with zero construction aside). Scarcity on a BGP requires
   supply proportional to population or stock, which is exactly what
   desideratum 3 refuses. The impossibility argument in Section 1 stands.
2. **The GE damping loop is absent there.** Its closure preserves scarcity
   in levels but never notes that scarcity feeds back to births and damps
   the renewal multiplier — the model-owned part of the population statement
   (Section 2) and the main reason the model, not demography alone, is
   needed for the number.
3. **Half-life, saturation asymmetry, the \(q^*\)-as-fertility-fit
   diagnosis, the welfare-across-populations rule, and the trio units** are
   all absent there and all kept. In particular its claim that the
   fixed-population protocol "makes no defensible statement about
   population" misses that the required-immigration-offset row is the same
   number as the population row in different units.
4. **The local open-city appendix** it proposes (with a freshly estimated
   origin-destination migration elasticity) exceeds the runway, and a
   one-market model calibrated to national moments cannot be reinterpreted
   as one metro without re-anchoring supply and the entrant pool. Default
   remains the sensitivity-row demotion; the appendix is an option only if
   an elasticity is credibly imported and the author accepts the metro
   reinterpretation.
5. **The age-structured fixed point** \(\mu=[I-A]^{-1}\iota\) it recommends
   is, with all entrants at age 18, exactly what the stationary solve
   already computes — the scalar closure is its reduced form. It becomes
   substantive only if immigrants enter at multiple ages with their own
   composition; worth a stated caveat, not implementation.

**Net.** No change of destination across three takes. The union is the
stronger spec: quota closure with external \(\kappa\) and national
re-anchoring as a precondition (second take), plus the damping
decomposition, half-life, trio units, and CE-unit sensitivity band (this
take). The decision remains the author's.
