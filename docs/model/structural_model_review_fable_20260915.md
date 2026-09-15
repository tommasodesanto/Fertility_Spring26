# Structural review of the housing–fertility OLG model

Reviewer: Fable 5.1, repository-reading reviewer. Date: September 15, 2026.
Brief: `docs/prompts/MASTER_PROMPT_model_review_20260915.md`.

What was read: `CALIBRATION_STATUS.md` (September 15 top, the shock-timing
section at line 3416, the August 15 architecture decision at 7312),
`docs/model/POST_PRESENTATION_ISSUES.md` (M01–M39),
`docs/model/ACTIVE_DECISION_LEDGER.md`, `latex/september_14_presentation.tex`
with its algorithm appendix, the July design notes (equivalence-scale
specification, sequential fertility, fecundity, bequests, income risk), the
September notes on choice timing, nested choice, population closure and the
property-tax audit, and the live package
`code/model/intergen_eqscale_seq_optimized/` with the transition drivers under
`code/model/tools/`. No model code was run. No file other than this one was
created or changed. Every code claim carries a path and line; Appendix A lists
the prompt's twelve code claims with verdicts and the corrections they need.

Notation: \(u\) flow utility, \(\mathcal C\) the aggregator, \(m\) children at
home, \(n\) children ever born, \(\bar h(m)\) the parents' space floor,
\(e(m)\) the equivalence scale, \(\kappa_H\) the housing taste scale
(`tenure_choice_kappa`), \(\kappa_1,\kappa_C\) the fertility scales
(`kappa_fert`, `kappa_fert_continuation`), \(\xi\) the first-birth cost
(`first_birth_fixed_cost`), \(\mu\) the per-period child exit probability.

## Summary judgement

The household block is coherent and better grounded than the slides say. The
preference form is the Scholz–Seshadri–Khitatrakun scale plus a Stone–Geary
space floor for parents. The fertility block is a priced biological clock with
a fecundity schedule within three points of Léridon's four-year conception
probabilities. The mortgage block is the standard one-period-bond
approximation. Four things are structurally wrong or undefended, and all four
sit outside the household problem:

1. **Supply is a reversible static schedule** \(H_0P^\eta\). A price fall
   destroys stock within the period, and a shrinking population has no
   stationary limit because \(P\to0\). This is why the transition has no
   terminal steady state and why the tax experiment shows housing services
   falling at impact. Fix: stock-flow supply with depreciation and
   price-elastic construction (H4, D2).
2. **Child maturation is memoryless.** With \(\mu=4/18\), 22 percent of
   newborns leave within four years, 8 percent of children born at parental
   age 25 are still dependents at 65, and the dependents of dying households
   vanish. Fix: a parent-age-dependent exit hazard that reaches one by age 62
   plus a newborn flag. This closes the orphan accounting without joint death
   or a care pool (F3, F4).
3. **The earnings process is two half-processes.** A Floden–Lindé AR(1)
   scaled to after-tax units is added to a PSID gross-earnings fixed effect
   from a different decomposition. Not double counting; not one process
   either (E1).
4. **The housing taste shock is a numerical device at a bound whose value
   flips the sign of the policy's ownership response.** It must be zero or
   estimated (F2).

Everything else is keep-and-explain or an author decision with two defensible
positions. The largest omission not on the list is a time cost of children,
which is why the model cannot produce the negative income gradient of
fertility (X3).

## Preferences

### P1. CRRA over an aggregator plus a linear child term

**(a) At stake.** The aggregator decides how children change the demand for
space at a given budget; that is the policy's transmission. The child term
decides why anyone has children and how many; that is what the policy is
measured on.

**(b) Defensible?** The implemented form is
\[
u=\frac{1}{1-\sigma}\left[\frac{c^{\alpha}\,(s-\bar h(m))^{1-\alpha}}{e(m)}\right]^{1-\sigma}+\psi\,m,
\qquad e(m)=\Big(\tfrac{2+0.7m}{2}\Big)^{0.7},\qquad \bar h(m)=h_P\,\mathbf 1\{m>0\},
\]
`kernels.py:64` (renter), `kernels.py:127` (owner), `solver.py:2253-2278`.
The scale exponents 0.7 and 0.7 are hard-coded literals (`solver.py:2278`),
not parameters. The per-child rooms slope `hbar_child_rooms` is restricted to
zero and \(h_P=2.3\) sits at its upper bound 2.3 in the September 13 vintage
(`output/model/e5f_final_night_20260913/corrected_initial/parameters.csv`).
\(\psi\) multiplies children at home, not children ever born
(`solver.py:2254`), so parents get nothing from adult children, consistent
with the child-blind bequest.

Each piece has a precedent: the utility scale is Scholz, Seshadri and
Khitatrakun (2006); the parents' housing floor is the Stone–Geary device of
the housing literature; CRRA over a Cobb–Douglas composite is the workhorse of
every housing lifecycle model. The combination of a scale and a floor has no
precedent I know of, and it has a rationale that the July 20 specification
note already wrote down: a Barten-type scale cancels from the Cobb–Douglas
allocation, so \(e(m)\) alone gives no rooms response; the floor was re-added
(the "floor arm") for that reason. The paper should say exactly this: \(e(m)\)
prices children in welfare and saving, \(\bar h(m)\) prices them in space.

The linear \(\psi m\) is the weak part. Benefits are linear in \(m\). Costs
are concave in \(m\): the scale increments for the first three children are
0.234, 0.216, 0.203, and the floor is a one-time step. So the net flow benefit
of the \(k\)-th child rises in \(k\), and a household without taste shocks
wants zero children or the biological maximum. The observed distribution of
family sizes (one-child share 21 percent among mothers) is generated entirely
by \(\kappa_C\), the fecundity decline and the seven-period window. Smoothing
a discrete choice this way is legitimate, but the identification slide's
"\(\kappa_C\to\) one-child families" reads to a referee as "the intensive
margin is a noise parameter". Linear is not without precedent: Doepke and
Kindermann (2019) give a linear one-time utility from each birth, but pair it
with a fixed goods cost, a fixed utility cost and a time cost per child that
rises with the wage, so their marginal cost of children is not concave.
Sommer (2016) and Barro and Becker (1989) put curvature in the child term
(Appendix B).

**(c) Minimal change.** Replace \(\psi m\) by \(\psi\,v(m)\) with \(v\)
concave, \(v(m)=\log(1+m)\) or \(m^{1-\varepsilon}/(1-\varepsilon)\), with the
curvature fixed externally so that the free-parameter count is unchanged.

**(d) Discriminating observable.** The number of children among mothers by
household income, and the three-plus share. Linear benefits with concave costs
imply an intensive margin that rises with income; a concave benefit with an
absolute space floor implies a flatter one.

**(e) Recommendation.** Author decision. Position 1: keep linear and state that
\(\kappa_C\) and the clock carry the intensive margin. Position 2: concave
\(v(m)\), which makes the intensive margin an economic object. Cost low; blocks
recalibration.

### P2. Sign and role of \(\mathcal C_{sm}\)

**(a) At stake.** The mechanism needs children to raise the value of space at
the margin where a young household decides whether to buy an owner-only size.
That margin is a discrete comparison of housing products, not a derivative.

**(b) Derivation at the implemented form.** Write
\(X=c^{\alpha}(s-\bar h)^{1-\alpha}\) and \(\mathcal C=X/e(m)\). Then
\[
u_s=(1-\alpha)\,e(m)^{\sigma-1}\,\frac{X^{1-\sigma}}{s-\bar h(m)},
\qquad
\frac{u_s}{u_c}=\frac{1-\alpha}{\alpha}\,\frac{c}{s-\bar h(m)} .
\]
Three facts follow. First, on \(u\) both channels raise the marginal utility of
space when \(\sigma\ge1\): the scale through \(e^{\sigma-1}\), the floor
through \(s-\bar h\). For \(\sigma<1\) the scale channel reverses and at
\(\sigma=1\) it vanishes. The slide's bullet is true at \(\sigma=2\), on
\(u\), not on \(\mathcal C\); M24 is right that \(\mathcal C_{sm}\) is
ambiguous because \(e\) lowers \(\mathcal C_s\). Second, the marginal rate of
substitution between space and goods does not depend on \(e(m)\) at all. At
expenditure \(E\) and unit price \(q\) the interior demand is
\(s=\bar h(m)+(1-\alpha)(E-q\bar h(m))/q\), so a child moves rooms by
\(\alpha\,\Delta\bar h\) and nothing else: with \(\alpha=0.733\) and
\(h_P=2.3\), 1.69 rooms for an unconstrained renter. The entire space channel
is the floor. Third, \(e(m)\) works elsewhere: at \(\sigma=2\) it multiplies
the marginal utility of the composite by 1.23 for the first child, so a
household that expects a birth has higher future marginal utility and saves
more before it; this is where a down-payment accumulation motive comes from.
It also sets the welfare cost of a child that the fertility logit weighs
against \(\psi\).

**(c) Minimal change.** None to the model. Write the three facts. The realized
rooms response of 1.01 against the floor-implied 1.69 is the constrained
version: renters at the six-room cap and owners facing the 6 percent sale cost
cannot expand. That gap is the mechanism, and it should be shown.

**(d) Observable.** The rooms response to a first birth by prior tenure and
wealth: unconstrained renters near \(\alpha h_P\), capped renters and
low-wealth owners below it.

**(e) Recommendation.** Keep. Cost none.

### P3. The first-birth cost \(\xi\)

**(a) At stake.** \(\xi\) shifts the first-birth gain by \(\pi_a\xi\) at every
fertile age (`solver.py:2718-2729`). It is the only intercept of the
first-birth decision: \(\psi\) is fixed by the 2.1 normalization, \(h_P\) by
the rooms moments, \(\kappa_1\) by the age spread.

**(b) Defensible?** It identifies a level; \(\kappa_1\ne\kappa_C\) identifies
dispersions. In a binary logit with gain \(\Delta\) and scale \(\kappa\), a
larger \(\kappa\) pushes the attempt probability toward one half from either
side, so \(\kappa_1\) cannot lower attempts for households whose gain is
positive; only a cost can. So it is a childlessness parameter by construction
(M30 says why: \(h_P\) was taken by the rooms moments), and it is not
disguised, just unexplained. The closest precedent is Doepke and Kindermann
(2019), who carry a fixed utility cost \(\phi_u\) for every child alongside a
fixed goods cost and a time cost (Appendix B). Their cost is per child and sits
next to a wage-scaled time cost; ours is first-child only and stands alone, so
it cannot borrow their justification without the other two pieces.

**(c) Minimal change.** Replace it with something measured. The candidate is a
per-period earnings cost of children at home, \(y^d\,(1-\tau_c)\) while
\(m>0\): mothers' earnings fall persistently by about 20 percent after a first
birth in Denmark (Kleven, Landais and Søgaard 2019), more in the United States
(Kleven et al. 2019; verify the U.S. figure before use). This cost scales with
the wage, so it lowers first births more for high-\(z\) households and
produces the negative income gradient of fertility (Jones and Tertilt 2008
document the U.S. cross-section; the mechanism is Becker's time cost).
Dropping \(\xi\) with no replacement leaves childlessness to \(\kappa_1\) and
fails.

**(d) Observable.** Childlessness by permanent income or education (CPS
children ever born). \(\xi\) is income-neutral in utility units and predicts a
flat or positive gradient; an earnings penalty predicts a negative one.

**(e) Recommendation.** Author decision, lean change: add the measured child
penalty, then test whether \(\xi\) is still needed for the childlessness
moment. Cost low; blocks recalibration.

### P4. Two child states, \(n\) and \(m\)

**(a) At stake.** \(n\) enters choices only through \(\kappa(n)\) and
\(\xi\mathbf 1\{n=0\}\); \(m\) enters utility, the floor and the scale. \(n\)
is also the measurement state for completed fertility, childlessness and the
one-child share.

**(b) What breaks if they collapse.** Under memoryless maturation 22 percent of
first children leave within four years, so a mother at 30 with \(m=0\) and
\(n=1\) is common; collapsing would give her the first-birth scale and cost
again and would remove the three ever-born moments from the model. The cost of
keeping both is state space only: \(n\in\{0,1,2,3+\}\) with \(m\le n\), ten
cells.

**(c) Minimal change.** Keep both. Drop \(n\) from the state after the last
fertile age, since nothing after 46 depends on it once the bequest is
child-blind, if the solver does not already do so.

**(d) Observable.** None needed.

**(e) Recommendation.** Keep. Cost none.

### P5. Bequests

**(a) At stake.** The bequest object decides whether the old hold the large
homes for their own use or as the cheapest estate, hence how a holding tax
moves them, and where housing wealth goes at death.

**(b) Defensible?** The live form is child-blind De Nardi:
\(\theta_0(\theta_1+w^e)^{1-\sigma}/(1-\sigma)\) with \(\theta_n=0\) pinned
(`e5_profile.py:156`; form at `solver.py:7193-7227`). Child-blind is the
convention (De Nardi 2004; Lockwood 2018; Kopczuk and Lupton 2007 find an
extensive margin for having a motive, not a gradient in the number of
children), and the July 14 memo made the case. Two problems are in the
implementation, not the form. (i) \(w^e\) values the house gross of the 6
percent selling cost (`solver.py:2455-2459` against `solver.py:2442`), so
dying in the house is cheaper than selling it, an artificial reason to hold a
large home to death. (ii) The estate has no receiver. It is a utility argument;
the wealth leaves the economy, and the bequest flow is computed only as a
moment (`solver.py:6001-6062`). Greaney, Parkhomenko and Van Nieuwerburgh
make the same choice: households are born with zero wealth and bequests
leave the economy (their section 2.1.9, verified). For a paper about
intergenerational allocation it is a strange
choice: inheritances at ages 45–65 are how large homes and housing wealth
pass between generations in the data, and a wealth-to-earnings ratio of 4.85
against 6.15 with \(\beta\) at its cap is what a model without inheritances
looks like.

The bigger omission for the policy is not the bequest form but the absence of
every downsizing force the data show: widowhood, long-term care, health.
Venti and Wise (2004) document that the elderly rarely reduce housing equity
except after a spouse's death or a move to care (verify before citing). In the
model the old release large homes only at death, so a holding tax can move
them only along the intensive margin at a 6 percent cost.

**(c) Minimal change.** Net the selling cost in \(w^e\). Route estates to
households as a lump sum by age (uniform over 45–65 or proportional to
wealth), financed by the estate flow; no new state. Author decision on a
single age-rising forced-move shock calibrated to the mobility of the old
(ACS mover rates by age), which is the cheapest downsizing force; De Nardi,
French and Jones (2010) medical-expense risk is the full one.

**(d) Observable.** Ownership and rooms by age above 70, which in the model do
not fall; the age profile of sales; the housing share of estates.

**(e) Recommendation.** Change the two implementation points; author decision
on the downsizing shock. Cost medium; blocks recalibration.

## Fertility and children

### F1. Timing within the period

**(a) At stake.** Whether a household can adjust space in the same four years
as a birth, and whether the down payment bites before or after the birth.

**(b) Defensible?** The code order is fertility shock, attempt, conception,
housing shock, housing choice (`solver.py:2719-2768`, `kernels.py:672-756`;
the September 6 timing note). Housing is chosen knowing the realized family.
For a four-year period this is the right baseline: the PSID rooms response is
measured from one year before to three years after the birth, inside one
model period, so same-period adjustment is the object the target measures.
The reverse order (commit housing, then conceive) makes households buy the
family home before trying, the Hacamo pattern, and makes a birth in a small
home a real cost for up to four years. Both orders produce anticipatory
purchases across periods; they differ only within the period. The September 6
nested-choice experiment found the joint-versus-sequential difference at fixed
continuation values to be of order \(10^{-6}\) births per household; it did
not test commitment.

**(c) Minimal change.** None. One robustness run with housing committed before
conception bounds the timing sensitivity of the rooms response and of the
policy.

**(d) Observable.** Whether ownership rises before or after the first birth in
the PSID event study, once the measurement revision is settled.

**(e) Recommendation.** Keep and explain. Cost none.

### F2. Two extreme-value shocks

**(a) At stake.** The fertility shocks carry unobserved taste heterogeneity and
unplanned births and are identified by the level and age spread of first
births. The housing shock is over rent plus every owner size jointly
(`kernels.py:672-756`, `solver.py:2358`), a product-level taste shock rather
than a tenure shock.

**(b) Defensible?** No. \(\kappa_H=0.005\) is recorded as externally fixed in
the September 13 vintage (`corrected_initial/parameters.csv`), was the June 28
lower bound that searches returned (M19), and the package's own profile pins
0.0 (`e5_profile.py:156`). The fertility scales are 0.34 and 0.40 in the same
utility units, seventy times larger. At scales where the housing shock does
economic work, 0.05 and above, old-age ownership goes to 0.90 (June sweep) or
ownership collapses toward 42 percent (September 6 panel). So the data say the
housing choice is deterministic up to numerical noise, and the shock exists to
smooth the excess-demand map for the root finder; the September 7 failure at
the affordability threshold is the deterministic kink. Deterministic is also
the reference paper's convention: Greaney et al. put Gumbel shocks on
residence and workplace only and choose tenure and size deterministically
(their sections 2.1.4 and 2.1.6, verified), so a product-level taste shock
cannot be called "DUE-style". The problem is that the
device changes the answer: the August 20 vintage at \(\kappa_H=0\) gives a
positive ownership response to the tax, the September 1 vintage at 0.010 a
negative one (property-tax audit, confounded with closure and horizon).

**(c) Minimal change.** Two admissible options. Either \(\kappa_H=0\) with the
kink handled at the distribution level (interpolate wealth mass across the
affordability threshold rather than smoothing choices), or estimate
\(\kappa_H\) against a moment that measures the unpredictability of tenure
transitions: the PSID four-year-ahead ownership Brier score (0.117, SE 0.002)
is the right kind of moment and was dropped from the target ledger on July 24
(`docs/model/e5_target_review_20260724.md:74-76`). Either way the tenure-taste
isolation run (property-tax audit, item 2) must precede any reported policy
sign.

**(d) Observable.** The Brier score, or the share of tenure transitions
unexplained by wealth, income, age and children.

**(e) Recommendation.** Change. Cost low (declare zero) to medium (Brier
moment). Blocks recalibration.

### F3. Stochastic maturation, deterministic adult ageing

**(a) At stake.** How long a child occupies space, the age profile of parents'
space demand, and with survival the orphan flow.

**(b) Defensible?** The asymmetry is standard wherever child ages are not
tracked, and independence across siblings is harmless. Memorylessness is the
problem. With \(\mu=4/18\) per period (`parameters.py:888-921`;
`independent_count` runtime-asserted in `corrected_initial/candidate_result.json:1016`):
22 percent of newborns leave within four years, 37 percent of children are
still home after sixteen years, 8 percent after forty. The first distorts the
margin the paper is about, since a young parent's expected space need in the
birth period is cut by a fifth. The last is the orphan flow: every dependent
lost at parental death belongs to a parent aged 66 or more (M13), that is, a
child who should have left. An Erlang-2 law with the same mean cuts newborn
exits to 7 percent but leaves 25 percent at home after 24 years, so a
two-stage clock does not fix the tail.

**(c) Minimal change.** No new state: let the exit hazard depend on the
parent's age, \(\mu(a)\), low through the fertile window and rising to one by
62. The last birth is at 46 and mortality starts at 66 (SSA 2023 schedule,
`run_e1_chain.py:364-370`), so every child leaves before any parent can die.
Stochastic maturation for young parents, which the author wanted to keep, is
preserved. Add a one-bit newborn flag, since births are at most one per period,
that exempts the newborn from the first draw. Calibrate \(\mu(a)\) to the ACS
profile of own children at home by parent age.

**(d) Observable.** Dependents at home by parent age (ACS), which the model
overstates above 60 and understates at 25–35.

**(e) Recommendation.** Change. Cost low to medium; blocks recalibration (the
rooms-by-children moment is defined on \(m\)).

### F4. Dependents of a dying household

**(a) At stake.** About 1.2 percent of dependents per period, 5 percent of
births, vanish (`run_e5f_open_population_transition.py:847-851, 887`). The
accounting error is small; the conceptual error is that person counts and
household counts are not linked, so population statements are not statements
about persons.

**(b) Defensible?** Under F3's fix the flow is exactly zero: no household with
dependents dies. That is the least damaging option because it corrects the
assumption that created the flow. Joint death makes 5 percent of births die at
ages 20 to 40, which is not a demographic fact; a care pool adds a fiscal
object with no counterpart in the data; reassignment needs a matching rule.
All three treat a symptom.

**(c) Minimal change.** F3.

**(d) Observable.** The orphan flow itself, which should be zero.

**(e) Recommendation.** Change, through F3.

### F5. Exogenous fecundity, chosen attempts

**(a) At stake.** The split separates biology from choice and keeps
\(\kappa_1\) identified, since both smooth hazards.

**(b) Defensible?** Yes, and standard. Sommer (2016) uses an age-dependent
infertility shock, absorbing once it arrives (her section 3.4, verified);
Doepke and Kindermann (2019) make the per-period birth probability a
function that includes natural fecundity and imperfect birth control (their
section III.B, verified); de la Croix and Pommeret (2021) fit a
conception-by-age schedule. The live schedule
\(\pi_a=1-0.02e^{0.134(a-18)}\) (`run_e1_chain.py:380`, `e1_profile.py:34`;
`parameters.py:791-813`) gives 0.97, 0.90, 0.83, 0.71, 0.50 at ages 22, 30,
34, 38, 42 and 0.25 at 45, against Léridon's four-year figures 0.91, 0.84,
0.64 at 30, 35, 40. It was chosen by eye; the July 20 note calls the
least-squares fit a half-day item. Two things are missing. Attempts are
costless, so the only reasons not to try are the taste shock and the costs of
success. And every birth is chosen, whereas about 45 percent of U.S.
pregnancies were unintended in 2011 (Finer and Zolna 2016). Unintended births
are concentrated among the young and the poor, which is where the model needs
early first births and gets them from \(\kappa_1\).

**(c) Minimal change.** An exogenous birth hazard \(q_a\) for non-attempters,
calibrated to the unintended share by age. No state, one external schedule.
It lowers the price elasticity of total births by roughly the intended share
and makes \(\kappa_1\) a taste parameter rather than a stand-in for
contraception failure.

**(d) Observable.** First-birth hazards by age and income among households
the model classifies as non-attempters; the NSFG unintended share by age.

**(e) Recommendation.** Keep the split; author decision on \(q_a\). Cost low.

## Housing and finance

### H1. The constraint set

**(a) At stake.** The down payment at purchase is the mechanism. The rest of
the constraint set decides how binding it is: amortization and payment
limits make owning a cash-flow burden, free equity extraction makes housing
wealth liquid.

**(b) Verified and assessed.** Down payment on pre-income wealth with
\((1-\phi)\), \(\phi=0.80\) (`solver.py:2917-2922`, `kernels.py:640-662`).
Collateral limit \(b'\ge-\phi Ph\) re-applied to stayers every period at no
cost (`kernels.py:979-985`), so any owner can cash out to the limit each
period. \(\lambda_d=0\) (`parameters.py:146`). One-period debt at the saving
rate, no amortization, no spread (`parameters.py:166`, `solver.py:3824`).
Three consequences. Housing equity is perfectly liquid for stayers, so an
owner has no buffer-stock reason to hold liquid wealth, and the old pay a
holding tax out of equity, which mutes Coven's cash-flow channel. No
amortization means a young owner never faces forced saving; a thirty-year
mortgage retires roughly a tenth of the principal in its first four years.
Interest-rate risk is correctly absent under a constant \(R_b\). Nothing
present is unneeded except the age taper (H2). The block as it stands is
Greaney et al.'s: loan-to-value at origination, additional debt only while
the ratio is below \(\phi\), no amortization or refinancing (their equations
2.2–2.3, verified). The two papers the policy is compared with are richer:
Kaplan, Mitman and Violante (2020) have long-term amortizing mortgages with an
origination cost, refinancing at that cost, and loan-to-value plus
payment-to-income limits at origination; Coven et al. have a
payment-to-income limit of 0.36 at origination and an amortization rate of
0.0173 (Appendix B). Both devices bear on how a holding tax squeezes
cash-flow-constrained buyers and old owners, the margin the comparison is
about.

**(c) Minimal change.** A required amortization share per period, geometric so
that debt is zero at retirement; the collateral test at origination only
(purchase, or explicit refinancing at a cost); no cash-out for stayers; a
payment-to-income test at origination if the Coven comparison is to be
like-for-like. This is the Kaplan–Mitman–Violante contract in reduced form.
It removes the rollover taper: a stayer's balance follows the schedule
regardless of price.

**(d) Observable.** Liquid wealth of young owners (the July young-liquid-wealth
tension), loan-to-value by age, the share of owners who extract equity.

**(e) Recommendation.** Change. Cost medium; blocks recalibration.

### H2. The underwater rollover

**(a) At stake.** With one-period debt and an every-period collateral test, a
price fall would force stayers to inject equity or sell. The rollover prevents
that (`kernels.py:982-984`); the taper forces deleveraging between 42 and 62
(`parameters.py:147-148, 643-647`); a sale with shortfall carries unsecured
debt into the renter state that decays under the same taper
(`solver.py:3471-3474`). One array serves both as the unsecured-credit scale
(inert at \(\lambda_d=0\)) and as the rollover share.

**(b) Defensible?** The economics is right: U.S. mortgages are not called after
a price fall. The taper is ad hoc, and unstated it is a referee's discovery.

**(c) Minimal change.** It disappears under H1. If H1 is not adopted, the
rollover belongs in the paper's appendix with the taper stated as an
assumption.

**(d) Observable.** None in the model either way.

**(e) Recommendation.** Change through H1; otherwise keep and disclose. Cost
low.

### H3. The rental cap

**(a) At stake.** The cap is the tenure segmentation: the owner menu is
\(\{2,4,6,8,9.5,11\}\) rooms, renters up to 6 (`parameters.py:187-191`), so
family sizes are owner-only. It is load-bearing for the ownership fit and for
the fertility lever (July 1 size-mapping audit).

**(b) Defensible?** The cap itself has the precedent the slide claims:
Greaney et al. cap rental size at 5.39 against a largest owner size of 10.78
and calibrate the level to the ratio of the 90th percentile of rental units
to the 10th percentile of owner-occupied units (their section 3.2, Table 3,
verified). So "6 rooms, Greaney et al." should become "the 90th percentile of
rented units, as in Greaney et al.", measured in the AHS. Beyond the level, a
hard cap is the limit of a size-dependent rental premium, and it fails at the
edge: 7.5 percent of renters live in seven-plus-room
units and renters hold about 8 percent of that stock (AHS 2023,
`code/data/ahs_supply_snapshot/output_ahs_family_unit_menu_national/ahs_stock_by_room_tenure.csv`).
More important, 31.8 percent of renters aged 25–45 sit exactly at the cap
(property-tax audit), a mass point that creates the kinks \(\kappa_H\)
smooths. The economic reason large units are owned is the rental externality
of Henderson and Ioannides (1983): tenants under-maintain and the cost rises
with the unit, so landlords price large units above user cost. That argument
delivers a wedge, not a wall. (Their paper is not in the local library; the
argument is cited from memory and should be checked before use.)

**(c) Minimal change.** Rent per room \(r_t(h^R)=r_t+w(h^R)\) with \(w\)
increasing, calibrated so that renter shares by size match the AHS bins (the
three bins in the snapshot pin a three-point schedule). Under outside-landlord
pricing \(w\) is a landlord operating cost, which is exactly how Kaplan,
Mitman and Violante price rentals: their rent is the user cost plus a
per-period operating cost \(w\) per unit rented out, from a competitive rental
sector that owns the units (their equation 11, verified). Renter rooms are
already continuous, so the code change is small; the cap becomes a
sensitivity.

**(d) Observable.** Renter share by unit size, which a cap matches only as a
step.

**(e) Recommendation.** Change. Cost medium; blocks recalibration.

### H4. Supply and rent

**(a) At stake.** Supply decides how a demand shift divides between prices and
quantities and what the long run looks like as population shrinks. The rent
rule decides who bears the anticipated capital loss.

**(b) Defensible?** No, on two counts. The stationary rule is
\(H^s=H_0(\text{user cost}/\bar r)^{\eta}\) (`solver.py:1765`) and the dated
rule is \(H^s_t=H^s_0(P_t/P_0)^{\eta}\) evaluated at the current price only
(`code/model/tools/run_dynamic_population_transition.py:125-141, 471`): no
lag, no stock state, no depreciation dynamics. Two elasticities coexist in the
code: 1.75 in the profile chain that builds the parameter object
(`e1_profile.py:39`, `run_e1_chain.py:390`) and 0.63 as a command-line
override on the dated rule, recorded as externally fixed
(`run_e5f_transition_calibration.py:2014-2030`;
`output/model/e5f_original_queue_20260913a/inherited_2023_tax/manifest.json`).
Nothing forces the two onto one number, so the value used by a stationary
anchor solve depends on whether the flag was passed; that must be checked
before any terminal state is reported. At impact the stock is reversible: the
tax's 17 percent price fall removes
\(0.63\times17\approx11\) percent of housing within one period, while the
physical stock can fall at most by depreciation, 8 percent in four years, and
only with zero construction. This is the "housing services fall" that the
September 1 note misread as a mechanism. In the long run a level supply curve
with a population going to zero sends \(P\to0\); there is no stationary
per-household equilibrium along a contracting population, which is why the
transition has no terminal steady state (M35; the August 15 finding that
\(B/E\) peaks at 0.87 even at a near-zero price). The number 0.63 has no
derivation from Baum-Snow and Han, whose headline is about 0.5 for floor
space (M29; `docs/model/e5f_independent_quantitative_audit.md:382`), and the
older 1.75 has no citation at all (`e1_profile.py:39`). The rent rule
\(r_t=(R_b+\delta_H+\tau^p_t)P_t-P_{t+1}\)
(`run_e5f_perfect_foresight_transition.py:233-256`; stationary
\(r=(q+\delta_H+\tau^p)P\), `parameters.py:255`) is the arbitrage condition of
a deep-pocketed competitive landlord under perfect foresight and is fine once
the landlord is named (H5). Landlords bearing no risk follows from perfect
foresight, not from the rule. One number on the slide is wrong: live
depreciation is \(1-0.989^4\), 1.1 percent a year (`run_e1_chain.py:389`),
not the 2 percent the parameter table states; the module default is 2
percent (`parameters.py:164`).

**(c) Minimal change.** Stock-flow supply,
\(H_{t+1}=(1-\delta_H)H_t+I_t\), \(I_t=I_0P_t^{\eta}\), with the elasticity
on construction. This is the reference paper's own block: Greaney et al. have
irreversible construction, construction demand equal to \(\dot H+\delta H\),
a supply elasticity on that flow, and steady-state construction that exactly
offsets depreciation (their sections 2.2 and 3.2, verified); Kaplan, Mitman
and Violante have a competitive construction sector with elasticity 1.5
(their equation 13, verified). Coven et al. use the same static form as this
model, \(H=cP^{\rho}\) with \(\rho=0.232\) for California (their equation 5,
verified), so for the Coven comparison the static rule is like-for-like; the
objections above stand on their own. Short-run supply becomes
inelastic, so the price fall at impact is larger and the transfer to young
buyers stronger, which moves the model toward Coven's first stage. Along a
contracting population the stock shrinks through depreciation and a
stationary per-household equilibrium exists whenever the population's
per-period decline is slower than depreciation, which holds at any realistic
completed fertility. This is the change that makes D2 solvable.

**(d) Observable.** The four-year response of permits to prices versus the
response of the stock; the price–rent ratio along the transition.

**(e) Recommendation.** Change. Cost medium (one aggregate state, a new
terminal condition). Blocks the transition and \(H_0\), not the household
block.

### H5. The landlord

**(a) At stake.** Rents go to an unmodeled agent: no landlord sector or
rental-income identity exists in the package (zero hits for a landlord or
resource constraint across the live code). The tax is levied on rented and
owner-occupied stock alike and rebated per household
(`solver.py:7356-7371, 5054-5065`), so renters pay it in full through the user
cost and get their share back. The welfare accounting depends on who the
unnamed agent is.

**(b) Defensible?** The consistent reading is a small open economy: \(R_b\) is
exogenous (live 2 percent a year, `run_e1_chain.py:388`; the module default
is 4 percent, `parameters.py:163`), there is no bond-market clearing
condition anywhere in the live code, the landlord is a foreign
investor who earns \(R_b\) and bears the capital loss on the rented stock when
the tax lowers \(P\). That is coherent and is what the exogenous rate already
implies. It is also Greaney et al.'s construction: rented units are owned by
perfectly competitive real estate investment trusts and rents satisfy
\(r=(\delta+q+\tau^h)p-\dot p\), the continuous-time twin of the rule here
(their section 2.3, equation 2.27, verified). In the June 2025 draft of
Coven et al. held locally there is no landlord agent either: rent is a
user-cost no-arbitrage price with the tax inside it, so the tax passes to
renters there exactly as here (their equations 7 and 16, verified). The
September 1 note attributes a different treatment to the August 2026 version
(owner-occupied tax base, rental stock clearing separately, no
pass-through); that version is not in the library, so whether rental
incidence differs from Coven is unverified either way and must be checked
against the current version before it is stated. The open-economy reading
matters for welfare twice: rents
are a leakage abroad, and part of the policy's price fall is paid by
foreigners, a free lunch in domestic accounting. It matters for the rebate
only through the base.

**(c) Minimal change.** Name the trust, as Greaney et al. do. Report the
trust's capital loss and the rent leakage as lines in the policy accounting.
Domestic ownership of the rental stock through a fund whose returns are
rebated in proportion to wealth is the full fix and is only needed if welfare
numbers become headline results.

**(d) Observable.** None separates the readings; the choice is stated.

**(e) Recommendation.** Keep and state. Cost low.

### H6. The owner premium \(\chi\)

**(a) At stake.** In the model renting and owning the same size cost the same
user cost. Without \(\chi\) no one would own, given the sale cost and the down
payment. \(\chi=1.05\), applied as \(\chi(h-\bar h(m))\)
(`kernels.py:957-970`; the slide writes \(\chi h\) with the floor subtracted
after), is therefore the only reason for ownership below the cap and stands in
for every missing owner advantage: non-taxation of imputed rent and the
mortgage interest deduction (Gervais 2002; Sommer and Sullivan 2018), the
rental externality (H3), the hedge against rent risk (Sinai and Souleles
2005), which is absent under perfect foresight.

**(b) Defensible?** A premium scales services; a rent wedge raises the price
of renting. The ownership rate identifies their sum, and only the size
dimension separates them. Since H3 introduces a size-dependent wedge anyway,
one device suffices: \(w(h^R)\) with a positive intercept delivers the level
of ownership and its size gradient with \(\chi=1\).

**(c) Minimal change.** Replace \(\chi\) by the intercept of the rental wedge;
calibrate the intercept to the ownership rate and the slope to renter shares
by size. Do it jointly with H3.

**(d) Observable.** The price–rent ratio, which a wedge raises and a premium
does not, and ownership by size.

**(e) Recommendation.** Author decision, lean change with H3. Cost medium;
blocks recalibration.

## Demography, government, closure

### D1. Household, person and dependent accounting

**(a) At stake.** Every population object depends on what a household is.

**(b) The accounting as implemented.** The implicit unit is a couple: earnings
are PSID head-plus-spouse (`e6b_profile.py:75-78`), the scale base is two
adults, children ever born per household are matched to CPS children per
woman, and entrants are births divided by 2.1
(`output/model/e5f_original_queue_20260913a/spec.json`, conversion 0.476),
that is, two adults per new household plus a 5 percent leakage. Children ever
born are literal with a top-coded three-plus bin adjusted to the CPS mean
above three (`run_e5f_open_population_transition.py:740-780`). Births enter a
queue and become households five periods later, twenty years rather than
eighteen (`:230-238`), as childless renters with PSID 18–24 wealth ratios
(`calibration.py:60-97`, `run_dynamic_population_transition.py:345-363`),
with the earnings state drawn from its ergodic distribution. Survival is the
SSA 2023 schedule from 66 and death is certain at 82
(`run_e1_chain.py:364-370`, `solver.py:2486-2487`). A death removes the
household, its dependents (F4) and its estate (P5).

What the model gets wrong: (i) the couple never dissolves, so widowhood, which
is how most old households become one-person households in large homes, is
absent; (ii) the 5 percent leakage inside 2.1 is an external constant while the
model's own leakage, the orphan flow, is another 5 percent of births, similar
by accident; (iii) the 2.1 normalization makes the 2007 economy a stationary
population by construction (row 1 of `corrected_initial/target_fit.csv` is a
normalization, not a moment), while childlessness and the one-child share are
for the 1960–66 birth cohort and the timing moments are 2003–06 period
statistics. In a stationary model cohort and period coincide; in the data they
did not.

**(c) Minimal change.** One table: births to queue to entrants; survivors to
maturation; deaths to estates; the unit conventions; (i)–(iii) stated as
assumptions. Widowhood is the one structural addition worth its cost for the
intergenerational figure and can be an income and scale shift at an exogenous
age-specific rate without a new state.

**(d) Observable.** One-person households among owners aged 75 and over, and
their rooms.

**(e) Recommendation.** Keep and write; author decision on widowhood. Cost low.

### D2. Population closure

**(a) At stake.** Whether the model can make a population statement at all and
whether the transition has an end.

**(b) Defensible?** Two branches coexist and the prompt describes the wrong
one. The retained September 13 branch behind the slides propagates households
endogenously with the births/2.1 queue, with no imposed age masses, no
headship law and no migration (`spec.json`: `historical_age_conditioning:
false`, `person_headship_transition: false`; `history_manifest.json:694`).
The headship-on-persons law exists only in the person-demography branch
(`build_e5f_coherent_person_cohort_path.py:153-169`,
`household_person_coupling.py:33-99`). In either branch the post-2023 law is
common to baseline and policy, so the counterfactual is internally coherent:
the policy changes births, births change entrants twenty years later. The
incoherence is elsewhere. With zero migration and completed fertility below
2.1 no positive stationary population exists (\(B/E\le0.87\) at any price,
August 15), so the terminal steady state the transition solver uses is a
fiction, and with a level supply curve the shrinking economy has no
stationary per-household limit either (H4). Where the headship law is used it
freezes household formation at 2023 rates, so housing costs cannot move the
formation margin that the data show (Ermisch 1999; Couillard's
living-arrangement margin).

**(c) Minimal change.** Two coherent options. Option 1: the August 6 closure
(iv), a stationary level equilibrium with a policy-invariant outside inflow
anchored to net international migration; population responds through births
with the renewal multiplier. Option 2: zero migration, the author's stated
preference, with a balanced-contraction terminal condition: population
shrinks at the constant rate implied by terminal completed fertility while
per-household quantities and prices are stationary, which requires H4's
stock-flow supply. Option 2 is the only way "zero migration" becomes a
well-posed forecast. Both keep the 2007 stationary normalization.

**(d) Observable.** None for the closure itself. For the formation margin,
headship rates by age against rents across metros.

**(e) Recommendation.** Author decision, two positions. Cost low for option 1,
high for option 2 (new terminal solver). Neither blocks the initial
calibration; both block transition and policy.

### D3. PAYGO pensions and the rebate

**(a) At stake.** The pension rule sets old-age income along an ageing
transition; the rebate rule sets who gets the tax back.

**(b) Defensible?** A fixed payroll tax of 17.9 percent with the benefit
adjusting each period (`parameters.py:256-257, 741-751`; dated joint root of
price, pension and rebate in `e5f_rebated_surprises.py:274`) is one of the two
standard closures.
The other, a fixed replacement rate with the tax adjusting, would put the
burden of ageing on the young, the constrained group, and strengthen the
mechanism in the long run; fixed tax is the conservative choice and matches
U.S. practice up to the trust fund. The equal lump-sum rebate to every
household is Coven's rule (all property-tax revenue rebated to residents as an
equal lump-sum transfer, their section 3.6, verified) and enters the period
budget after the
housing choice, not the down-payment test (`solver.py:226-233`,
`kernels.py:640-662`), so it helps buyers only through accumulation. Renters
pay the tax through rent and get the rebate back; owners of above-average
housing pay net. Nothing is wrong. The interaction to report: in a shrinking
population pensions fall, the old get poorer, and that is a downsizing force
the model otherwise lacks (Mankiw and Weil 1989 is the reference for
demographic housing demand).

**(c) Minimal change.** None. Report the fiscal lines by age and tenure along
the path.

**(d) Observable.** The pension replacement rate along the transition.

**(e) Recommendation.** Keep. Cost none.

### D4. Geography

**(a) At stake.** \(H_0\), \(\chi\), \(h_P\), the rooms moments and the
ownership rate are 42-metro objects (`working_contract.json`, geography
fields); fertility, wealth and the rooms response are national; the
population closure is national. A pooled market calibrated to two populations
has parameters that belong to neither.

**(b) Defensible?** No, and the author agrees (M33). The mechanism is general;
the property tax is state or local; the population statement is national. The
national economy is the only one for which every block has data: ACS national
rooms and ownership exist, the PSID and CPS are national, Census population is
national. A metro economy would need metro fertility and wealth moments that
the CPS and PSID cannot deliver for 42 metros with useful precision.

**(c) Minimal change.** Remeasure the four ACS rows nationally under the same
sample rules; keep the 42-metro versions as a sensitivity that shows how much
of the mechanism is an expensive-metro phenomenon. \(H_0\) re-estimates;
\(\chi\) and \(h_P\) move. No code change.

**(d) Observable.** The national-versus-metro gap in the four rows.

**(e) Recommendation.** Change to national; author decision if the paper is to
be about expensive metros, in which case the fertility moments must become
metro too. Cost medium (data); blocks recalibration.

## Earnings

### E1. Three permanent types times a persistent AR(1)

**(a) At stake.** The mechanism is a liquidity constraint on the young.
Permanent heterogeneity decides who is constrained forever and who never is,
and it is the only source of old-age wealth dispersion in a model with no
return heterogeneity.

**(b) History and assessment.** The three types were added on July 27 (E6b)
to reach the old-age wealth tail (p90/p50 from 2.07 to 3.91 at fixed
parameters) after the AR(1) alone could not. They come from a PSID
minimum-distance decomposition of head-plus-spouse gross earnings, 1984–2019:
fixed 0.393, AR(1) with annual \(\rho=0.886\) and stationary variance 0.332,
transitory 0.310 (`e6b_profile.py:13-39, 75-78`). The AR(1) in the live model
is not that one. It is Floden and Lindé (2001), annual \(\rho=0.9136\),
innovation variance 0.0426, with the innovation scaled by \(1-0.181\) for
Heathcote–Storesletten–Violante tax progressivity (`externals.py:12-25`),
converted to four years as \(\rho_4=0.697\), \(\sigma_{\varepsilon,4}=0.298\),
stationary log-variance 0.173 (`local_panel.py:1063-1071`), on five
Rouwenhorst points crossed with the three types (`e6b_profile.py:50-54`). So
the live process is an after-tax-scaled literature AR(1) plus a gross-earnings
PSID fixed effect. This is not double counting. It is an incoherent sum:
persistent-plus-fixed variance 0.566 in the model against 0.725 in the
decomposition the fixed effect came from, in different tax units (the fixed
effect after the same HSV scaling would be 0.264, not 0.393). What the types
buy: the tail, the wealth-to-earnings level, and a childlessness gradient
across types that E6b reversed relative to the CPS (July 27), which is the
symptom of P3's missing income-dependent child cost. What the paper loses
without them: the p90/p50 target and much of wealth over earnings, unless the
AR(1) is made far more persistent (annual \(\rho\) near 0.97), which is the
RIP representation of the same cross-sectional dispersion (Guvenen 2009;
Storesletten, Telmer and Yaron 2004). The age profile \(e_a\) has no source
(`parameters.py:195-196`).

**(c) Minimal change.** One decomposition, one tax treatment: use the E6b PSID
estimates for all three components and scale persistent and fixed by the same
after-tax factor, or model the tax explicitly. Parameter-file change. Then run
the author's comparison: types removed at fixed parameters, then refit under
the same contract, reporting wealth dispersion, ownership by income and the
policy response.

**(d) Observable.** The age profile of the cross-sectional variance of log
earnings: a fixed effect is a level at 25, an AR(1) a rising profile. And the
old-age wealth tail.

**(e) Recommendation.** Change (make the process coherent); author decision on
removal only after the comparison. Cost low; blocks recalibration.

### E2. Earnings risk and fertility

**(a) At stake.** The precautionary channel of fertility (Sommer 2016) is the
option value of postponing an irreversible commitment under income risk.

**(b) Defensible?** The model has it: children are irreversible, the scale
raises marginal utility for the whole dependency spell, and the value functions
integrate over \(z\). Nothing should enter the fertility decision "beyond
wealth", because that is the channel. What mutes it is frequency: at four
years the AR(1) has persistence 0.70, so most risk is near-transitory at the
decision horizon, and there is no unemployment or disaster state.

**(c) Minimal change.** None. Report the fertility response to a rise in
\(\sigma_\varepsilon\) at fixed parameters as a diagnostic and compare it with
Sommer's.

**(d) Observable.** Timing of first births by earnings volatility (occupation
or industry) in the PSID.

**(e) Recommendation.** Keep. Cost none.

## Cross-cutting

### X1. Interactions

- **H4 with D2.** Without stock-flow supply there is no balanced-contraction
  long run; choose them together.
- **F3 with F4 and P4.** One maturation change closes the orphan flow and
  restores the meaning of the \(n\)–\(m\) distinction.
- **F2 with H3 and H6.** The housing shock exists to smooth kinks created by
  the cap and the discrete owner menu; a size wedge removes the mass point at
  the cap; the wedge and \(\chi\) are one device.
- **H1 with H2 and P5.** Amortization removes the taper and changes what the
  estate contains; the estate must then be valued net of debt and sale cost.
- **P1 with P3 and E1.** The income gradient of fertility is where the linear
  child term, the missing child cost and the permanent types meet; E6b's
  reversed childlessness gradient is their joint symptom.
- **D4 with everything estimated.** A national remeasurement moves \(H_0\),
  \(\chi\), \(h_P\) and the weights; do it before, not after, the refit.
- **P5 with \(\beta\) at its cap.** Estates delivered to households raise
  mid-life wealth and relieve the patience the wealth target now demands.

### X2. Smallest defensible set

Change, seven items, in this order: F3/F4 maturation and E1 earnings
coherence first (parameter-level, cheap, must precede any refit); H1/H2
mortgage block and F2 housing shock second (household block, one refit); H4
supply with D2's terminal condition third (transition); D4 geography with the
refit; P5's two implementation fixes with it.

Explain better, no change: P2, P4, F1, F5, H5, D1, D3, E2.

Author decisions with two positions: P1 (linear versus concave child term),
P3 (\(\xi\) versus a measured child penalty), H3 with H6 (cap and premium
versus one rental wedge), D2 (outside inflow versus balanced contraction),
F5's unintended births, P5's downsizing shock.

A hostile referee's first four questions will be H4, E1, F2 and X3's item 1.
None is answered by better exposition.

### X3. Not on the list

1. **No time cost of children.** The cost of a child is a homothetic scale
   plus an absolute space floor, so relative to income it falls with income,
   and the model cannot produce the negative income–fertility gradient (Becker
   and Lewis 1973; Jones and Tertilt 2008). Both quantitative fertility
   references carry the time cost: Sommer's parents give up labor time to
   child quality with a minimum quality floor, and Doepke and Kindermann's
   pay the mother's foregone wage or the childcare price (verified, Appendix
   B). E6b's gradient reversal is the symptom. This is the largest omission
   for a fertility paper and is P3's replacement.
2. **No household-formation margin.** Entrants appear as households at
   twenty with PSID wealth; housing costs cannot delay formation, which the
   data show (Ermisch 1999; Couillard 2025 models the living-arrangement
   margin explicitly). Where the headship law is used it freezes the margin at
   2023 rates.
3. **No widowhood and no health-driven downsizing** (P5, D1).
4. **The taste shocks must be stated as mean-zero type-I extreme value.** The
   code computes \(\kappa\log\sum\exp\) with no Euler constant
   (`utils.py:194-198`), which is exact under that convention (M11 agrees);
   under location-zero shocks the value of every choice occasion is
   understated by \(\kappa\gamma\), which is not innocuous across states with
   different menus.
5. **Inconsistent period conversions.** Depreciation compounds,
   \(1-0.98^4\); the property tax is \(4\times0.01\)
   (`parameters.py:163-165`). Harmless numerically, ugly to a referee.
6. **\(\beta R>1\) at the cap** (0.99 by 1.02 annually). The model needs
   implausible patience to hold six years of earnings because the old have no
   expense risk, no return above \(R_b\) and no inheritances. That is
   structural, not a search failure, and it says the wealth target is being
   hit by the wrong margin.
7. **Owner services** are \(\chi(h-\bar h)\) in code and \(\chi h-\bar h\) on
   the slide (`kernels.py:957-970`). State one.
8. **\(h_P\) at its upper bound with the per-child slope restricted to
   zero.** The model wants more space per parent than the bound allows while
   overshooting the first-birth rooms response, so the two rooms moments pull
   \(h_P\) in opposite directions and the three-plus versus one-to-two gap has
   no lever. That is the identification statement M25 asks for.
9. **One price per room at every size.** The down payment on an 11-room home
   is 5.5 times that on a 2-room home; hedonic prices per room fall with size.
   A concave \(P(h)\) schedule is cheap and measurable.
10. **Fecundity is live but undocumented in the calibration layer.** The
    schedule is set in `run_e1_chain.py:380` and `e1_profile.py:34`, not in
    `local_panel.py` or the profiles; a reader of the calibration layer
    concludes conception is certain. Move it to the profile that owns the
    other externals.

## Final table

| Label | Recommendation | Cost | Blocks recalibration |
|---|---|---|---|
| P1 | Decide | Low | Yes |
| P2 | Keep | None | No |
| P3 | Decide (lean change) | Low | Yes |
| P4 | Keep | None | No |
| P5 | Change | Medium | Yes |
| F1 | Keep | None | No |
| F2 | Change | Low–medium | Yes |
| F3 | Change | Low–medium | Yes |
| F4 | Change (via F3) | None extra | Yes |
| F5 | Keep (decide \(q_a\)) | Low | No |
| H1 | Change | Medium | Yes |
| H2 | Change (via H1) | Low | Yes |
| H3 | Change | Medium | Yes |
| H4 | Change | Medium | Transition only |
| H5 | Keep | Low | No |
| H6 | Decide (lean change) | Medium | Yes |
| D1 | Keep | Low | No |
| D2 | Decide | Low / High | Transition only |
| D3 | Keep | None | No |
| D4 | Change | Medium | Yes |
| E1 | Change | Low | Yes |
| E2 | Keep | None | No |

## Appendix A. The prompt's code claims, checked

Paths are relative to `code/model/intergen_eqscale_seq_optimized/` unless
another root is given.

1. **No unsecured borrowing.** Confirmed: `parameters.py:146` `lambda_d = 0.0`;
   `parameters.py:677` builds zero caps at every age.
2. **Underwater rollover with a taper from 42 to 62.** Confirmed:
   `parameters.py:147-148, 643-647`; `kernels.py:982-984`; the shortfall on
   sale carries into the renter state, `solver.py:3471-3474`. Correction: the
   taper array also scales new unsecured credit, inert at \(\lambda_d=0\).
3. **\(\kappa_H\) estimated at its lower bound 0.005.** Partly wrong. The
   September 13 vintage records it as externally fixed at 0.005
   (`corrected_initial/parameters.csv`); the package profile pins 0.0
   (`e5_profile.py:156`); the module default is 0.01 (`parameters.py:128`);
   the bound declarations disagree (`production_profile.py:49` gives
   \([0,0.12]\)). The value 0.005 is the June 28 lower bound that searches
   returned (M19). It has never had a moment; the Brier candidate was dropped
   on July 24.
4. **\(\xi\) subtracted once inside the attempt value at \(n:0\to1\).**
   Confirmed: `solver.py:2718-2729`, multiplied by \(\pi_a\); default zero
   (`parameters.py:83`); live 0.266.
5. **\(\psi_0\) is a nested root, not an SMM coordinate.** Confirmed:
   `code/model/tools/run_e5f_transition_calibration.py:609-712`
   (bracket plus false position), hard-gated to 2.1 at `:1670-1675` and
   `run_e5f_open_population_transition.py:48`; `parameters.csv` row
   `psi_child` marked "normalized to 2.1". The free SMM set is nine
   parameters: \(H_0,\beta,\chi,\xi,h_P,\kappa_1,\kappa_C,\theta_0,\theta_1\),
   with \(\beta\), \(\kappa_1\), \(\kappa_C\), \(\theta_1\) and \(h_P\)
   flagged at or near bounds. The slide's "ten" counts \(\psi_0\).
6. **Earnings: five-point Rouwenhorst AR(1) times three permanent types.**
   Confirmed: `externals.py:12-25` (Floden–Lindé, HSV scaling),
   `local_panel.py:1063-1071` (four-year conversion), `e6b_profile.py:13-54`
   (types \([0.278,0.822,2.435]\), weights \([1/6,2/3,1/6]\), PSID variance
   0.393). Module default (three points, persistence 0.85,
   `parameters.py:198-202`) is overridden. Types are drawn once and fixed;
   entry wealth depends on the full earnings state; fecundity, survival and
   \(\psi\) do not.
7. **Population: exogenous survival, no migration, imposed historical
   masses, headship after 2023, dependents lost.** Half wrong. Survival and
   the orphan loss are confirmed (`run_e1_chain.py:364-370`,
   `run_e5f_open_population_transition.py:847-851, 887`). The retained branch
   behind the September 14 figures imposes no age masses and uses no headship
   law: `output/model/e5f_original_queue_20260913a/spec.json` has
   `population_law: original_household_birth_vintage_queue`,
   `historical_age_conditioning: false`, `person_headship_transition: false`.
   Imposed masses and headship belong to the person-demography branch
   (`build_e5f_coherent_person_cohort_path.py:153-169`). Entrants arrive
   twenty years after birth (four waiting slots,
   `run_e5f_open_population_transition.py:230-238`), as childless renters with
   PSID 18–24 wealth ratios (`calibration.py:60-97`).
8. **Targets: 42 metros for housing, national otherwise.** Confirmed in
   `output/model/e5f_matched_pf_20260909a/design_research/working_contract.json`
   (geography fields) and `observer_contract/README.md:7`; 13 live rows in
   `corrected_initial/target_fit.csv`, row 1 a normalization, row 13 the
   recent-parent ownership gap the slides omit.
9. **Supply \(H_0P^\eta\), \(\eta=0.63\), no lag.** Confirmed for the dated
   rule (`run_dynamic_population_transition.py:125-141`, current price only)
   and the stationary rule (`solver.py:1765`). Correction: 0.63 is a
   command-line override on the dated rule; the profile chain that builds the
   parameter object still carries 1.75 (`e1_profile.py:39`,
   `run_e1_chain.py:390`). No derivation of 0.63 exists on file (M29); 1.75 is
   uncited. Live depreciation is 1.1 percent a year (`run_e1_chain.py:389`),
   not the slide's 2 percent.
10. **Selling cost 6 percent and rental cap 6 rooms cited to Greaney et al.**
    Values confirmed (`parameters.py:168, 191`); no citation in code. The AHS
    snapshot bins are 1–4, 5–6, 7+ rooms; renters hold about 8 percent of the
    7+ stock. See Appendix B for what Greaney et al. actually assume.
11. **Property tax 1 percent rebated equally; experiment 2 percent.**
    Confirmed: revenue on rented plus owned stock, equal transfer per
    household (`solver.py:7356-7371, 5054-5065`); the transfer enters income
    for every household (`solver.py:226-233`); the experiment's manifest
    records annual taxes \([0.01,0.02]\), equal rebate, balanced PAYGO at
    fixed payroll 0.179 and supply elasticity 0.63, and the driver asserts
    that only the tax rate differs between arms
    (`code/cluster/run_e5f_inherited_2023_tax_long.py:112-120`).
12. **\(\beta\) at the 0.99 cap.** Confirmed: raw bound 0.9995
    (`e5_profile.py:145`), cap enforced in
    `code/model/tools/run_e5f_rebated_initial_overnight.py:119-122`, saved
    value 0.99.

Further corrections to the model as presented: fecundity is the age-declining
schedule, set in `run_e1_chain.py:380` (one verifier reading only the
calibration layer concluded conception is certain; that reading is wrong);
maturation is the binomial `independent_count` law by runtime assertion
(`corrected_initial/candidate_result.json:1016`), not the module-default
clock; owner services net the floor before the premium; the estate is gross of
the selling cost and has no receiver; the down-payment test uses pre-income
wealth and the rebate enters after the housing choice.

## Appendix B. Literature checks

Each claim was checked against the paper's own text in the local library by a
read-only verifier; "verified" means a verbatim passage was found. Papers not
on disk are marked and are cited above from memory with that caveat.

**Greaney, Parkhomenko and Van Nieuwerburgh, Dynamic Urban Economics (2025;
Zotero MAXJ699L).** Rental size is capped: largest rental unit 5.39 against
largest owner unit 10.78, calibrated to the ratio of the 90th percentile of
rental units to the 10th percentile of owner-occupied units (section 3.2,
Table 3). Transaction fee \(\psi=0.06\), "standard value" (section 2.1.5).
Households are born with zero wealth; bequests leave the economy (section
2.1.9). Construction is irreversible, construction demand is \(\dot H+\delta
H\), a supply elasticity applies to it, and steady-state construction offsets
depreciation (sections 2.2, 3.2). Rented units are owned by competitive REITs
and \(r=(\delta+q+\tau^h)p-\dot p\) (section 2.3, equation 2.27).
Loan-to-value at origination with additional debt only while the ratio is
below \(\phi\); no amortization or refinancing is described (equations
2.2–2.3). Gumbel shocks are on residence and workplace, scales 2.08 and
0.22; tenure and size are chosen deterministically (sections 2.1.4, 2.1.6).
Interest rate \(q=0.02\), the 1962–2024 average ten-year real rate (section
3.2).

**Kaplan, Mitman and Violante (2020, JPE; Zotero GBH534TV).** Rent equals the
user cost plus a per-period operating cost \(w\) per unit, from a competitive
rental sector that owns the units and pays property taxes (equation 11).
Mortgages are long-term, carry a fixed origination cost, amortize over the
buyer's remaining life, can be refinanced at the origination cost, and face
loan-to-value and payment-to-income limits at origination (equations 5–7).
A competitive construction sector gives an aggregate supply elasticity
\(a/(1-a)\), calibrated to 1.5 (equation 13). The shock is an aggregate
housing-preference shock; no idiosyncratic moving shock.

**Doepke and Kindermann (2019, AER; Zotero FXNH3B64).** Utility is linear in
consumption and additively separable in a one-time felicity \(v_g b\) from a
birth (section III). Each child carries a fixed monetary cost \(\phi_c\), a
fixed utility cost \(\phi_u\), and a time cost paid through the mother's
foregone wage or the childcare price. The per-period birth probability
\(\gamma(\cdot)\) reflects natural fecundity and imperfect birth control
(section III.B).

**Sommer (2016, JME; Zotero KJD3R4Z4).** Utility
\(c^{1-\gamma}/(1-\gamma)+\zeta(nq)^{1-\kappa}/(1-\kappa)\), concave in
children times quality (p. 33). Quality is produced with parental time and
goods, with a minimum quality for families with children; labor supply is
time not spent on child rearing (section 3.2). An age-dependent infertility
shock arrives each period and is absorbing (section 3.4). Headline: with
children as consumption commitments, higher earnings risk lowers family size
and delays childbearing (abstract).

**Coven, Golder, Gupta and Ndiaye (2025; `docs/reference/coven2025_property_tax.txt`,
the June 19, 2025 draft; the August 2026 version cited by the September 1
memory note is not on disk and may differ on the rental block).**
Rents follow a user-cost no-arbitrage relation \(R_i=(\tau_i+r+\delta-\gamma)P_i\)
with no landlord agent (equations 7, 16). All property-tax revenue is rebated
to residents as an equal lump-sum transfer (section 3.6). Supply is
\(H_i=c_iP_i^{\rho_i}\) with \(\rho=0.232\) for California from Baum-Snow and
Han (equation 5, Table 4). Lock-in devices: a 5 percent transaction cost, a
15 percent capital-gains tax, and a years-owned state. A payment-to-income
limit of 0.36 at origination and an amortization rate of 0.0173. No fertility
or children anywhere in the state vector. Results: California prices fall 11.2
percent; ownership rises from 61 to 67 percent overall and from 35 to 43
percent at ages 25–44 (section 5; the introduction says nine points for the
young, an internal inconsistency in the paper).

**Not in the local library.** Sommer, Sullivan and Verbrugge (2013);
Henderson and Ioannides (1983). Claims attributed to them above are from
memory and are marked as such.

**Cited from memory, not checked here.** Scholz, Seshadri and Khitatrakun
(2006) (a local PDF exists under `docs/reference/`, not re-read); De Nardi
(2004); Lockwood (2018); Kopczuk and Lupton (2007); De Nardi, French and
Jones (2010); Venti and Wise (2004); Favilukis, Ludvigson and Van Nieuwerburgh
(2017); Floden and Lindé (2001); Heathcote, Storesletten and Violante (2017);
Guvenen (2009); Storesletten, Telmer and Yaron (2004); Gervais (2002); Sommer
and Sullivan (2018); Sinai and Souleles (2005); Mankiw and Weil (1989);
Ermisch (1999); Becker and Lewis (1973); Barro and Becker (1989); Jones and
Tertilt (2008); Kleven, Landais and Søgaard (2019); Kleven et al. (2019);
Finer and Zolna (2016); de la Croix and Pommeret (2021). Each is used for its
headline result only.
