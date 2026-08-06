# Fresh-take prompt: population closure for a housing–fertility lifecycle model

Use this document as the complete context. Take nothing in it as settled: the
menu below reflects one prior discussion and may be badly framed. Your job is
an independent economic assessment. Challenge the framing itself if warranted.

## The model, in brief

One-market quantitative lifecycle model of housing and fertility (US,
national interpretation intended but debated). Households live 17 four-year
periods from age 18; persistent earnings risk; sequential births with an
attempt/taste-shock margin and declining fecundity; children ever born (n)
and children at home (m) tracked separately; renters choose continuous
housing up to a cap, owners a discrete ladder with a 20 percent down payment;
tenure deterministic. Housing supply is an upward-sloping curve in LEVELS,
\(H^S = H_0 (r/\bar r)^\xi\); rents tied to prices by a user-cost identity
\(r=(q+\delta+\tau_H)p\); wages and the interest rate exogenous. Calibrated
completed fertility is ~1.92–1.97 — below replacement. Ten estimated
parameters, twelve moments, stationary equilibrium, SMM.

Policy experiments: property-tax and housing-grant reforms with a balanced
government budget. Household-level fertility responses are small (<0.7
percent) under every design tried; the live question is what population-level
statement the model can honestly make.

## The closure problem

With endogenous below-replacement fertility, a closed stationary population
does not exist: city-born matured children (flow \(B_0 \approx 0.053\) per
unit population per period) cover only ~86 percent of the entrant flow
\(E_0 \approx 0.062\) required for stationarity. The current construction
closes the gap with OUTSIDE ENTRANTS:

\[ S E_0 = q\,[M + S B_0], \qquad q = \textstyle\sum_z \omega_z
\Lambda(\lambda[V_z - \bar V^{out}]) \]

An outside option with representative value \(\bar V^{out}\), taste scale
\(\lambda\) (set to 2, not estimated), and entry probability \(q\). The
empirical anchor: the outside-origin share of entrants is pinned to 0.169
(ACS 2012–2023 cumulative across-CBSA arrivals, ages 18–22), which at the
baseline implies \(q^* \approx 0.97\) — the margin is nearly saturated. The
baseline backs out \(M\) and \(\bar V^{out}\) so population scale \(S=1\):
stationarity is guaranteed regardless of whether fertility is above or below
replacement. Bequests are warm-glow (estates transfer to no one), so there is
no intergenerational wealth flow to complicate the accounting.

## The author's desiderata

1. No full transition dynamics (heterogeneous agents + discrete choices +
   forward-looking prices; too costly, and without a terminal steady state
   there is nothing to shoot toward).
2. Population scale MUST be able to respond in counterfactuals — a
   fixed-population closure is rejected; the author wants a genuine "macro"
   statement (e.g., what happens to prices when fertility changes).
3. The aggregate housing scarcity margin should stay meaningful — long-run
   supply should not become mechanically proportional to population.
4. The construction should not feel ad hoc. The disembodied outside value has
   "always been suspect" to the author.
5. Two months of runway; the model is otherwise mature.

## Counterfactual protocols discussed so far (treat all as suspect)

(i) **Open-city/local**: hold \(\bar V^{out}\), \(M\), \(\lambda\) fixed; the
policy changes household values; \(q\) responds via the logit; \(S\) solves
the identity. Concern: \(\lambda\) is unidentified and drives ~80 percent of
the resulting population responses (+9–10 percent); near-saturation makes
responses asymmetric; interprets a national policy as a single treated city.

(ii) **National/symmetric**: treat the outside as more of the same economy;
under a nationwide policy the outside value moves in parallel, the gap
\(V-\bar V^{out}\) is unchanged, entry rates stay constant by optimality, and
population moves only through births, with multiplier
\((1-s^{out})/s^{out} \approx 4.9\) pinned by the migration share (grant:
births +0.67 percent → stationary population +3.3 percent). Concern (author):
recomputing the outside value inside the experiment feels like moving the
goalposts; fixing flows seems to ignore that higher prices should deter
entrants.

(iii) **Residual valve**: outside inflow adjusts to hold \(S=1\) in
counterfactuals too. Rejected by the author — population must respond.

(iv) **Strip the option**: delete the outside value and logit entirely;
entrants = matured city-born + an exogenous immigration inflow \(M\)
calibrated to the 0.169 share and held fixed in counterfactuals
(quota interpretation). Observationally equivalent to (ii); the question is
whether removing the choice-theoretic layer is a feature (honesty about what
is identified) or a loss (no entrant margin at all, no local experiments).

(v) **Demographic balanced path**: endogenous population growth rate n,
age pyramid tilted by \((1+n)^{-j}\), one outer fixed point. Delivers "policy
changes the population growth rate." Concern: a balanced path forces housing
supply per capita, assuming away the aggregate scarcity margin (desideratum
3); comparisons along declining paths are awkward; the author has a standing
aversion to BGPs.

(vi) **Full transitions**: excluded by desideratum 1.

## Questions

1. Assess the menu. Which protocol (or combination) best supports (a)
   household-level fertility-policy results and (b) a disciplined
   population-level statement? Where is the menu itself wrongly framed?
2. Propose closures NOT on the menu. Examples of directions that were not
   explored: stochastic-aging/perpetual-youth recasts; endogenous marriage or
   household-formation margins as the entry object; land/urban-boundary
   models where scarcity survives population scaling; semi-open capital-
   market analogies; anything from the quantitative OLG, fertility-macro, or
   urban literatures with a track record.
3. What do the relevant literatures actually do? Quantitative OLG with
   below-replacement demographics (e.g., Japan), immigration in GE lifecycle
   models, endogenous-fertility macro, open- vs closed-city urban models,
   spatial QSM population closures. Cite only what you can verify; label
   working papers as such.
4. If entry should respond to prices (the author's instinct), what would
   identify the entry elasticity? Is there a credible empirical moment, or is
   the honest answer a preregistered sensitivity band?
5. State the minimal implementation of your recommended closure and what new
   assumptions it introduces, explicitly.

Answer as an independent senior macro/urban economist. Economics first,
implementation second. Do not defer to the prior discussion's conclusions —
including the claim that (ii) and (iv) are equivalent, and the claim that a
balanced path necessarily sacrifices the scarcity margin. If you conclude the
current outside-option construction is actually fine, say so and say why.
