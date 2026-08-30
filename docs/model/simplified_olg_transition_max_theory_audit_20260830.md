# Max-reasoning theory audit: simplified two-generation OLG note

Date: 2026-08-30
Source audited: `latex/simplified_olg_transition_note.tex` and the four-page
rendered PDF through the 18:29 final candidate
Scope: mathematical audit only; no source file was edited

## Bottom line

The smallest two-generation OLG is capable of carrying the paper's two intended
results:

1. a local, compensated reallocation result between an old incumbent and a
   financing-constrained young owner; and
2. a conditional fertility response to a relaxation of a housing-access limit.

The final candidate gets the architecture right and is mathematically
defensible. I find no fatal issue in Definition 3, Proposition 1, or the
transition corollary. Equations (3)--(14), the old-distribution law in (15), the
cohort identities in (16), the steady-state restrictions (18)--(20), the
marginal-value formulas (22)--(24), and the lifetime/derivative algebra in
(25)--(28) survive re-derivation.

The revision resolves the four blockers in the first draft: it defines old
renter constraints, fixes the tax and sales conventions, supplies the inherited
cohort state and perfect-foresight/fixed-policy language, replaces the
transition Pareto claim with a valuation-gap corollary, and deletes the
arbitrary path plot. In particular, **equal rebates to the cohort entering at
`t+1` are no longer a welfare problem**, because the revised corollary is only
an accounting identity. At a steady state, `\ell_t=\ell_{t+1}=0` and moving a
unit within the fixed stock does not change the aggregate property-tax base, so
the proposition changes no government rebate and involves no unborn household.

The final proof now ties “finance the transaction” to the model's closing
constraint. It advances the pledgeable amount

\[
L=(1-\phi)P^*\,\mathrm dh
\]

at closing and receives `RL` from the buyer next period. Then the binding
constraint remains binding at `h+\mathrm dh`, the contract has zero present
value, and

\[
(\lambda_i+\mu_i)L-\beta V^O_{a,i}RL
=\mu_iL
=\lambda_i\zeta_i^{O,F}\,\mathrm dh,
\qquad \lambda_i=\beta RV^O_{a,i}.
\]

This is the exact budget-feasible source of the financing surplus. The final
pair setup also states that the old incumbent's estate-composition and
retention constraints are slack. Consequently Proposition 1 is a valid local
constrained-Pareto result, not merely a marginal-value analogy.

The final prose also correctly separates the welfare and fertility results.
That separation matters because equation (26) is a derivative at fixed
lifetime wealth. For reference, if the young household instead pays
`x\,\mathrm dh` after closing, `r_t^m+x` replaces `r_t^m` in the numerator's
resource-cost term and the exact threshold becomes

\[
\chi\leq
\frac{(1+k)\kappa(r_t^m+\zeta_t^m)^2}
{\alpha(r_t^m+x)}.
\]

The simple condition (28) still suffices whenever
`0\leq x\leq\zeta_t^m`, but the note does not need this extra funded-incidence
result.

The safest five-page architecture is:

- a **transition valuation-gap lemma** (an accounting identity, not a welfare
  theorem);
- a **steady-state compensated-reallocation proposition** (the clean headline
  welfare result);
- a short transition remark stating the additional two-date transfer technology
  under which the valuation gap becomes a local Pareto surplus; and
- the fertility derivative with its exact sign threshold, followed by the
  simpler primitive sufficient condition already displayed.

If the author later wants a full transition Pareto proposition, it can be made
valid without welfare weights over unborn agents using the exact dated
construction in Section 5. It requires explicit bridge credit and two-date side
transfers. An untaxed "occupancy transfer" is not the right repair under the
current primitives: it would make the old household a landlord or introduce a
hybrid housing trust, neither of which exists in this model.

**Recommended figure count: exactly one.** Keep the revised marginal-value
wedge as a schematic implication of the steady-state proposition. The deleted
transition-path plot should stay deleted; a timing diagram is optional prose
support, not needed for either result in a four-page note.

## 1. Primitive and cash-flow audit

### 1.1 Domains and institutional statements that must be explicit

The log problems require

\[
c>0,\qquad n>0,\qquad h-\kappa n>0,
\qquad c^O>0,\qquad h^O>0,\qquad e>0.
\]

These domains are not cosmetic. They rule out zero fertility and make this a
continuous completed-fertility toy model, not a model of childlessness or birth
hazards. The final candidate now states them and correctly avoids an
extensive-margin fertility claim.

The final financial-sector sentence now says that a competitive sector trades
a claim delivering one unit of next-period goods at price `q`, finances
mortgages, and buys stepped-up estate housing at the market price. This is the
correct distinction: `q` prices future goods, while retained housing enters the
estate at `P_{t+1}`.

The following ownership rules are also required for the algebra already in the
note:

- young renters remain renters when old;
- young owners become old incumbents;
- an old incumbent can sell part of a divisible house and self-occupy the
  retained part, but cannot become a renter or a household landlord;
- retained housing is stepped up and liquidated by the estate/financial sector
  after death;
- the capital-gains tax applies to household sales only; competitive
  intermediaries mark housing to market without an additional household
  capital-gains liability; and
- a household that buys at `t-1` has basis `P_{t-1}` when it is old at `t`.

The last two bullets are needed for `\ell_t=\tau^g[P_t-P_{t-1}]_+` to be the
complete tax state in a two-period life.

### 1.2 Competitive user cost

An intermediary purchasing one unit at `t`, paying the current property tax,
renting it for `r_t`, and selling it at `t+1` has date-`t` profit

\[
r_t+qP_{t+1}-(1+\tau_t^p)P_t.
\]

Zero profit gives

\[
r_t=(1+\tau_t^p)P_t-qP_{t+1}.
\]

Equation (3) is therefore correct. It assumes that `r_t` is the rent paid for
one unit of services and that the property tax applies to the full occupied
stock.

### 1.3 Buyer user cost

A young buyer pays the down payment and property tax at `t`, repays the mortgage
at `t+1`, and, at the margin used in the note, sells the extra unit when old.
The date-`t` cost is

\[
\begin{aligned}
&[(1-\phi)+\tau_t^p]P_t
+q\{R\phi P_t-(P_{t+1}-\ell_{t+1})\}\\
&\quad=(1+\tau_t^p)P_t-q(P_{t+1}-\ell_{t+1})
=r_t+q\ell_{t+1}.
\end{aligned}
\]

Thus equations (5) and (12) correctly incorporate the down payment, mortgage
principal, property tax, resale value, and next-period capital-gains tax. The
result relies on a zero-spread mortgage, `R=q^{-1}`.

### 1.4 Incumbent user cost and step-up

An old incumbent entering with `H` can sell a unit now for `P_t-\ell_t`, or
retain it, pay `\tau_t^pP_t`, and put `P_{t+1}` in the stepped-up estate. The
opportunity cost of retention is

\[
(P_t-\ell_t)+\tau_t^pP_t-qP_{t+1}=r_t-\ell_t.
\]

Equation (6) is correct. The sign has the intended lock-in interpretation:
realizing the current gain makes sale privately less attractive, so retention
has a lower private user cost. For the interior FOC used later, the active set
must satisfy `r_t-\ell_t>0`; otherwise the upper retention bound is generally
active.

### 1.5 Old household problem

Let liquid estate wealth be `x` and total estate goods be
`e=x+P_{t+1}h^O`. Starting from

\[
c^O+qx+\tau_t^pP_th^O
=a+(P_t-\ell_t)(H-h^O)+T_t
\]

and substituting `x=e-P_{t+1}h^O` gives exactly

\[
c^O+qe+(r_t-\ell_t)h^O
=a+(P_t-\ell_t)H+T_t,
\qquad e\ge P_{t+1}h^O.
\]

Equation (7) therefore does not double-count retained housing. The inequality
is the nonnegative-liquid-estate restriction.

The old-renter problem must be written as

\[
\max_{c^O,h^O,e>0}u^O(c^O,h^O,e)
\quad\text{s.t.}\quad
c^O+qe+r_th^O=a+T_t,
\qquad 0<h^O\le \bar h^R.
\]

An old renter has no housing component in the estate, so only `e>0` is needed.
The final candidate now states these constraints explicitly rather than
importing the incumbent's `h^O\le H` and estate-composition constraints.

### 1.6 Young owner problem and down-payment timing

Equation (10) is internally consistent. The current budget pays only the down
payment plus property tax,

\[
[(1-\phi)+\tau_t^p]P_th,
\]

while next-period net financial wealth subtracts mortgage repayment,

\[
a_{t+1}=Ra'-R\phi P_th.
\]

Discounting the repayment gives `qR\phi P_th=\phi P_th`, so the full purchase
price is recovered. The separate collateral condition

\[
(1-\phi)P_th\le b_i
\]

is also consistent with the stated closing timing: current income and the
ordinary rebate arrive too late to be pledgeable. Any welfare intervention
that relaxes this condition must therefore say that its transfer is available
*before closing*; an ordinary after-closing lump-sum rebate is not enough.

## 2. Envelopes, FOCs, and access wedges

### 2.1 Old incumbent

Put multipliers `\lambda^O`, `\rho`, and `\sigma` on the old budget, liquid
estate constraint, and retention constraint. A convenient Lagrangian is

\[
\begin{aligned}
\mathcal L^O={}&u^O(c^O,h^O,e)
+\lambda^O\{a+(P_t-\ell_t)H+T_t-c^O-qe-(r_t-\ell_t)h^O\}\\
&+\rho(e-P_{t+1}h^O)+\sigma(H-h^O).
\end{aligned}
\]

The FOCs and envelopes are

\[
\lambda^O=\frac1{c^O},\qquad
\frac{\omega_B}{e}=q\lambda^O-\rho,
\qquad
\frac{\gamma}{h^O}=\lambda^O(r_t-\ell_t)+\rho P_{t+1}+\sigma,
\]

\[
V_a^O=\lambda^O,
\qquad
V_H^O=\lambda^O(P_t-\ell_t)+\sigma.
\]

When both `e>P_{t+1}h^O` and `h^O<H` are slack, `\rho=\sigma=0` and

\[
\frac{u_h^O}{u_c^O}=\frac{\gamma c^O}{h^O}=r_t-\ell_t,
\qquad
\frac{V_H^O}{V_a^O}=P_t-\ell_t.
\]

Equation (8) is correct. The final intergenerational-allocation setup now
carries *both* slackness assumptions forward. Slack retention alone would not
be enough if the estate-composition constraint bound.

### 2.2 Young renter

With `\lambda_t^R` on the budget and `\eta_t^R` on the rental cap, interior
saving gives

\[
u_c^Y=\lambda_t^R,
\qquad
\lambda_t^R=\beta R V_{a,t+1}^R,
\qquad
u_h^Y=\lambda_t^R r_t+\eta_t^R.
\]

Since `u_c^Y=1/c`, division by `\lambda_t^R` yields

\[
\frac{\alpha c_t^R}{h_t^R-\kappa n_t^R}
=r_t+\frac{\eta_t^R}{\lambda_t^R}.
\]

Thus `\zeta_t^R=\eta_t^R/\lambda_t^R` has units of current goods per unit
of housing and is exactly the renter-cap shadow value.

### 2.3 Young owner

Let

\[
D_t=[(1-\phi)+\tau_t^p]P_t.
\]

With the signs used in the note, the owner housing FOC is

\[
u_h^Y-\lambda_t^OD_t
+\beta\{V_{H,t+1}^O-R\phi P_tV_{a,t+1}^O\}
-\mu_t(1-\phi)P_t-\eta_t^O=0.
\]

Interior saving gives

\[
\lambda_t^O=\beta R V_{a,t+1}^O,
\qquad
\frac{\beta V_{a,t+1}^O}{\lambda_t^O}=q.
\]

If the next-period old household's estate-composition and retention constraints
are both slack, then

\[
V_{H,t+1}^O=(P_{t+1}-\ell_{t+1})V_{a,t+1}^O.
\]

Substitution gives

\[
\begin{aligned}
\frac{u_h^Y}{u_c^Y}
&=D_t+\phi P_t-q(P_{t+1}-\ell_{t+1})
+\frac{\mu_t}{\lambda_t^O}(1-\phi)P_t
+\frac{\eta_t^O}{\lambda_t^O}\\
&=r_t+q\ell_{t+1}+\zeta_t^{O,F}+\zeta_t^{O,S}.
\end{aligned}
\]

This verifies equations (11)--(12), including the financed-share convention:
the down payment is `(1-\phi)`, not `\phi`.

### 2.4 Fertility FOC

For either tenure,

\[
\frac{\vartheta_t}{n}
-\frac{\alpha\kappa}{h-\kappa n}
-\lambda_t^m\chi=0.
\]

Using `\lambda_t^m=1/c` and

\[
\frac{\alpha c}{h-\kappa n}=r_t^m+\zeta_t^m
\]

gives

\[
\frac{\vartheta_tc}{n}
=\chi+\kappa(r_t^m+\zeta_t^m).
\]

Equation (14) is correct. It is a conditional interior FOC; it is not a result
about the extensive margin of childlessness.

## 3. Cohorts, government, and positive steady states

### 3.1 Cohort law

With one model period equal to a generation and no survival beyond old age,

\[
Y_{t+1}=\nu\bar n_tY_t,
\qquad
O_{t+1}=Y_t
\]

is correct. The final note defines the relevant averages over `F` and `G_t`;
explicitly, they are

\[
\bar n_t=\int n_t(i)\,\mathrm dF(i),
\qquad
\bar h_t^Y=\int h_t(i)\,\mathrm dF(i),
\]

and should reserve `\bar h_t^O` for an average over the inherited old-state
distribution. The fixed type distribution of each newly entering young cohort
is a primitive; the children do not inherit their parent's individual state.

### 3.2 Old-state distribution and sales

Let `G_t` be a probability distribution across old households over net
financial wealth, inherited housing, and tenure. For any test function `f`, the
law generated by the date-`t` young policies is

\[
\begin{aligned}
\int f\,\mathrm dG_{t+1}
=\int \bigg[&\mathbf 1\{m_t(i)=R\}
 f\big(Ra_t^{\prime R}(i),0,R\big)\\
&+\mathbf 1\{m_t(i)=O\}
 f\big(Ra_t^{\prime O}(i)-R\phi P_th_t^O(i),h_t^O(i),O\big)
\bigg] \mathrm dF(i).
\end{aligned}
\]

Equation (15) in the final candidate is this distribution law in compact
pushforward notation. Because `O_{t+1}=Y_t`, `G_{t+1}` is correctly normalized
to integrate to one; cohort mass is carried separately.

Define sales per *old household*, including zeros for renters, by

\[
\bar s_t^O
=\int \mathbf 1\{H>0\}\,[H-h_t^O(a,H)]\,\mathrm dG_t(a,H,m).
\]

Then total taxable old-household sales equal `O_t\bar s_t^O`. If the author
instead wants `\bar s_t^O` to mean sales conditional on being an incumbent,
the budget must use `O_t\pi_t^O\bar s_{t\mid I}^O`, where `\pi_t^O` is the old
incumbent share. The final candidate uses the first convention, assigning zero
sales to old renters.

### 3.3 Government budget

Under full occupancy, the property-tax base is the fixed stock `\bar H`. With
the sales definition immediately above, the exact balanced-rebate identity is

\[
(Y_t+O_t)T_t
=\tau_t^pP_t\bar H+\ell_tO_t\bar s_t^O.
\]

Equation (17) is correct with the final averaging convention. The final
institutional text also states that stepped-up estate liquidations and
intermediary mark-to-market trades do not generate `\ell_t`.

### 3.4 Positive steady state

A positive steady state must include a stationary old distribution `G^*` and
constant tax rates as well as `P^*`, `T^*`, and `Y^*=O^*>0`. Its restrictions
are

\[
1=\nu\mathcal N(P^*,T^*;\vartheta),
\]

\[
\bar H=Y^*\{\mathcal H^Y(P^*,T^*;\vartheta)
+\mathcal H^O(P^*,T^*;\vartheta)\},
\]

\[
2Y^*T^*=\tau^pP^*\bar H,
\qquad
G^*=\Gamma(P^*,T^*;\vartheta).
\]

The first three are equations (18)--(20); the final candidate now states the
distribution fixed point explicitly. Every positive steady state necessarily has

\[
\mathcal N(P^*,T^*;\vartheta)=\frac1\nu.
\]

This is worth saying in words: a change in `\vartheta` does not permit a
different replacement fertility level in a positive stationary population.
Instead, prices, transfers, scale, and policies must adjust so fertility again
equals `1/\nu`.

The three aggregate equations jointly *restrict* price, transfer, and scale;
they do not prove that a solution exists or is unique. The correct statement is
"For values of `\vartheta` for which a positive solution exists, let
`\mathcal E(\vartheta)` denote one such equilibrium." Because `P_t=P_{t-1}`
at a stationary point, `\ell^*=0`; that conclusion is correct.

## 4. Perfect-foresight transition and state sufficiency

### 4.1 Minimal state

The predetermined state at date `t` is

\[
S_t=(Y_t,O_t,G_t,P_{t-1}).
\]

`P_{t-1}` is needed because it is the common tax basis of the old owners who
bought when young. No longer basis history is needed in this two-period life.
The newly entering young draw fresh types from `F`, so no parental state is
carried into their distribution.

Given the anticipated sequence `\{P_s,T_s\}_{s\ge t}`, the fixed tax rules,
and `\vartheta_t`, household policies determine `\bar n_t`, `\bar h_t^Y`,
old sales, and `G_{t+1}`. Housing clearing and the government budget determine
the contemporaneous equilibrium restrictions; the cohort laws determine
`Y_{t+1}` and `O_{t+1}`.

### 4.2 Date-zero inheritance

If the economy is in `\mathcal E^0` through date `-1` and the preference
changes at date zero, the complete initial condition is

\[
Y_0=O_0=Y^0,
\qquad
G_0=G^0,
\qquad
P_{-1}=P^0.
\]

The final definition now states the cohort masses and `P_{-1}` and says that the
old enter with the state inherited from `\mathcal E^0`. Naming that inherited
distribution `G_0=G^0` would make the same condition fully checkable; no new
economic state is required.

### 4.3 Correct equilibrium definition

A state-complete definition lists

\[
\{P_t,T_t,Y_t,O_t,G_t,
m_t(\cdot),c_t(\cdot),h_t(\cdot),n_t(\cdot),a_t'(\cdot)}_{t\ge0}
\]

such that:

1. households take the entire price and transfer sequence as known and solve
   their problems at every date;
2. `G_{t+1}` is generated by the explicit policy law in Section 3.2;
3. housing clears and the balanced-rebate identity holds at every date;
4. the intermediary zero-profit condition determines `r_t` from `P_t` and
   `P_{t+1}`;
5. the cohort laws advance `Y_t` and `O_t`;
6. the date-zero state is the one displayed above; and
7. the sequence converges to a selected positive equilibrium `\mathcal E^1`.

The final definition now states that the property- and capital-gains-tax
schedules remain fixed and labels the path perfect foresight. These are the
right closure statements: equation (3) alone would not say whether households
know the path.

This is a definition conditional on existence. The note proves neither
existence, uniqueness, convergence, nor monotonicity. That is acceptable for a
five-page mechanism note, but the scoping sentence must be explicit.

### 4.4 What is predetermined and what may jump

At `t=0`, `Y_0`, `O_0`, and `G_0` are predetermined. `P_0`, `T_0`, and current
policies are forward-looking equilibrium objects and may jump, but the model
does not prove that `P_0\ne P^0` or determine the direction of the jump. The
first endogenous cohort response is

\[
Y_1=\nu\bar n_0Y_0.
\]

Those are the only transition-timing statements established without solving
an example.

## 5. Welfare claim: exact repair

### 5.1 Terminology

The revised steady-state result uses the weakest defensible object: a **local
constrained Pareto test** between one current old household and one current
young household, whose lifetime utility includes her own old age. It is Pareto,
not utilitarian: both must be weakly better off and at least one strictly better
off. No unborn household enters the comparison.

For the optional transition extension, the analogous object is a local
two-date constrained Pareto test. Every household affected at dates `t` and
`t+1` must be weakly better off and at least one must be strictly better off;
again no cardinal welfare weights are needed.

The word "constrained" records what is held fixed: the price sequence, cohort
masses, fertility, tenure assignments, and ordinary saving choices. It also
records the intervention set. To turn the transition valuation gap into a
Pareto result, the institution must be allowed to:

1. execute a taxable marginal sale from an old incumbent to a young owner;
2. issue pledgeable pre-closing bridge credit at `t`, repaid at `t+1` at `R`;
3. trade one-period claims at price `q`;
4. make balanced lump-sum side transfers at `t` and `t+1`; and
5. pay a real reallocation cost `K_{ijt}` measured in date-`t` goods.

Items 2 and 4 are required only for a transition Pareto theorem. The final
Proposition 1 now includes item 2 explicitly. Because capital-gains taxes
vanish at the steady state and the property-tax base is unchanged, it needs no
dated rebate construction. An untaxed occupancy transfer is not a harmless
shortcut: under the current model it
would leave the old household owning housing occupied by the young, creating a
household landlord or a hybrid tenure institution. Use the taxable sale, or add
the new institution openly.

### 5.2 Transition valuation-gap identity

For an old incumbent whose estate-composition and retention constraints are
slack,

\[
MV_{jt}^O=r_t-\ell_t.
\]

For a young owner with interior saving, a binding down-payment constraint, a
slack owner-size limit, and slack next-period old constraints,

\[
MV_{it}^Y=r_t+q\ell_{t+1}+\zeta_{it}^{O,F}.
\]

Therefore

\[
MV_{it}^Y-MV_{jt}^O
=\zeta_{it}^{O,F}+\ell_t+q\ell_{t+1}.
\]

This identity is fully established by the current household FOCs. It is safe
to call it the **transition valuation gap** even if the paper does not adopt
the expanded compensation institution.

### 5.3 Exact two-date compensated construction

Let `\varepsilon>0` be a small unit of the old incumbent's retained housing.
Keep active sets fixed.

**Step 1: finance the young purchase.** The authority lends the young household

\[
L_t=(1-\phi)P_t\varepsilon
\]

before closing and requires repayment `RL_t` at `t+1`. The loan has zero
present value:

\[
L_t-qRL_t=0.
\]

Because the down-payment constraint remains binding and the size cap remains
slack, the pledgeable loan raises the young household's housing by exactly
`\varepsilon`. The current resource value of the loan cancels its future
repayment by the saving Euler equation:

\[
\lambda_{it}L_t-
\beta V_{a,i,t+1}^ORL_t=0,
\qquad
\lambda_{it}=\beta RV_{a,i,t+1}^O.
\]

What remains is the envelope value of relaxing the down-payment constraint:

\[
\mu_{it}L_t
=\lambda_{it}\zeta_{it}^{O,F}\varepsilon.
\]

Thus the young household gains
`\zeta_{it}^{O,F}\varepsilon` in date-`t` goods-equivalent units before any
additional compensation levy.

**Step 2: compensate the old seller.** The old incumbent sells an additional
`\varepsilon` and reduces retention by that amount. At an interior retention
choice,

\[
\frac{1}{u_{c,jt}^O}\frac{\mathrm dV_{jt}^O}{\mathrm d\varepsilon}
=-\left[\frac{u_{h,jt}^O}{u_{c,jt}^O}-(r_t-\ell_t)\right]=0.
\]

The extra net sale proceeds and the change in estate resources exactly
compensate the marginal utility loss to first order. No welfare weight is used.

**Step 3: write the two government accounts.** The fixed property-tax base is
unchanged. The extra old sale realizes current tax `\ell_t\varepsilon`. Since
the young buyer's next-period retention constraint is slack, the marginal
extra unit is sold at `t+1` rather than retained to death; that sale realizes
`\ell_{t+1}\varepsilon`. Let `N_s=Y_s+O_s`. The mandated equal rebates change
by

\[
N_t\,\mathrm dT_t=\ell_t\varepsilon,
\qquad
N_{t+1}\,\mathrm dT_{t+1}=\ell_{t+1}\varepsilon.
\]

This is the step that would be needed to upgrade the revised transition
valuation-gap corollary into a Pareto theorem.

**Step 4: hold every other household harmless.** At each date, the compensation
authority may levy at most the incremental rebate from any household. If it
levies exactly `\mathrm dT_s`, that household's rebate plus side transfer is
zero relative to the baseline. This applies in particular to every young
household entering at `t+1`: it is held exactly at its baseline allocation, so
no welfare weight over an unborn agent is required. Reclaiming all incremental
rebates produces date-`t` present value

\[
\ell_t\varepsilon+q\ell_{t+1}\varepsilon.
\]

The authority can additionally charge the initially constrained young
household up to `\zeta_{it}^{O,F}\varepsilon` in present value without making
it worse off, because that is the bridge-credit gain from Step 1.

**Step 5: pay the real cost and distribute the residual.** After paying
`K_{ijt}\varepsilon`, the available present-value surplus is

\[
\left[
\zeta_{it}^{O,F}+\ell_t+q\ell_{t+1}-K_{ijt}
\right]\varepsilon.
\]

If the bracket is strictly positive, choose the side levies so that no rebate
recipient loses relative to baseline, charge the young no more than its bridge
gain, pay the cost, and leave a positive residual with at least one household.
Every affected household is then weakly better off and one is strictly better
off. By differentiability, a strict first-order inequality supports a genuine
improvement for sufficiently small `\varepsilon`.

This proves the transition formula displayed in the first draft, but only after
the transfer technology and the two government accounts are made explicit. The
revised candidate correctly demotes that formula to the accounting identity in
equation (24).

### 5.4 Recommended theorem hierarchy

For the main five-page note, use the following hierarchy.

**Lemma (transition valuation gap).** Under the stated active sets,

\[
MV_{it}^Y-MV_{jt}^O
=\zeta_{it}^{O,F}+\ell_t+q\ell_{t+1}.
\]

This is pure household accounting and requires no welfare language.

**Proposition (steady-state compensated reallocation).** At a positive steady
state, permit a small taxable sale, zero-present-value bridge credit, and
balanced lump-sum compensation. Because `\ell_t=\ell_{t+1}=0`, the local
Pareto surplus is

\[
\mathcal S_{ij}^*=\zeta_i^{O,F}-K_{ij}.
\]

If `\zeta_i^{O,F}>K_{ij}`, a sufficiently small transfer from the interior old
incumbent to the constrained young owner is a constrained Pareto improvement.

**Remark (transition upgrade).** With the explicit two-date side-transfer
technology in Section 5.3, the transition valuation gap becomes the compensated
surplus

\[
\mathcal S_{ijt}
=\zeta_{it}^{O,F}+\ell_t+q\ell_{t+1}-K_{ijt}.
\]

This ordering makes the clean steady-state financing result the headline and
does not ask the reader to accept a cross-cohort welfare theorem from one line
about rebates.

The result is sufficient, not necessary. `\mathcal S\le0` does not establish
conditional efficiency, because other household pairs, active sets, or
reallocation instruments may still generate gains.

## 6. Lifetime resource identity and fertility response

### 6.1 Old allocation and saving Euler equation

Let

\[
A=1+\gamma+\omega_B,
\qquad
k=\beta A.
\]

When the old household's housing and estate-composition constraints are slack,
Cobb-Douglas expenditure shares give

\[
c_{t+1}^O=\frac{X_{t+1}}{A},
\]

where `X_{t+1}` is total old resources. The young saving FOC gives

\[
\frac1{c_t}=\beta R\frac1{c_{t+1}^O},
\qquad
c_{t+1}^O=\beta Rc_t,
\qquad
qX_{t+1}=kc_t.
\]

For a young owner,

\[
X_{t+1}
=Ra_t'-R\phi P_th_t
+(P_{t+1}-\ell_{t+1})h_t+T_{t+1}.
\]

Multiplying by `q` and substituting into the current budget gives

\[
(1+k)c_t+\chi n_t
+[r_t+q\ell_{t+1}]h_t
=y_i+b_i+T_t+qT_{t+1}.
\]

For a renter the same derivation replaces the bracketed owner cost by `r_t`.
Thus equation (25),

\[
(1+k)c+\chi n+r_t^mh=\mathcal W_t,
\qquad
\mathcal W_t=y_i+b_i+T_t+qT_{t+1},
\]

is correct. It requires interior young saving and slack old housing/estate
constraints. It is a lifetime *private budget identity*, not a social resource
constraint.

### 6.2 Derivative at a fixed access limit

Let

\[
D=\mathcal W_t-r_t^m\widehat h-\chi n,
\qquad
S=\widehat h-\kappa n,
\qquad
c=\frac{D}{1+k}.
\]

The interior fertility FOC can be written

\[
F(n,\widehat h)
=\frac{\vartheta_t}{n}
-\frac{\alpha\kappa}{S}
-\frac{(1+k)\chi}{D}=0.
\]

Its derivatives are

\[
F_n=-\frac{\vartheta_t}{n^2}
-\frac{\alpha\kappa^2}{S^2}
-\frac{(1+k)\chi^2}{D^2}<0,
\]

\[
F_{\widehat h}
=\frac{\alpha\kappa}{S^2}
-\frac{(1+k)\chi r_t^m}{D^2}.
\]

Therefore

\[
\frac{\mathrm dn}{\mathrm d\widehat h}
=\frac{
\dfrac{\alpha\kappa}{S^2}
-\dfrac{(1+k)\chi r_t^m}{D^2}}
{\dfrac{\vartheta_t}{n^2}
+\dfrac{(1+k)\chi^2}{D^2}
+\dfrac{\alpha\kappa^2}{S^2}}.
\]

Equation (26) is exactly correct, including every sign.

### 6.3 Exact sign threshold and the simpler sufficient condition

At the constrained allocation, the housing FOC implies

\[
\frac{D}{S}
=\frac{(1+k)(r_t^m+\zeta_t^m)}{\alpha}.
\]

For `r_t^m>0`, substitution into the numerator gives the exact local sign
condition

\[
\frac{\mathrm dn}{\mathrm d\widehat h}\ge0
\quad\Longleftrightarrow\quad
\chi\le
\frac{(1+k)\kappa(r_t^m+\zeta_t^m)^2}
{\alpha r_t^m}.
\]

Since `\zeta_t^m\ge0`,

\[
\frac{(r_t^m+\zeta_t^m)^2}{r_t^m}\ge r_t^m.
\]

Consequently the displayed condition (28),

\[
\chi\le\frac{(1+k)\kappa r_t^m}{\alpha},
\]

is a valid but stronger sufficient condition. It should be preceded by the
exact threshold. That is more informative and makes clear why a tighter access
wedge can enlarge the region in which more housing raises fertility.

The derivative holds `\mathcal W_t`, tenure, prices, and the active set fixed.
It is not the sign of a funded policy or a general-equilibrium comparative
static. The final paragraph should not say that any institution reaching a
positive financing gap necessarily raises fertility. It does so only if its
intervention locally raises `\widehat h` without an offsetting wealth/price
change and the sign condition holds.

## 7. Figure audit

### Deleted transition-path figure

**Verdict: deletion was correct.** The first-draft figure did not establish a
result. It assumed `P^1<P^0`, `Y^1<Y^0`, an impact price decline, and monotone
convergence. The model establishes only the timing facts that the date-zero
cohorts and old distribution are inherited, while
`Y_1=\nu\bar n_0Y_0`. The revised prose now states those facts without drawing
an invented equilibrium path.

A sign-free timing diagram would be logically harmless, but it would merely
repeat the definition. It is not worth a second figure in a four-page note.

### Current Figure 1: marginal-value wedge

**Verdict: retain.** It is a schematic implication of the steady-state
proposition, not independent evidence or a proof. Its vertical gap is correctly
`\zeta_i^{O,F}` when `\ell=0`, and the caption correctly labels it gross of the
real implementation cost.

The ordering `h_i>h_j^O` is not a mathematical error: the result compares
marginal valuations and does not require the old household to occupy more
housing in levels. Reverse the ordering only if the author wants the picture to
carry that additional empirical story, which the model has not proved.

### Exact count

Use exactly **one** figure: the current steady-state marginal-value wedge.
Do not restore a transition path or add a numerical transition plot until an
explicit transition is solved.

## 8. What can be frozen and what must change

### Safe to freeze after small wording clarifications

- young and old preferences, with positivity domains;
- competitive user cost (3);
- capital-gains tax definition (4), conditional on the stated one-period basis;
- buyer and incumbent user costs (5)--(6);
- old-owner budget accounting (7);
- young renter and owner problems (9)--(10), including `(1-\phi)` as the down
  payment;
- access-wedge definitions and housing FOCs (11)--(12), with the future-old
  retention constraint slack;
- observable access-value identity (13);
- fertility FOC (14);
- old-distribution pushforward (15);
- cohort and housing identities (16);
- government budget (17), with the revised all-old sales convention;
- positive steady-state equations (18)--(20), conditional on existence;
- preference shock (21);
- marginal-value and valuation-gap algebra (22)--(24);
- lifetime resource identity (25);
- fertility derivative (26);
- exact sign threshold (27); and
- primitive sufficient condition (28).

### Remaining nonfatal wording improvements

All mathematical claim blockers are resolved. The 18:29 build also resolves
the last presentation issue: Proposition 1's proof is contiguous and the sole
figure appears after the analytic results at the bottom of page 4.

### Minimality discipline

No additional state, second shock, construction sector, spatial market,
sequential-birth architecture, dynastic inheritance kernel, or quantitative
calibration is needed for these two results. The exact repairs above preserve
the intended two-generation model and fit within a five-page note if the full
two-date compensation construction is kept to one proof paragraph or moved to
a short appendix.

## Verification

The 18:29 source and repository PDF compile without LaTeX warnings, errors,
overfull boxes, or underfull boxes. The result is four pages. All four pages
were rendered and inspected. The audit modified only this report; it did not
modify the repository TeX or PDF.
