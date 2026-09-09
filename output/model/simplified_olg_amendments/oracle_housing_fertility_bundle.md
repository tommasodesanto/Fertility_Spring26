Please answer the following follow-up to your most recent one-date planner response. The benchmark is fixed; the new tasks concern housing restrictions and fertility valued through parents only.

### File: docs/prompts/oracle_simplified_olg_housing_fertility.md
```md
# Follow-up: housing allocation, utility restrictions, and parental fertility

Continue the existing household model and your most recent ONE-DATE planner
response in this chat. The author has now read and discussed it. We need a
compact, economically transparent argument for the paper, with the three tasks
below in order. Please reason carefully and check the proofs before answering;
do not merely produce another provisional menu. The deadline is close, so
prioritize a useful theorem and its clear limits over a sprawling catalogue.

The author is worried that the planner's powers keep changing in response to
inconvenient results. Freeze the benchmark below. Do not alter it to obtain a
desired sign. The only new extension authorized here is to allow fertility
to respond, and then to let the same static planner choose fertility while
valuing ONLY the currently living parents' utility, not unborn persons' utility.

## Maintained model and precise boundaries

Use the exact household budgets already attached in this chat. Preferences are
\[
u^y(c,h,n)=\log(c-\chi n)+\alpha\log(h-\kappa n)+\vartheta\log n,
\qquad
u^o(c^o,h^o,e)=\log c^o+\gamma\log h^o+\omega_B\log e.
\]
All taste weights and child costs are positive. Beta weights a young parent's
own old-age utility. It is not a weight on children's utility. Incoming wealth
is exogenous; warm-glow estates do not define dynastic altruism. Preserve
income/wealth heterogeneity and the logistic ownership taste. Tenure is chosen
young and retained old. Rental and owner physical caps remain, and may bind.
Current young income is available at purchase; the mortgage finances at most
the share phi of the house price. Do not reinstate the old cash-only timing.

For clarity about the issue we just discussed: competitive OLD households
cannot borrow at all. They CAN sell their existing home, downsize or buy a
larger home with available resources. The mortgage from youth is repaid on
entering old age. Old saving satisfies a^e >= 0. For OLD OWNERS,
\[
c^o+a^e+(1+q\tau^p)P_t h^o=a+P_tH+y^o+T_t,
\qquad e=q^{-1}a^e+P_{t+1}h^o.
\]
Old renters instead have e=q^{-1}a^e and receive no housing-sale proceeds.
Thus selling is already possible; no result may rely on owners being unable
to sell, on an inherited-home size limit, or on a new reverse-mortgage market
in the COMPETITIVE model. Your counterexample instead exploited the direct
PLANNER's permission to relax nonnegative old financial saving. We checked
that counterexample, including mortgage accounting, both-tenure averaging,
stationarity, and finite caps chosen above both market and planner demands.

For tasks 1 and 2, retain the last submitted dated allocation benchmark:
same currently living households, initial fertility, retained tenure, total
current consumption C^eq and total housing Hbar; equal weights on living
households' remaining utility; future real opportunities of these households
and old net estate payments fixed. The planner chooses all current c and h.
It can redistribute financial claims and relax BOTH young financing and old
nonnegative saving, while preserving obligations and aggregate feasibility.
The intervening discussion proposed retaining the old borrowing rule but did
not adopt that change. Do not silently impose e/P as an additional old cap.
Actual tenure caps remain physical constraints.

The reference is one date of a positive stationary competitive equilibrium.
That supplies equal cohort masses N and the common normalized type/tenure law
Q. The welfare comparison is one-date, not stationary-cohort welfare or a
permanent transfer scheme. Owner/renter shares coincide across ages because
of this equilibrium's tenure persistence; this is not an extra assumption
and does not imply equal home sizes. The main comparison is INDIVIDUAL first,
aggregate separately. More young housing in total does not mean each young
household receives more, and a local improving transfer does not establish
the full optimum's aggregate direction.

## 1. Utility and primitive restrictions for housing toward the young

The author expects old households to occupy weakly larger homes than young
households. This is an economic expectation to investigate, NOT a fact already
proved for every competitive equilibrium. We have now clarified this chain:
\[
\alpha\ge\gamma,\quad h_i^{o,eq}\ge h_i^{y,eq},\quad n_i>0
\quad\Longrightarrow\quad
\frac{\alpha}{h_i^{y,eq}-\kappa n_i}
>\frac{\alpha}{h_i^{y,eq}}
\ge\frac{\gamma}{h_i^{o,eq}}.
\]
This is a paired marginal comparison, with a feasible housing transfer only
when the young household can receive it below its cap. Separately, the full
fixed-fertility optimum has
\[
h_i^{y,SP}=\min\{h_{d_i}^{\max},\kappa n_i+\alpha/\lambda\},\qquad
h_i^{o,SP}=\min\{h_{d_i}^{\max},\gamma/\lambda\}.
\]
With alpha >= gamma, positive fertility, the common Q, and
Hbar < 2N integral h_d^max dQ, it gives H_Y^SP > Hbar/2. Hence the MARKET
aggregate condition H_O^eq >= H_Y^eq delivers H_Y^SP > H_Y^eq. Verify and use
this simple lemma. Do not replace it with your much longer price-certificate
theorem unless the longer theorem establishes something economically useful.

Explain precisely what alpha >= gamma means, how far it can be weakened by
children's space needs, and which conclusions require a paired housing
ordering versus merely an aggregate ordering. Establish useful conditions
on preferences, initial wealth, old income, mortgage terms and estate tastes
under which equilibrium generates the relevant ordering. Exploit saving and
household optimality rather than only loose all-type cash bounds. Preserve
both old estate regimes and meaningful housing caps. Income ratios may help,
but do not silently require extraordinarily high old income or remove all
heterogeneity. A transparent restricted benchmark can teach the mechanism;
identify explicitly what survives beyond that benchmark.

State the actual inequalities first, then explain their economic content and
whether they seem demanding as THEORY restrictions. We expect high patience
in the quantitative model, but beta there and the two-age beta here do not
map one-to-one. Do not claim calibration compatibility without a mapping or
invent calibration values. Distinguish equal-weight redistribution from an
effect causally attributable to financing constraints.

An important correction to your last answer: with your notation
E=1+alpha+vartheta, K=1+gamma+omega_B, D=E+beta K and
bar w0=integral(y^y+b)dF, your lower-price certificate (18) necessarily implies
\[
\beta<\frac{\nu\vartheta\bar w_0/\chi-E}{K}.
\]
Indeed the minimum inside (18) is bounded above by
vartheta w0/[D(chi+p_lower kappa)]. Thus saying beta is unrestricted was
misleading: at fixed other primitives your sufficient conditions impose a
finite ceiling. Independent checks also give positive stationary equilibria
above that ceiling. This is a limitation of the sufficient bound, not a
necessary condition for the desired housing effect. Address it explicitly
and try to avoid reproducing that artifact in the new main statement.

## 2. When does the housing reallocation raise privately chosen fertility?

Start from task 1's allocation with fertility initially fixed. Now distinguish
two experiments carefully:

First, give a parent a specified feasible current consumption/housing bundle
(c,h), preserve that parent's continuation opportunities, and let the parent
choose n within the bundle. This is a conditional response to an assigned
bundle, not automatically an equilibrium policy with all choices free.
Second, if claiming a grant, tax, credit or other market intervention has the
same effect, specify its financing and which choices and prices respond.
Do not infer a policy sign from an unexplained free housing gift.

The existing note already contains
\[
\frac{\vartheta}{n}
=\frac{\chi}{c-\chi n}+\frac{\alpha\kappa}{h-\kappa n}.
\]
Use and verify this equation to give both a local and a finite-change
fertility test. In particular, for fixed vartheta, x=c-chi n and s=h-kappa n,
the candidate local identity is
\[
dn=
\frac{(\chi/x^2)\,dc+(\alpha\kappa/s^2)\,dh}
{\vartheta/n^2+\chi^2/x^2+\alpha\kappa^2/s^2}.
\]
An increase in housing may come with lower consumption, so higher welfare or
more housing alone is not the fertility theorem. Characterize the consumption
loss the space gain can offset. Try to translate this into meaningful utility
or resource-cost restrictions rather than leaving a borrowing multiplier in
the condition. If a price or allocation object is unavoidable, say so and
separate the exact condition from primitive sufficient restrictions.

Apply the response test to the ACTUAL full fixed-fertility planner bundles,
where possible: which young households gain both consumption and housing,
which face a tradeoff, and when does individual or average fertility rise?
Distinguish evaluating the new bundles at original n from the resulting
optimal n. Applying a conditional fertility choice after task 1 is not the
same as resolving a planner that chooses c, h and n jointly (task 3).
Retained tenure is the maintained allocation comparison; for a genuine policy
with tenure choice free, flag and account for selection before claiming an
aggregate fertility sign.

## 3. The same static planner also chooses fertility, valuing parents only

Now explicitly change ONLY the fertility choice: the planner chooses
{c_i,h_i,n_i} for current young and {c_j^o,h_j^o} for current old. Its objective is
\[
\int_Y[\log(c_i-\chi n_i)+\alpha\log(h_i-\kappa n_i)
             +\vartheta\log n_i]di
+\int_O[\log c_j^o+\gamma\log h_j^o]dj,
\]
with fixed estate and incumbent continuation terms omitted as constants.
Use the SAME current total goods/housing, tenure caps, welfare weights and
financial permissions as above. Each child's goods and space requirements
are already embedded in c-chi n and h-kappa n: do not double-count them in
aggregate resources. Parents enjoy children through vartheta log n. There is
NO separate welfare weight for unborn people and NO dynastic continuation
utility of descendants.

Formulate this problem fully, check concavity and feasibility, characterize
the optimum including binding tenure caps, and compare its n and h with the
competitive reference and with the sequential exercise in task 2. Seek one
simple proposition or useful sufficient inequalities for increased fertility,
individual first and aggregate separately. Do NOT assume that freeing n at
the planner optimum must increase n relative to its fixed reference value.
Define any needed goods/housing resource multipliers plainly before using
them; derive their equations and distinguish solved primitive restrictions
from conditions still involving endogenous allocations or multipliers.

The cohort mass TODAY is fixed. Freely changing fertility means tomorrow's
new cohort need not retain its reference size. It is inconsistent to impose
replacement fertility 1/nu again just to keep the new allocation stationary.
Preserve existing adult commitments; do not promise that unborn cohorts'
entire future allocation stays unchanged when their numbers change. State
the exact static scope without opening a new population-welfare problem.
Changing n must not add implicit value to the unborn through the objective.

## Short bridge to the transition question, after the three core tasks

The quantitative project studies today's economy along a transition following
an exogenous fertility-preference decline, then compares a policy along that
path. This run should not solve or redesign that quantitative model. In a
short final paragraph, identify which fertility conditions above apply at
each date at the same exogenous preference path, and what extra market-
clearing and funding work is needed to infer a POLICY path. Distinguish
finite-date implications from statements requiring convergence. The identity
Y_{t+1}=nu bar n_t Y_t, O_{t+1}=Y_t is already available: fertility ordered
date by date gives cohort ordering from a common initial condition, but it
does not prove a policy raises fertility. Positive stationary endpoints both
have bar n=1/nu; differing population levels can depend on fertility along
the path. Do not replace this with a new stationary-cohort welfare criterion.

## Deliverable and discipline

Deliver one integrated, readable note: (1) housing comparison and interpretable
restrictions, (2) parents' fertility response to the specified reallocation,
(3) the static parent-utility planner with fertility free. Give the main
statements and their economic meaning first, followed by compact proofs.
Aim for roughly 2,000–3,000 words excluding equations, with only a short
technical appendix if genuinely needed. We want a problem-set-style argument,
not an inflated literature review or six alternative frameworks. Optional
functional-form changes belong in a brief clearly labeled discussion only;
none is adopted in this run. If a desired general statement fails, provide
an analytic reason and the strongest useful conditional statement. Do not
substitute a numerical point plus continuity for parameter restrictions, and
do not describe sufficient conditions as mild, calibrated or beta-unrestricted
without checking what their combined inequalities actually imply.

Do not spend the whole response repeating the old counterexample. Use it as
a verified limitation and make progress on the three tasks. Keep our notation,
use lower-case y and o for ages, explain each new object before using it, and
state precisely what is proved versus what remains open. Do not silently
change planner powers, mortgage timing, utility, taxes, tenure or inheritance.
```
