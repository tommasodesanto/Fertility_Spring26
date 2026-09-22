# Timing review: multi-year earnings, entry wealth, and house purchase

This is a bounded primary-source review of the timing and aggregation objects
that matter for the four-year model period. It is a red-team memo, not an
adopted earnings specification. The two earlier proposals remain withdrawn or
unadopted: a one-state AR(1) collapse and the endpoint persistent state plus
averaged iid risk. The proposed joint block-income/next-persistent-state
kernel has not been implemented or adopted.

## Finding

The defensible coarse-period precedent is to define the income object at the
model frequency and to observe the current-period state before the household
chooses. Chambers, Garriga, and Schlagenhauf (2009) use three-year model
periods, report a
five-state productivity process whose values and transition matrix are already
for that three-year horizon, and let current income enter the same-period
purchase budget. De Nardi (2004) does the analogous measurement exercise at a
five-year frequency by aggregating complete five-year PSID cells before
estimating the process. Sommer, Sullivan, and Verbrugge (2013) use annual
periods and make the timing explicit: the within-period income shock is
observed before households choose deposits, mortgage debt, and housing.

These papers establish a timing convention and a measurement principle. They
do not establish an annual-to-four-year shortcut. Directly estimating a
multi-year cell and then treating that cell as current income still makes the
cell available at the coarse-period decision; it is a deliberate information
assumption, not an exact representation of annual decisions inside the cell.
In particular, the identity

\[
 z_{t+4}=\rho^4 z_t+\sum_{j=1}^{4}\rho^{4-j}\eta_{t+j}
\]

is an identity for the persistent latent state. It is not an identity for the
four-year earnings total

\[
 Y_t^{(4)}=\sum_{j=0}^{3}
 \exp\{h_{t+j}+z_{t+j}+u_{t+j}\},
\]

or for its mean \(Y_t^{(4)}/4\). The total depends on the full within-block
path, including the persistent path and transitory shocks. If a realized
\(Y_t^{(4)}\) is placed in the state at the start of the block, the household
has been told the future annual shocks. That changes purchase eligibility and
saving behavior relative to the standard coarse-time process.

The project therefore needs one explicit choice before any joint refit:

1. **Direct four-year object.** Construct complete four-year income cells using
   the intended income concept, estimate a process at four-year frequency, and
   treat the observed current cell as current resources. This is closest to De
   Nardi's measurement discipline and Chambers et al.'s three-year
   implementation, but it still deliberately makes the completed cell current
   information at the four-year decision.

2. **Annual process with a four-year approximation.** Simulate annual paths,
   aggregate them to four-year totals or means, and fit a finite process while
   preserving the information set. A state observed at the block start may
   contain the current annual flow and the beginning persistent state; the
  remaining annual path is integrated in the continuation value. It may not
  contain an ex-post four-year total unless the model deliberately assumes
  advance knowledge of that total. This synthetic approximation therefore also
  departs from annual sequential decisions; its advantage is that the departure
  is explicit and can be checked against the simulated annual paths.

3. **Exact within-block aggregation.** Add annual subperiod assets and the
   within-block purchase/earnings timing. This is the only route that makes a
   house purchase during the block and a four-year realized income sum jointly
   exact without a stronger information assumption.

The current endpoint approximation can remain a labeled diagnostic, but its
7.84% excess continuous variance is not by itself evidence that the household
model fails. It is evidence that the endpoint and the four-year flow have been
combined by an approximation whose timing and distribution need to be stated.

## Evidence from the primary sources

### Chambers, Garriga, and Schlagenhauf (2009), *International Economic Review*

Primary sources: [published article record and DOI](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1468-2354.2009.00544.x),
[author-compatible full PDF](https://cooperative-individualism.org/chambers-matthew_accounting-for-changes-in-the-homeownership-rate-2009-aug.pdf).
The article is pp. 677–726; the cited page numbers below are the printed pages.

* The model uses three-year periods. The demographic parameterization starts
  individuals at age 20 and labels age 83 as model period 23 (p. 701); the
  capital-gain and income sections state explicitly that the reported shocks
  are adjusted to the three-year horizon (pp. 701–702).
* The current productivity state is in the household state, and current
  after-tax income enters the same-period budget. In the renter-to-owner
  problem, the household chooses consumption, next-period assets, and a house
  while the right-hand side contains current income (Eq. (6) and Eq. (9),
  pp. 690–693). The continuation expectation is over the next productivity
  state. Thus, the current state is available to the current purchase choice;
  the next state is not.
* The purchase timing is unusually explicit. A footnote says the mortgage
  payment is made in the period in which the mortgage is written and that the
  household can purchase a home and consume its service flow in the same period
  (p. 689, footnote 14). The renter-to-owner budget charges the buyer's
  transaction cost, the down payment, the current mortgage payment, maintenance,
  and current resources in the same equation (Eq. (9), p. 693).
* The income process is not an endpoint-only construction. The authors start
  from the Storesletten--Telmer--Yaron persistent-plus-transitory specification,
  discretize it with Tauchen's method into five states, and report the state
  values and transition matrix for the three-year horizon (p. 702). The
  transition is therefore a period-frequency object; the paper does not tell
  the reader to raise an annual persistence to the third power and then append
  an independently averaged transitory shock.
* Entry wealth is explicitly coupled to income. Each income state receives an
  initial asset level chosen to match the nonhousing-wealth-to-earnings ratio
  for the PSID age-20-to-23 cohort (p. 702). This is a direct precedent for
  preserving a wealth marginal while coupling entrant wealth to the observed
  income rank. It is not a justification for redrawing income independently of
  entry wealth.

The paper's author choices are also instructive limitations. It uses a single
current three-year productivity state, not four annual subperiods; it does not
claim to recover the distribution of an ex-post sum of annual earnings. The
entry asset distribution is a calibration device, not a universal law of entry.

### Sommer, Sullivan, and Verbrugge (2013), *Journal of Monetary Economics*

Primary sources: [author-hosted PDF](https://www.kamilasommer.net/RentPriceRatio.pdf),
[journal DOI](https://doi.org/10.1016/j.jmoneco.2013.04.017). The article is
pp. 854–870.

* The model period is one year (p. 857). Young households are born as renters,
  and the model excludes intergenerational transfers of wealth and human
  capital (p. 857). This is a clean housing benchmark for the entry convention,
  but it does not support importing a nonzero entry-wealth law into the project.
* The information set is unambiguous: households alter housing, deposits, and
  collateral debt at the beginning of the period **after observing** the
  within-period labor-income shock (p. 857). The Bellman equation then has the
  current state \(w\), current resources in the budget (Eq. (7)), and an
  expectation over \(w'\) (Eq. (6), pp. 859–860). The down-payment/equity
  constraint is Eq. (9), p. 859.
* This means current income is available for the purchase choice in the model,
  but future income is not. A four-year version that reveals a realized
  four-year total before the purchase would be a different information
  structure, even if its endpoint Markov state had the same persistence.
* The paper's housing market is useful for the project's accounting checks:
  purchase and sale costs, deposits, collateral debt, and the minimum equity
  requirement all appear in the same-period budget and constraint (pp. 858–860).
  It does not model four-year income aggregation and supplies no four-year
  conversion formula.

### De Nardi (2004), *Review of Economic Studies*

Primary sources: [published PDF](https://users.nber.org/~denardim/research/denardi.pdf),
[journal record](https://academic.oup.com/restud/article-abstract/71/3/743/1536436),
and [replication package](https://users.nber.org/~denardim/research/eflife.zip).
The article is pp. 743–768.

* The model period is five years (p. 747). Individuals enter at age 20,
  consume and work, and the current productivity state is in the recursive
  state. The current earnings term enters the budget constraint, while the
  next productivity state is integrated in the continuation value (Eqs. (1)--(3),
  pp. 748–749). This is the same coarse-time information convention as the
  housing papers.
* The measurement is the important precedent. Appendix A.3 says that PSID
  earnings were divided into equally spaced five-year periods, only complete
  five-year spells were retained, and income was aggregated within each cell
  before estimating the log-earnings process (pp. 765–766). This is direct
  estimation of a five-year object, not a conversion of an annual AR(1) by
  persistence exponentiation.
* Entry wealth is zero by construction (p. 747). The paper also transmits a
  parent's productivity state to the child's first productivity draw, while
  physical wealth is handled through bequests. These are author choices for a
  family OLG model, not evidence that every model should use zero entrant
  wealth or a parent-productivity inheritance link.

The downloadable `eflife.zip` contains Fortran source and Markov/data files.
The paper's code is useful for historical replication, but it does not provide
an exact four-year aggregation or a housing-purchase implementation.

### Sommer (2016), *Journal of Monetary Economics*: annual family timing

Primary sources: [journal article PDF](https://www.kamilasommer.net/Fertility.pdf),
[JME record and DOI](https://www.sciencedirect.com/science/article/pii/S0304393216300745).
The published article is pp. 27–38; the accessible author PDF has the same
model equations and calibration section.

This is an annual fertility model rather than a housing model, so it is a
timing cross-check, not a purchase precedent. The model period is one year and
the wage process has a deterministic age profile, a persistent AR(1), and an
iid transitory shock (Eq. (1), p. 30). The current wage enters the budget and
the next wage is integrated in the Bellman expectation (Eqs. (1)--(2), p. 31).
Households start at age 18 with zero assets (p. 30). The calibration uses
annual \((\rho,\sigma_{\mathrm{persistent}},\sigma_{\mathrm{transitory}})
=(0.95,0.21,0.17)\) and discretizes persistent and transitory components
separately (Table 2 and Section 4.1, pp. 31--32). Those numbers are annual
inputs; the paper does not certify a four-year conversion.

## Joint aggregation and information traps

Suppose annual logs follow

\[
 z_{j+1}=\rho z_j+\eta_{j+1}, \qquad
 \log y_j=h_j+z_j+u_j,
\]

with iid transitory \(u_j\). A mathematically exact joint object for a four-year
block is the conditional kernel

\[
 K(dY,dz'\mid z_0)=
 \Pr\!\left(\sum_{j=0}^{3} e^{h_j+z_j+u_j}\in dY,
 z_4\in dz'\mid z_0\right).
\]

That kernel is useful for simulation and moment construction. It does not by
itself determine when \(Y\) or \(z'\) becomes observable. If \(Y\) is an ex-post
four-year total and is sampled into the state before the purchase decision, the
household has been told the future annual flows and can use them to satisfy a
down-payment condition. If \(Y\) is deliberately defined as current
beginning-of-block cash flow, that particular leakage is absent, but drawing
\(z'\) before the decision still reveals the future endpoint persistent state.
Standard coarse-period models instead reveal the current state, choose, and
then draw the next state in the continuation expectation.

There are three distinct timing choices; they must not be conflated:

* **Beginning-of-block decision:** observe \(b_t\), sale proceeds if applicable,
  the current income state, and current prices; choose the house and financing;
  future annual earnings are not known. This is the convention in the cited
  recursive housing models.
* **End-of-block decision:** first realize and accumulate four years of income,
  then purchase. This can use a realized \(Y_t^{(4)}\), but the state is now
  end-of-block wealth, and the model must say what housing services and rent
  households receive during the block.
* **Within-block decisions:** retain annual subperiod assets and housing status.
  This is the exact solution if purchase timing within four years is an object of
  interest; it is not a one-state endpoint approximation.

The proposed joint `(block income, next persistent state)` kernel is therefore
safe only after this timing choice is fixed. If block income is an ex-post total,
drawing it before the decision leaks future flows; if it is current cash flow,
the endpoint state remains the separate leakage. A joint kernel can preserve
the covariance between a realized block total and the next persistent state,
but the total and endpoint must remain inside continuation integration unless
the model explicitly assumes both are observed at purchase. It should not be
promoted to a beginning-of-block observed state by construction.

## Purchase accounting contract needed here

The project should write the following definitions before fitting or comparing
earnings candidates:

1. **Income unit:** state whether \(y_t\) is a four-year total, a four-year
   average, or a one-period cash flow. If it is a total, every budget and tax
   object must use the total; if it is an average, the factor four must appear
   in all flow accounting. Gross versus disposable income, labor income versus
   transfers, and reference-person versus household income must be fixed.
2. **Asset timing:** state whether \(b_t\) is beginning-of-period liquid wealth
   and \(b_{t+1}\) is after current consumption and housing transactions. A
   purchase test such as
   \(b_t+S_t+y_t/R\ge(1-\phi)p_tH_t\) is coherent only if \(y_t\) is current
   within-period resources and the discounting convention is explicit. Adding
   the same \(y_t\) again to the post-purchase budget counts it twice.
3. **Information:** state exactly which income shock is observed at the purchase
   decision. Future annual realizations may affect expected continuation value,
   but cannot increase current down-payment capacity unless the model explicitly
   gives households a credit contract against them.
4. **Entry law:** state the joint law of entrant wealth and the observed income
   state. Chambers et al. provide one empirical joint-calibration choice: asset
   levels vary by income state to match a wealth-to-earnings ratio. Sommer et al.
   and De Nardi show that zero-transfer or zero-wealth entry is also possible as
   an author choice. An independent redraw of income that changes entry wealth
   is a changed initial-population experiment, not a pure earnings-risk test.
5. **Current purchase and service flow:** if purchase, mortgage writing, payment,
   and housing services all occur in the same period, say so. Chambers et al.
   explicitly do this. If the four-year period hides annual rent and saving,
   that hidden within-block accounting must be declared rather than inferred.

## Recommendation and stop condition

The shippable earnings block should not claim that the endpoint approximation is
an exact four-year aggregation. Before any joint refit, choose either direct
four-year cell estimation or an annual-path simulation whose beginning-of-block
information set is preserved. Then pin the income unit, current-shock timing,
purchase eligibility, beginning/end wealth convention, and entrant wealth-income
coupling in the contract. The current-income purchase adapter can be evaluated
as a separately labeled timing experiment, but its \(y/R\) term must not be
combined with a state that already reveals a realized four-year future total.

No source reviewed here supplies an arbitrary four-year parameter vector or
licenses future-income revelation at a house-purchase decision. No model run,
implementation change, target change, or specification adoption follows from
this memo.
