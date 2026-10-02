# Purchase financing and the calibration plateau: independent assessment

October 1, 2026. Read-only: no model solve, no change to code, targets, weights or status
files. Scripts and outputs are in this folder; identities and limits are in
[README.md](README.md).

Two evidence bases are used, and neither is the verified latest point.

- **The search sample**: the 1,277 exploratory candidates of the live search
  (`normalized_calibration_v2`, best sampled loss 30.08).
- **The 31.28 point**: the earlier selected point (chain 7, case `0173_nm`, loss 31.284, price
  0.7192) whose policy and distribution arrays are saved locally. Same economics, same
  targets; its ten parameters are within a few percent of the latest sample.

Statements about the 31.28 point are not results at the latest estimates.

## Bottom line

1. **The purchase block has no accounting error.** Income enters resources once. I find no
   double counting and no case for an extra income debit.
2. **The eligibility screen is redundant.** Positive consumption and the end-of-period debt
   floor already imply it, with slack of at least 15% of the house price. The operative
   rule is therefore: no test at closing, and 20% equity by the end of the first four-year
   period. That is an advance against the four years of income that follow closing, at the
   riskless rate. Neither Greaney, Parkhomenko and Van Nieuwerburgh nor Sommer, Sullivan
   and Verbrugge does this: in both, the down payment comes from resources in hand at the
   purchase date.
3. **This matters for tenure, not for the three large misses.** At the 31.28 point, 72% of
   purchases by young never-parent renters close with less than 20% down and 31% with
   nothing down. Yet moving purchase financing from 80% to 100% changes early fertility by
   −0.002, the first-birth rooms response by +0.006 and mean rooms by at most +0.05,
   against gaps of −0.28, −0.24 and +0.25.
4. **The plateau is three separate problems.** Early fertility is outside the model's
   support: its ceiling is 0.67 against a target of 0.81. The first-birth rooms response is
   truncated by the six-room rental cap and by who has children in the model. Mean rooms is
   traded against the ownership rows through the price, under weights inherited from older
   targets.
5. **The search has not converged, one bound is active, and the system is effectively
   under-identified by one parameter.** A better search is worth two to four loss points
   and cannot move early fertility. Whether a higher parenthood housing floor closes the
   first-birth rooms gap is open; one profile run answers it.
6. **The intended mechanism fails for a structural reason, not a financing-timing reason.**
   In the model a first birth lowers ownership in the birth period by 3 points; in your
   PSID event study ownership is 16 to 25 points higher one to two years after the birth.
   Parents can rent six rooms at user cost, so a child does not require a purchase, and a
   purchase constraint cannot hold back a birth.

## 1. The purchase equations

Notation follows the request: $b$ is liquid wealth before the transaction, $Q$ the house
price, $y$ the four-year after-tax income of the period, $R=1.02^4$, $\phi=0.8$, and
$\kappa\equiv\delta+\tau_H=0.09785$ the four-year maintenance and property-tax rate.

### 1.1 Feasible set

A renter who buys chooses consumption $c$ and end-of-period wealth $b'$ in

$$\mathcal F(b,y;Q)=\Big\{(c,b'):\ c>0,\quad b'\ge-\phi Q,\quad c+\kappa Q+b'=R\,(b-Q)+y\Big\}.$$

The set is non-empty if and only if

$$b+\frac{y}{R}\;>\;\Big(1-\frac{\phi}{R}+\frac{\kappa}{R}\Big)Q=0.3513\,Q .$$

The screen asks for $b+y/R\ge(1-\phi)Q=0.20\,Q$. The difference between the two thresholds is
$\phi(1-1/R)+\kappa/R=0.1513>0$, so every feasible purchase passes the screen. The screen
never binds on its own.

In end-of-period dollars a buyer must set aside, out of $Rb+y$,

$$(R-\phi+\kappa)\,Q=0.3803\,Q=\underbrace{0.2000\,Q}_{\text{equity}}+\underbrace{0.0824\,Q}_{\text{interest on the whole price}}+\underbrace{0.0978\,Q}_{\text{maintenance and tax}} .$$

Code checked line by line in `code/model/refactor_lab/engine/`: resources
`Rv = Rg * b + yj` (`household.py:298`, with `b` the post-transaction balance for a buyer);
screen `dp_choice = ctx.dp_arr - income_for_purchase` and
`bmo_purchase = max(ctx.bmo - income_for_purchase, b_grid[0])` (`household.py:894-901`),
where `income_for_purchase` is the same `income_at_state(...)` divided by `R`; buyer floor
`ctx.bmo = -phi * Q` with unsecured add-ons set to zero (`household.py:423-437`);
consumption must exceed `1e-10` (`kernels.py:174-176`). The lowest wealth node is −12, below
$-Q$ for every house, so the grid edge does not bind. From the saved policies the identity
$c+\kappa Q+b'-R x$ returns the same income at every wealth node (spread 0.000) for the
five upper income states.

### 1.2 A $100 house

| Wealth at closing $b$ | Debt at closing | Owed at year 4 | Paid out of $y$ by year 4 (debt reduction + $\kappa Q$) | Loan-to-value at closing |
|---:|---:|---:|---:|---:|
| 0 | 100 | 108.24 | 28.24 + 9.78 = 38.03 | 100% |
| 10 | 90 | 97.42 | 17.42 + 9.78 = 27.20 | 90% |
| 20 | 80 | 86.60 | 6.60 + 9.78 = 16.38 | 80% |

With $b=0$ the screen needs only $y\ge21.65$; the floor and budget need $y-c\ge38.03$.

### 1.3 Established, assumed, unvalidated

**Established mathematically.**

- No resource is counted twice. The screen is an inequality test; it creates no resources.
  The earlier claim of double counting used an example in which consumption is negative
  once the floor is applied.
- $x=b-Q$ is the balance immediately after closing. Renters cannot borrow, so $b\ge0$ and
  loan-to-value at closing is $1-b/Q\in[0,1]$.
- Interest is charged: $R$ applies to the whole post-transaction balance for the period.
- In a stationary equilibrium nobody is below $-\phi Q$, so the incumbent rule
  $b'\ge\min(b,-\phi Q)$ reduces to $b'\ge-\phi Q$ and the estate floor $-0.94\,Q$ is slack.
  Both matter only when prices fall or the purchase limit differs from the incumbent limit.

**Economic assumptions.**

- The budget dates the purchase one interest period before income: $Q$ accrues $R$, $y$ does
  not. Income received up to closing is already inside $b$. The $y$ in the screen is income
  of the four years after closing, so the screen is a present-value solvency test and not a
  cash-at-closing test.
- The lender ledger that makes this coherent is an 80% first lien plus a bridge loan of
  $(1-\phi)Q-b$ at the same rate, due at the end of the period. It is riskless inside the
  model because period income is known at the decision and death occurs between periods.
- One rate for saving and borrowing; interest-only after the first period; no
  payment-to-income test; house trades only at period starts.

**Empirically unvalidated.**

- $\phi=0.8$ is an end-of-first-period limit on net debt under this timing, not an
  origination limit. At the 31.28 point only 7% of young buyers end the purchase period at
  the limit; the mean net loan-to-value at that date is 0.53.
- The model's origination loan-to-value among young first-time buyers (mean 0.86; 72% above
  0.80; 31% at 1.00) is untargeted. The CFPB's *Market Snapshot: First-time Homebuyers*
  (March 2020) reports that the median combined loan-to-value of first-time buyers has
  exceeded 80% since 2002; a National Association of Realtors survey release reports a
  typical first-time-buyer down payment of 2% in 2005 and 2006. These are a report and a
  survey release, not estimates matched to the model's sample.
- Whether lenders advance against four years of income.

## 2. Timing alternatives

One dated example: $Q=100$, wealth at closing $b=5$, income $y=60$ over the following four
years, $\kappa Q=9.78$.

| Rule | Test at closing | Position at year 4 | Four-year consumption | Status |
|---|---|---|---:|---|
| Implemented | none that binds ($5+60/R=60.4\ge20$) | owes $95R=102.83$; must end at 80 | 27.39 | baseline |
| Same, stated in $x$, or with the screen deleted, or as "first lien 80 + bridge 15" | identical | identical | 27.39 | equivalent representation |
| Stock test (Greaney et al. eq. 2.2) | $5\ge20$ fails | no purchase this period | — | economic change |
| One-year income test | $5+15/R=18.9\ge20$ fails | no purchase this period | — | economic change |
| Either test with $b=20$ | passes | owes $80R=86.60$; ends at 80 | 43.62 | same as implemented at $b=20$ |
| Floor moved to $-R\phi Q=-86.6$ | none | owes 102.83; may end at 86.6 and stay there | 33.98 | economic change: a permanent 86.6% limit |
| Interest removed in the purchase period | none | owes 95; must end at 80 | 35.22 | economic change: a 7.83 subsidy |
| Income debit | — | income removed that the budget needs to repay the lender | — | not implied by the budget |

**The two papers, in calendar time.**

- Greaney, Parkhomenko and Van Nieuwerburgh (February 16, 2025): at a purchase,
  $b\ge-\phi p h$ after the transaction (eq. 2.2), so a renter needs $b\ge(1-\phi)ph'$ in
  hand; between purchases $\dot b\ge0$ whenever $b\le-\phi ph$ (eq. 2.3); income,
  consumption and interest flow continuously (eq. 2.4). Purchase dates are one year apart
  in their quantitative model. Their $\phi=0.8$ is labelled a standard value.
- Sommer, Sullivan and Verbrugge (2013): annual periods. The purchase $q(h'-h)$, the new
  mortgage $m'$, labour income $w$ and consumption share one budget (eqs. 7–8); the cap is
  $m'\le(1-\theta)qh'$ with $\theta=0.20$ (eqs. 3, 9); interest $(1+r^m)m$ is paid on
  inherited debt. Section 3.3 calls the 20% requirement a proxy for overall credit
  tightness, set conservatively, because loan approval, credit scores and income
  requirements are absent.

Both papers require the down payment from resources in hand at the purchase date. They
differ in what is in hand: nothing beyond accumulated wealth in the first, wealth plus the
income of the purchase year in the second. The second paper's allowance is one year of
income because its period is one year. Carried to a four-year period by the label "current
period income", it becomes four years.

The one-parameter family that nests the cases is a closing test

$$b+\lambda\,\frac{y}{R}\ \ge\ (1-\phi)\,Q,$$

with the budget and the end-of-period floor unchanged. $\lambda=0$ is the stock test,
$\lambda=1/4$ is one year of income, $\lambda=1$ is the implemented screen and is slack. On
the 31.28 point's distribution, the share of purchases by never-parent renters aged 18–42
that would fail is 71.6% at $\lambda=0$, 8.7% at $\lambda=1/4$ and 0% at $\lambda=1$.

Time aggregation cuts both ways. A household that would save for two years and then buy is
moved two years early by the implemented rule and two years late by the stock test. No
four-year convention reproduces it.

## 3. Could the borrowing rule explain the three misses?

Saved stationary cohorts at the 31.28 point, fixed price, all parameters fixed; "buyer" is
the limit applied to new purchases and "stayer" the limit applied to incumbents.

| Moment | Target | 80 / 80 | Buyer 90, stayer 80 | 90 / 90 | Buyer 100, stayer 80 | 100 / 100 |
|---|---:|---:|---:|---:|---:|---:|
| `early_fertility` | 0.8095 | 0.5312 | 0.5306 | 0.5306 | 0.5303 | 0.5291 |
| `first_birth_rooms` | 1.465 | 1.2232 | 1.2290 | 1.2250 | 1.2293 | 1.2295 |
| `mean_rooms` | 5.729 | 5.977 | 5.982 | 6.031 | 5.976 | 6.010 |
| `recent_parent_ownership` | 0.1276 | 0.1184 | 0.1110 | 0.1018 | 0.0752 | 0.0046 |
| `ownership_30_55` | 0.6763 | 0.6570 | 0.6668 | 0.6808 | 0.6875 | 0.7387 |
| Completed fertility | 2.1 | 2.1000 | 2.0977 | 2.0974 | 2.0953 | 2.0897 |

The general-equilibrium roots give the same picture.

**Too-tight financing is rejected as the cause.** Removing the equity requirement altogether
leaves the three rows where they are.

**Too-loose financing is untested.** The stock test would block 72% of young purchases on
impact, so it would move the ownership rows by points. Its effect on the three rows is
unknown. The pattern above runs against it: as the buyer limit tightens from 100% to 80%
with the stayer limit at 80%, first-birth rooms falls by 0.006 and mean rooms is unchanged.

Each miss has its own mechanism.

- **Early fertility.** The observer is children ever born at age 25½ in the 22–25 cell. The
  data moment is 45.7% mothers × 1.77 children per mother. The model has about 44% mothers
  with 1.2 children each. A woman can have at most one birth per four-year period and enters at
  18, so two children by 25 require a first birth at 18–21 and a second at 22–25, and three
  are impossible. With every 18–21 mother having a second birth the observer is 0.674. To
  reach 0.81, 40% of all women would need a first birth at 18–21 (the model has 26%), which
  the mean-age row forbids.
- **First-birth rooms.** After a first birth 57% of new parents rent, and 71% of those are at
  the six-room cap. The households that have a first birth are well housed before it: 5.4
  rooms if the birth does not occur, against 3.8 for the average childless household at
  risk and 4.7 for the PSID baseline. At ages 18–25 the first-birth probability is 0.00–0.01
  for the four lowest income states (40% of the childless) and 0.24–0.78 above.
- **Mean rooms.** The price is set by replacement fertility, and rooms follow the price
  (correlation −0.94 across candidates). A higher child benefit raises the price and lowers
  rooms but also lowers `recent_parent_ownership` and `ownership_30_55`, which carry
  weights of 27,056 and 2,339 against 128 for rooms.

One input discrepancy bears on the rooms rows. The model's earnings profile gives ages 18–25
72% of the working-age mean. The PSID cell means used to define the money unit give 23% at
18–21 and 45% at 22–25.

## 4. Search, identification, bound, measurement, or restriction

Regression of the ten scored moments on the ten parameters across the search sample
($R^2\ge0.99$ for nine rows, 0.94 for first-birth rooms; predicted and actual loss correlate
0.98 out of sample).

| Diagnosis | Evidence | Verdict |
|---|---|---|
| Optimizer | The loss gradient at the best point is not zero. The linear model predicts a reduction of about 2 within ±2 sample standard deviations of each parameter and about 4 within ±5, almost all from mean rooms. Chains are Nelder–Mead in ten dimensions with 35–78 calls, started from six neighbouring points that lie within about 3.5 sample standard deviations of each other. | Not converged; single basin; worth 2–4 loss points, not more |
| Binding bound | `h_P` is exactly 2.3 at the best point and in 25% of candidates; the loss falls as `h_P` rises (−0.28 per sample standard deviation). The bound is the sum of two inherited search intervals and has no documented empirical basis. | Active |
| Weak identification | Singular values of the weighted Jacobian per 1% parameter change: 5.64, 1.09, 0.79, 0.29, 0.19, 0.10, 0.018, 0.013, 0.003, 0.0001. The null direction is `child_benefit_curvature`. | Ten parameters, nine informative rows |
| Unreachable rows | `early_fertility` ranges 0.525–0.536 over all candidates; a five-standard-deviation step moves it 0.015 against a gap of 0.28. `first_birth_rooms` ranges 1.210–1.241; reach 0.034 against 0.24. | Structural for the first; cap and selection for the second |
| Measurement | Early fertility counts teen births and close spacing that the model cannot represent, while the first-birth age targets are already binned to model cells. First-birth rooms depends on the PSID design: 1.03 (household) to 1.47 (head or spouse at baseline). | Incompatible construction for early fertility; design choice for rooms |
| Objective | Weights for mean rooms, first-birth rooms and the ownership rows are inherited from older targets; early fertility uses a working weight of 100. With the documented sampling errors (0.0089, 0.050, 0.028) the three misses are 28, 4.7 and 10 standard errors. | The equal-sized misses are produced by the weights |

The age-25 count is the moment assigned to `child_benefit_curvature` in
`docs/model/calibration_identification_review_20260929.md`. The row does not respond to the
parameter (0.083 per unit, so at most 0.066 over the whole admissible range), and no other
row identifies it. That review found the same weakest direction at an earlier point under a
different scaling, so the result does not depend on the scaling used here.

## 5. Diagnostics

Ordered by cost. The first three need no model solve.

| # | Diagnostic | Hold fixed | Supports | Rejects |
|---|---|---|---|---|
| 1 | Repeat the purchase tabulations (origination loan-to-value; pass rates at $\lambda=0,\tfrac14,1$) on retained arrays of the selected point | everything | the 72% / 9% / 0% pattern carries over | materially different shares |
| 2 | Early-fertility ceiling at the selected point, and the data moment rebuilt on the model's support (at most one birth per four-year window from 18, earlier births assigned to the first window) from PSID birth histories | everything | ceiling below target; rebuilt moment near the model's range | ceiling above 0.81 |
| 3 | Add ownership to the existing matched first-birth branch observer and compare with the PSID ownership event study (+0.16 household design, +0.27 head-or-spouse design at +3/+4 years) | everything | model response near zero or negative: the birth-to-ownership link is missing | model response of similar size |
| 4 | Three policy evaluations at $\lambda=1,\tfrac14,0$, first at fixed price, then with the price re-solved; report all 14 rows, origination loan-to-value and first-birth probability by branch | all ten parameters | "financing timing drives the misses" if the three rows move by a third of their gaps | rows move by hundredths while ownership rows move by points |
| 5 | Bounded Gauss–Newton from the incumbent with a finite-difference Jacobian (11 evaluations per iteration, 3–4 iterations) | targets, weights, bounds | constrained optimum in the mid-to-high 20s with the two rows unmoved | large gains in the two rows |
| 6 | Profile in `h_P` above 2.3 (for example 2.6 and 3.0), other parameters re-optimized by the same method | targets, weights | the bound is the obstacle if the first-birth rooms gap closes without mean rooms rising | the cap is the obstacle if the response saturates |

Recalibration is needed only after a rule is adopted. Diagnostic 4 at fixed parameters
isolates the mechanical effect. If a different closing test is adopted, the ownership rows
shift by points, so the tenure and price parameters must be re-estimated and the change
labelled as economic. Diagnostics 2 and 3 are measurement decisions; a rebuilt target needs
a new contract name and fingerprint.

## 6. Joint or separate

Separate. The financing rule moves the ownership rows and little else. The link between the
two questions is indirect: a change in the closing test shifts `recent_parent_ownership`
and `ownership_30_55`, which shifts the price that balances them against mean rooms.

On the mechanism, every saved experiment gives the same answer. A current-period relaxation
from 80% to 90% raises ownership in the no-birth branch by 1.84 points, mostly through
two- and four-room houses, and in the birth branch by 0.62 points through six- and
eight-room houses; the value of waiting rises by 0.00027 and the value of a birth by
0.00010; first births fall by 0.006 points on a base of 22.2%. Parents cannot use the
two-room house and almost never choose the four-room one. At the median young income a
six-room purchase with no wealth leaves 0.4 units of four-year consumption at 80% and 0.85
at 90%, against about 1.4 for the same household renting. A stationary comparison with
replacement fertility imposed cannot add to this.

A positive response would require a child to make ownership more valuable than renting at
family sizes. The PSID ownership event study is the moment that measures that link. It is
not in the target system, and the cross-sectional `recent_parent_ownership` row is fitted by
selection: its own observer file calls it "not a causal forced-birth response".

## Open items that need an author decision

- The closing test: keep the implemented rule and describe it as an end-of-period equity
  requirement, or adopt a test at closing ($\lambda=\tfrac14$ or 0). The draft's model
  section still describes the income-excluding version.
- Whether `early_fertility` is rebuilt on the model's support, and what identifies
  `child_benefit_curvature` until then.
- The weight matrix.
- The `h_P` bound and the earnings age profile.
