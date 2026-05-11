# Alternative Housing / Finance Natural Experiments

Date: `2026-04-29`

Purpose: after the FHA/GSE loan-limit design produced a visible HMDA first
stage but a negative ACS/natality reduced form, broaden the search for other
policy variation that could identify effects of financial constraints, mortgage
access, housing user costs, or family-sized housing access on housing and
fertility.

## Recommendation

Do **not** spend another round trying to force the 2008-2009 loan-limit design.
The First-Time Homebuyer Tax Credit is now the best short-run U.S. screen in
the current public-data pipeline, but only as housing-mechanism evidence. The
implemented ACS screen shows a positive response in housing space for more
exposed cheaper metros, while ACS birth-in-past-year does not move positively.

The best next moves are:

1. Treat the implemented `2008-2010` First-Time Homebuyer Tax Credit screen as
   supportive evidence for a down-payment / family-sized-housing mechanism, not
   as evidence of a fertility response.
2. Scope a more serious pre-crisis U.S. mortgage-credit design
   using staggered state anti-predatory lending laws from `1999-2004`.
3. If an international appendix is acceptable, treat UK Right to Buy and UK Help
   to Buy as the cleanest ownership-subsidy analogues.

The U.S. tax-credit screen was feasible because the current ACS/HPI workflow
already covers the right window. The anti-predatory lending-law design is more
credible as a research design because it predates the Great Recession and has
staggered state/legal variation, but it requires older HMDA and outcome panels.

## Candidate A: First-Time Homebuyer Tax Credit

### Why It Is Worth Trying

This is a direct liquidity/down-payment intervention for entry into ownership.
The policy is national, but the fixed dollar subsidy creates strong local
heterogeneity:

\[
CreditShare_c = \frac{8000}{P_{c,pre}}.
\]

The policy should bind more in low-price markets and for deposit-constrained
first-time buyers. That maps more directly to the model's owner-entry constraint
than the FHA/GSE high-cost loan-limit design did.

The literature gives a first stage. Hembre's Regional Science and Urban
Economics paper finds that the tax credit increased first-time homebuyer
purchases, with larger effects in lower-home-value markets. The paper estimates
around `399,846` induced first-time homebuyers and reports that a doubling in
average home values reduced the local effect size by `18.8` percentage points.

### Empirical Design

Use:

\[
Y_{gct} = \alpha_{gc} + \lambda_t
  + \sum_{\tau \neq 2007}\beta_{\tau}
    1[t=\tau] \times CreditShare_c \times Eligible_g
  + X_{ct}'\gamma + \varepsilon_{gct}.
\]

Natural eligible groups:

- ACS household heads ages `25-34`;
- renters or nonowners before the policy;
- low-wealth / low-income households if observable;
- natality ages `25-34` as a secondary screen only.

Outcomes:

- ownership;
- rooms and `6+` room units;
- birth in past year in ACS;
- natality birth rates with lags `0-3`.

Controls:

- local HPI relative to `2007`;
- MSA unemployment;
- pre-2007 price growth;
- fixed effects by geography and year.

### Weaknesses

- It still sits on top of the Great Recession.
- Public ACS does not identify first-time buyers directly.
- In high-cost markets, the credit is too small to relax the relevant
  constraint, so estimates may mostly come from cheaper markets where the
  paper's high-cost spatial mechanism is less central.
- A positive fertility effect would be a timing effect unless completed
  fertility is observed.

### Verdict

This is the best **completed quick screen** so far. It is not the cleanest
design, but it is cheap and mechanically closer to the model than the loan-limit
expansion.

Implemented file:

- `code/empirical/mortgage_policy/run_fthc_tax_credit_screen.R`

Outputs:

- `code/empirical/mortgage_policy/output/fthc_exposure_2007.csv`
- `code/empirical/mortgage_policy/output/fthc_acs_screen_panel_2005_2012.csv`
- `code/empirical/mortgage_policy/output/fthc_acs_summary_bins.csv`
- `code/empirical/mortgage_policy/output/fthc_acs_event_study.csv`
- `code/empirical/mortgage_policy/output/fthc_acs_summary.md`

The implemented exposure is:

\[
CreditShare_c = \frac{8000}{MedianOwnerReportedHomeValue_{c,2007}},
\]

computed from ACS owner household heads ages `25-34` by CBSA. The continuous
treatment is scaled per `10` percentage points, and the binary treatment marks
top-quartile exposed CBSAs. The top-quartile cutoff is `0.0711`, so the high
exposure group is roughly the set of metros where the credit was at least
`7.11%` of the 2007 ACS owner-reported median home value for young household
heads.

Main read from the ACS screen:

- The high-exposure group is cheaper by construction: median owner-reported
  home value around `$112,500` versus `$225,000` in the other group.
- Ages `25-34` ownership does not show a strong first-order response. In the
  continuous specification, the post-2007 ownership coefficients are positive
  but small and statistically weak: `0.021` in 2008, `0.012` in 2009, and
  `0.010` in 2010 per `10` percentage points of credit share.
- Housing space moves more clearly. In the continuous specification, mean rooms
  rise by `0.167` in 2008 and `0.154` in 2009 per `10` percentage points of
  credit share; the `6+` rooms share rises by `0.046` in 2008 and `0.035` in
  2009.
- The top-quartile exposure specification gives a persistent positive response
  for `6+` room units: `0.022` in 2008, `0.018` in 2009, `0.029` in 2010,
  `0.027` in 2011, and `0.028` in 2012.
- ACS birth-in-past-year for women ages `25-34` does not move positively:
  top-quartile effects are `0.001` in 2008, `-0.007` in 2009, and `-0.008` in
  2010. The ages `35-44` placebo group is also not perfectly clean.

Interpretation: the FTHC evidence is useful for the model because a
fixed-dollar down-payment subsidy appears to change access to larger housing
units. It should not be sold as direct evidence that relaxing down-payment
constraints raised fertility. It also identifies cheaper housing markets, while
the paper's core mechanism is high housing costs in constrained locations.

## Candidate B: State Anti-Predatory Lending Laws

### Why It Is More Promising Than It First Sounds

States adopted anti-predatory lending laws at different times around
`1999-2004`, beginning with North Carolina in `1999`. These laws changed the
availability and composition of high-cost mortgage credit before the Great
Recession. The St. Louis Fed summary emphasizes that at least `23` states had
laws by the end of `2004`, with substantial variation in coverage and
restrictions.

This is attractive because it is:

- staggered across states;
- pre-crisis;
- directly about mortgage product availability;
- measurable in HMDA;
- supported by public-use replication material at ICPSR for Ho and
  Pennington-Cross.

The first stage is not one-signed by default. Some laws restricted subprime
flows; broader coverage laws may have increased applications by reducing fear of
predation. That is a feature rather than a bug if the treatment is decomposed
into legal components:

\[
MortgageSupply_{cst}
  = \alpha_c + \lambda_t
  + \beta_1 Coverage_{st}
  + \beta_2 Restrictions_{st}
  + X_{ct}'\gamma + \varepsilon_{cst}.
\]

### Empirical Design

Preferred first stage:

- replicate / extend Ho and Pennington-Cross style border-county design;
- use HMDA `1998-2005`;
- test applications, originations, and rejection rates for subprime/high-cost or
  high-rate proxy loans;
- split laws by coverage versus restriction indices.

Second-stage public-data screen:

\[
Y_{gst} = \alpha_{gs} + \lambda_t
  + \sum_{\tau \neq -1}\beta_{\tau}
    1[t-T_s=\tau] \times LawType_s \times Eligible_g
  + X_{st}'\gamma + \varepsilon_{gst}.
\]

Eligible groups:

- young renter / low-income ACS households;
- ages `25-34` natality rates;
- placebo ages `35-44` or older household heads.

Potential outcomes:

- ownership entry;
- rooms / crowding;
- fertility timing;
- state or border-county migration.

### Weaknesses

- Historical HMDA and outcome panels need a new build.
- Law indices must be handled carefully; a single `law on` dummy is too crude.
- State policy adoption may correlate with prior subprime-market abuse, so
  border-county and pretrend checks are essential.
- The mechanism is mortgage-market composition and consumer protection, not a
  clean increase in credit access everywhere.

### Verdict

This is the best **serious U.S. next design** if we want something more credible
than the crisis-era federal policies. It is a medium-sized empirical project,
not a one-hour check.

## Candidate C: UK Right To Buy

### Why It Is Clean Conceptually

Right to Buy gave eligible UK public-housing tenants the right to purchase their
homes at large discounts after its `1980` introduction. Disney, Gathergood,
Machin, and Sandi describe it as a large-scale natural experiment that raised UK
homeownership by more than `10` percentage points between 1980 and the 1990s,
using eligibility variation for identification.

This is very close to a subsidized tenure-transition experiment:

\[
owner\ entry\ price \downarrow
\quad \Rightarrow \quad
ownership \uparrow,\ housing\ investment \uparrow,\ fertility?
\]

### Potential Design

Use local authority exposure:

\[
Exposure_l = PublicHousingShare_{l,pre}
  \times DiscountEligibility_{t}.
\]

Then estimate:

- ownership responses;
- room / crowding responses from UK Census;
- local births by mother's age from ONS;
- long-run family size if accessible.

### Weaknesses

- It is the UK, not the U.S.
- The treated group is social-housing tenants, not marginal private renters.
- Fertility interpretation may mix wealth, tenure security, neighborhood
  composition, and selection out of social housing.

### Verdict

Excellent external-validity appendix if we are willing to go international. It
is probably cleaner for tenure than anything public in the U.S.

## Candidate D: UK Help To Buy

### Why It Maps To The Model

UK Help to Buy relaxed high-LTV and LTI constraints. The National Audit Office
summarizes the equity-loan scheme as giving buyers up to `20%` of the market
value of an eligible new-build property, or `40%` in London from February 2016,
interest-free for five years. IFS describes the equity loan as relaxing both the
deposit constraint and the loan-to-income constraint because the government loan
reduced the mortgage needed from private lenders.

This maps directly to:

\[
a \geq (1-\phi)p_iH - subsidy_{lt},
\]

or to a policy-specific loan-to-income cap relaxation.

### Potential Design

Exploit:

- April `2013` introduction;
- new-build eligibility;
- London `40%` equity loan from February `2016`;
- local pre-policy new-build share;
- local distribution of house prices below the cap.

Main outcomes:

- first-time buyer purchases;
- ownership among young households;
- new-build supply;
- local fertility with lags.

### Weaknesses

- International and recent.
- Completed fertility is unavailable for a long horizon.
- Local supply constraints may turn the subsidy into prices.
- The policy targets new builds, not all family-sized housing.

### Verdict

Good modern ownership-access design, especially for a model appendix comparing
LTV/LTI relief to supply constraints. Less suitable for completed fertility.

## Candidate E: German Rent Brake / Rent Control

### Why It Is Different

Germany's `2015` rent control generated staggered local treatment dates and
directly changed renter user costs. Mense, Michelsen, and Kholodilin use an
event-study / difference-in-differences design and find reduced regulated rents,
increased free-market rents, and lower mobility from rent-controlled areas.

This targets the renter user-cost channel:

\[
r_{it} h \downarrow
\quad \Rightarrow \quad
housing\ services,\ location,\ fertility?
\]

### Weaknesses

- It is not mortgage credit.
- Rent control also creates misallocation and market segmentation.
- Local fertility data may be possible but would require German data work.

### Verdict

Good if the paper wants a renter-cost validation, not if the appendix must be
about owner-entry constraints.

## Candidate F: Property Tax Limitations

California Proposition 13 and similar assessment-limit systems create a user
cost and mobility wedge for owners. Wasi and White show that Prop 13 increased
average tenure length in California relative to comparison states from `1970` to
`2000`, with larger lock-in effects where the subsidy was larger.

This can discipline mobility / housing mismatch:

\[
move\ cost_{owner} \uparrow
\quad \Rightarrow \quad
less\ relocation,\ possible\ mismatch\ with\ family\ size.
\]

Weakness: this is not a credit-constraint design, and fertility effects would be
indirect.

Verdict: useful background for the moving-cost / location block, not the next
first-order fertility design.

## Candidate G: CRA / GSE Affordable Housing Goals

CRA and GSE affordable housing goals give clean tract-level thresholds. Bhutta's
GSE affordable-housing-goals paper uses discontinuous tract eligibility and
changes in eligibility in `2005`, finding a small increase in GSE activity
without apparent crowd-out of FHA and subprime lending.

This is strong for mortgage supply, but public fertility and household outcome
data are too geographically coarse. It becomes interesting only with restricted
tract-level births or linked administrative data.

Verdict: cite as mortgage-supply evidence; do not make it the next public-data
project.

## Candidate H: Housing Voucher Lotteries / Moving To Opportunity

Voucher lotteries and MTO are clean randomized housing-assistance experiments.
Jacob and Ludwig use a randomized Chicago voucher lottery; Fuller et al. study
MTO and teen / young-adult parenting. These are credible housing-policy
experiments, but the mechanism is low-income rental assistance and
neighborhood relocation, not mortgage constraints or family-sized owner access.

Verdict: cite as housing-affordability / location evidence. Not a direct match
to the model's owner-entry margin.

## Implementation Status And Next Build

### Completed Short Screen: First-Time Homebuyer Tax Credit

The quick public-data screen is implemented in:

- `code/empirical/mortgage_policy/run_fthc_tax_credit_screen.R`

It reuses the existing mortgage-policy output architecture and builds local
exposure as:

\[
FTHCExposure_c = 8000 / P_{c,2007},
\]

where \(P_{c,2007}\) is the ACS owner-reported median home value for owner
household heads ages `25-34`.

It estimates ACS event studies for:

- ownership;
- mean rooms;
- `6+` room share;
- birth in past year.

The screen passes the housing-space bar but not the fertility bar. Do not use
the FTHC screen as a headline fertility estimate unless a more credible outcome
panel and stronger placebo behavior are built.

Useful next refinements, if pursuing this further:

- split by pre-policy renter status if the ACS extract can support it cleanly;
- replace ACS birth-in-past-year with natality rates only after auditing the
  birth-rate denominator and geography coverage;
- add pre-2007 price-growth controls and stricter trend balance checks;
- construct a newly-affordable-homes exposure if transaction-level or binned
  home-value data are available.

### Medium Build: Anti-Predatory Lending Laws

1. Download Ho and Pennington-Cross ICPSR public-use data and documentation.
2. Extract law dates and law-component indices.
3. Download / process HMDA `1998-2005`.
4. Reproduce the first-stage law effects on applications, originations, and
   rejection rates.
5. Only then merge to Census/ACS/natality outcomes.

## Sources

- St. Louis Fed, Ho and Pennington-Cross, state predatory lending laws:
  <https://www.stlouisfed.org/publications/regional-economist/january-2006/states-fight-predatory-lending-in-different-ways>
- ICPSR public-use replication data, Ho and Pennington-Cross:
  <https://www.icpsr.umich.edu/web/ICPSR/studies/1342>
- Hembre, first-time homebuyer tax credit:
  <https://www.sciencedirect.com/science/article/abs/pii/S0166046216303738>
- Treasury announcement of the expanded `2009` first-time homebuyer credit:
  <https://home.treasury.gov/news/press-releases/tg39>
- IFS, Help to Buy distributional effects:
  <https://ifs.org.uk/articles/who-benefits-help-buy-schemes>
- NAO, Help to Buy equity-loan progress review:
  <https://www.nao.org.uk/reports/help-to-buy-equity-loan-scheme-progress-review/>
- Disney, Gathergood, Machin, and Sandi, UK Right to Buy:
  <https://eprints.lse.ac.uk/119338/>
- Mense, Michelsen, and Kholodilin, German rent control:
  <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3494242>
- Wasi and White, Proposition 13 lock-in:
  <https://www.nber.org/papers/w11108>
- Bhutta, GSE affordable housing goals:
  <https://www.federalreserve.gov/econres/feds/gse-activity-and-mortgage-supply-in-lower-income-and-minority-neighborhoods-the-effect-of-the-affordable-housing-goals.htm>
- Jacob and Ludwig, housing voucher lottery:
  <https://www.nber.org/papers/w14570>
- Fuller et al., Moving to Opportunity and teen / young-adult parenting:
  <https://www.sciencedirect.com/science/article/pii/S2352827319300606>
- Bulman, Goodman, and Isen, lottery resources and fertility:
  <https://www.nber.org/papers/w30743>
