# Mortgage Policy Natural Experiments

Date: `2026-04-28`

Purpose: identify credible policy variation that can discipline the paper's
housing-credit mechanism: how mortgage constraints, housing user costs, and
housing wealth affect tenure, housing space, location, and fertility.

## Bottom Line

The best near-term design is the 2008-2009 FHA / GSE loan-limit expansion.
It is not a perfect regional rollout, but it creates geographically uneven
credit-supply exposure at the county / CBSA level and maps closely to the
model's purchase constraint.

The design should be treated as an empirical appendix / external validation
exercise, not as a live calibration target yet. The live model keeps `phi`
external; this policy variation is useful precisely because it can validate the
economic importance of mortgage access without forcing fertility moments to
internally identify `phi`.

## Model Mapping

In the discrete-time model, `phi` is the financed share, so the owner entry
constraint uses the down-payment threshold:

\[
(1-\phi)p_iH_k.
\]

A policy that expands FHA/GSE eligibility for larger mortgages is not literally
a change in `phi` for all households. It is closer to a state-dependent
relaxation of the feasible mortgage set:

\[
m \leq L_{ct},
\]

where \(L_{ct}\) is the local FHA or conforming loan limit. For households whose
desired mortgage lies between the old and new limits, the policy should raise
the probability of ownership and potentially allow larger units. If the model
counterfactual is ever implemented, the clean mapping is a location-time
eligible ownership wedge or loan-limit cap, not a global recalibration of
`phi`.

## Candidate 1: 2008-2009 FHA / GSE Loan-Limit Expansion

### Why This Is The Lead Candidate

The Housing and Economic Recovery Act of 2008 and the crisis-era stimulus
changes sharply increased federally backed mortgage eligibility in high-cost
areas. HUD's November 2008 release states that, beginning January 1, 2009, FHA
would insure single-family mortgages from a low-cost floor of `$271,050` up to a
high-cost maximum of `$625,500`; the earlier February 2008 stimulus ceiling had
temporarily reached `$729,750`. The same release explains the rule: FHA limits
were set at 115 percent of local median house price, bounded below by 65 percent
and above by 150 percent of the national conforming loan limit.

FHFA's conforming-loan-limit documentation gives the broader GSE structure:
loan limits are county-based, high-cost values can reach 150 percent of the
baseline, and counties in a CBSA can inherit the limit implied by the
highest-cost component county.

This is close to the model mechanism because it changes the purchase financing
frontier for young households in expensive housing markets.

### Empirical Object

Define treatment intensity at county or CBSA level:

\[
Treat_c = \log L^{FHA}_{c,2009} - \log L^{FHA}_{c,2007}
\]

or, better,

\[
Exposure_c =
\frac{\#\{old\ limit < loan\ amount \leq new\ limit\}_{c,pre}}
{\#\{loan\ amount \leq new\ limit\}_{c,pre}}.
\]

The second measure is preferable because it captures how many local transactions
were near the newly eligible range, not just the statutory change.

### First Stage

Use HMDA to test whether treated markets experienced:

- more FHA-insured purchase originations;
- more GSE-conforming purchase originations;
- more lending in the newly eligible loan-amount band;
- stronger effects for young / moderate-income borrowers when applicant-level
  information is available.

The first-stage regression can be a loan-bin event study:

\[
Orig_{bct} = \alpha_{bc} + \lambda_t
  + \sum_{\tau \neq -1}\beta_\tau
    1[t=\tau] \times NewlyEligibleBand_{bc}
  + X_{ct}'\gamma + \varepsilon_{bct}.
\]

This is cleaner than beginning with fertility because it verifies the mortgage
market channel directly.

### Housing And Fertility Outcomes

After the first stage is established, aggregate ACS/IPUMS outcomes by
CBSA-year, age group, tenure status, and income cell:

- ownership among ages `25-34` and `25-39`;
- rooms and family-sized-unit occupancy;
- center/periphery residence if the MMS PUMA bridge is available for those
  years;
- birth in the past year for women of childbearing age;
- number of own children / childlessness as lower-frequency outcomes.

The main reduced-form specification should be a DDD:

\[
Y_{gct} = \alpha_{gc} + \lambda_t
  + \beta(Treat_c \times Post_t \times Eligible_g)
  + X_{ct}'\gamma + \varepsilon_{gct},
\]

where \(Eligible_g\) is a pre-specified group likely to be affected by mortgage
access: young renters, young households near the income range needed to buy,
or women ages `25-39` in high first-time-buyer cells.

### Identification Concerns

This design is promising but not clean by default.

- The treatment is mechanically larger in expensive markets, which were also
  differentially exposed to the housing boom-bust.
- FHA/GSE policy arrived during the Great Recession, so local unemployment,
  foreclosure distress, and credit-supply collapse are first-order confounds.
- County/CBSA loan limits are partly based on local median prices; this makes
  treatment exposure endogenous to pre-policy housing-market tightness.
- Fertility should respond with a lag, while mortgage originations respond
  immediately.

Required checks:

- event-study pretrends for `2005-2007`;
- placebo outcomes for older households less likely to be first-time buyers;
- placebo loan-amount bins just below the old limit;
- controls for local HPI, unemployment, foreclosure rates, and pre-2007 price
  growth;
- separate estimates for owners versus renters, consistent with Dettling and
  Kearney and Lovenheim and Mumford.

## Candidate 2: HARP Refinancing

HARP is strong for an incumbent-owner cash-flow channel, not for the
first-time-buyer entry constraint. FHFA describes HARP as a refinancing program
introduced in early 2009 for borrowers who could not otherwise refinance because
of home value declines or mortgage-insurance constraints. Eligibility required
that the existing mortgage had been sold to Fannie Mae or Freddie Mac on or
before May 31, 2009, with LTV above 80 percent and current payments.

The literature already shows strong first-stage and balance-sheet effects.
Abel and Fuster use quasi-random access around the HARP cutoff and find that
refinancing reduced defaults and changed consumer debt use. Agarwal et al. show
that HARP relaxed equity constraints, generated average annual savings around
`$3,000`, increased spending, and affected local foreclosures and house prices.

Model mapping:

\[
payment_{it} \downarrow \quad \Rightarrow \quad disposable\ resources_{it}
\uparrow
\]

for owners with fixed mortgages. This maps to user cost / cash-flow and possibly
home-equity extraction, not to the young renter purchase constraint.

Use case for this project:

- secondary empirical appendix on whether owner cash-flow relief accelerates
  births;
- regional exposure design using pre-HARP share of underwater GSE loans, if the
  data can be obtained;
- borrower-level design only if linked mortgage-credit-birth data are available,
  which is unlikely in the current public-data workflow.

## Candidate 3: First-Time Homebuyer Tax Credit

The 2008-2010 first-time homebuyer credit is directly targeted at entry into
ownership. GAO summarizes the policy: the 2008 credit was up to `$7,500` and
repayable over 15 years; the 2009 Recovery Act expanded it to up to `$8,000`
with no repayment unless the home stopped being the taxpayer's principal
residence within three years. The 2009 version was enacted on February 17, 2009
and made retroactive to purchases beginning January 1, 2009.

This policy has a natural heterogeneous exposure:

\[
CreditShare_c = 8000 / median\_home\_price_c.
\]

It should matter more in cheaper markets because the dollar subsidy is fixed.
That is useful, but it is less clean for the current center/periphery mechanism:
the most treated places may be cheaper, less constrained markets rather than
high-cost markets where fertility is suppressed by housing costs.

Use case:

- robustness design for down-payment / liquidity at entry;
- DDD with young renters versus older owners, high versus low credit-share
  markets, before versus after 2009;
- fertility lags of one to three years.

Main problem: the policy is national and sits exactly on top of the Great
Recession, so identification needs strong group-level controls and pretrends.

## Candidate 4: GSE Affordable Housing Goals / Underserved Areas

The GSE affordable housing goals created discontinuous tract-level eligibility
rules and changed tract eligibility in 2005. Bhutta studies this design and
finds a small effect on GSE activity in lower-income and minority neighborhoods,
without apparent crowd-out of FHA and subprime lending.

This is clean for mortgage-supply identification but awkward for this project
unless restricted geocoded fertility or household data are available. Public
ACS/IPUMS PUMA geography is likely too coarse for a tract discontinuity.

Use case:

- literature support for mortgage-supply variation;
- possible future restricted-data design;
- not the first public-data implementation.

## Candidate 5: State Bank Branching Deregulation

State branching and interstate banking deregulation provide classic staggered
credit-supply variation. Jayaratne and Strahan use relaxation of bank branch
restrictions to study finance and growth, and later work connects interstate
branching to mortgage credit and homeownership.

This is less attractive for the current paper because it is broad credit-market
variation: it affects wages, local growth, business formation, and mortgage
supply simultaneously. It is useful background but not a tight housing-fertility
design.

## Literature That Matters For Interpretation

Dettling and Kearney's house-price work is the conceptual benchmark. They argue
that higher house prices create a negative price effect for nonowners but a
positive wealth effect for owners. Their estimates imply that short-run house
price increases reduce births among nonowners and raise births among owners.

Lovenheim and Mumford use housing-market wealth variation in PSID and find that
housing wealth increases fertility for homeowners, with little evidence of an
effect for renters.

Dettling and Kearney's 2025 mortgage-history paper is especially relevant for
this project: they argue that low-down-payment, long-term, fixed-rate FHA and VA
mortgages increased homeownership for young families and led to roughly 3
million additional births from 1935 to 1957. That paper supports the exact
mechanism we care about: access to ownership affects fertility.

Bulman, Goodman, and Isen's lottery-win paper is useful discipline. Financial
resources have large and persistent effects on homeownership and marriage, but
fertility is mostly accelerated rather than raised permanently. This warns
against over-interpreting any short-run birth response as a completed-fertility
effect.

## Practical Implementation Plan

### Stage 0: Feasibility And First Stage

Build a county / CBSA panel for `2005-2012`:

- FHA loan limits by county;
- conforming loan limits by county;
- HMDA purchase originations by loan type, amount bin, county/tract;
- FHFA or Zillow local house-price controls;
- local unemployment / foreclosure controls if available.

Produce:

- treatment maps;
- loan-limit exposure tables for the MMS large-metro sample;
- event-study plots for FHA/GSE originations in newly eligible loan-amount
  bands.

Stop here unless the first stage is visible.

### Stage 1: ACS Outcomes

Aggregate ACS/IPUMS outcomes for the MMS large-metro sample:

- young ownership;
- rooms by tenure and child status;
- center/periphery residence where the PUMA-MMS bridge is available;
- birth-in-past-year for women ages `20-44`, especially `25-39`.

Use group-level DDD designs so national recession shocks are absorbed as much as
possible.

### Stage 2: Fertility Validation

If county-age birth data are available:

- use CDC Natality / county-year birth rates by mother's age;
- estimate lags of `0-4` years after the loan-limit shock;
- compare births among age groups most likely to enter ownership against older
  placebo groups.

The core fertility interpretation should be:

\[
credit\ access \rightarrow ownership / rooms \rightarrow timing\ of\ births,
\]

not necessarily completed fertility.

### Stage 3: Model Connection

Use any credible estimate as an external validation target:

- does a mortgage-access relaxation raise ownership for young households?
- does it raise family-sized housing consumption?
- does it produce a delayed birth response?
- are effects larger for renters than existing owners?

Do not compare this reduced-form response directly to the current SMM loss. It
is a mechanism check, not a calibration objective.

## Implemented Candidate 1 Pipeline

The first-stage implementation now exists under:

`code/empirical/mortgage_policy/`

Candidate 1 has now been implemented through the HMDA first stage, the ACS
housing-outcome screen, and a preliminary natality merge scaffold. The GSE
channel is the first-stage object. FHA 2009 versus 2008 is too weak for this
design because the 2008 stimulus expansion had already moved FHA limits.

Use:

\[
\log L^{GSE}_{c,2009} - \log L^{GSE,proxy}_{c,2007}
\]

as the exposure measure, and track the HMDA first stage:

\[
1\{loan\ amount \in (L^{GSE,proxy}_{c,2007},L^{GSE}_{c,2009}]
  \text{ and sold to Fannie/Freddie}\}.
\]

This is a mortgage-market validation step, not yet a causal fertility estimate.

## Local Data Check

A quick initial workspace check found natality / birth files already available
under `/Users/tommasodesanto/Desktop/Projects/Datasets/natality_data`,
including `2005-2012` files and hybrid MSA/HPI panels. No obvious HMDA,
FHA loan-limit, conforming-loan-limit, or CLL files were already present under
the datasets root.

The first-stage implementation therefore downloaded and processed:

- HUD county FHA/GSE loan-limit files for `2008` and `2009`;
- CFPB HMDA first-lien owner-occupied 1-4 family LAR files for `2007`,
  `2008`, and `2009`;
- HMDA dictionaries for the labeled file schema.

Generated files are under:

`code/empirical/mortgage_policy/output/`

They are intentionally ignored by git because the raw HMDA files are large.

## First-Stage Result From The Implemented Pipeline

The initial HMDA first-stage diagnostic is encouraging. In counties with
above-median positive GSE loan-limit exposure, conventional purchase loans sold
to Fannie/Freddie in the newly eligible GSE loan-amount band rose from `0.48%`
of originations in `2007` to `2.73%` in `2008` and `3.84%` in `2009`.

This is the first-stage object to preserve. Broad counts of all loans in the new
band are contaminated by jumbo / portfolio lending that existed before the
policy. Product-specific GSE-sold loans in the newly eligible band are the
cleaner channel.

Implemented extra scaffold:

- `build_acs_gse_exposure_panel.R` builds ACS housing outcomes for household
  heads ages `20-54`, with the main screen for ages `25-34`;
- `build_natality_exposure_panel.R` merges the GSE exposure panel to the
  collapsed MSA/HPI/pop natality panel for `2005-2012`;
- `run_refined_acs_takeup_design.R` narrows treatment to positive-exposure
  CBSAs with strong HMDA takeup in the newly eligible GSE-sold band, then adds
  local HPI and unemployment controls;
- `run_refined_natality_takeup_design.R` applies the same HMDA takeup screen to
  the natality scaffold;
- the first run matched `230` of `238` natality geographies in the target age
  groups;
- the generated birth-rate summaries should be treated as a data-merge check,
  not an estimate. The panel construction needs a separate audit before any
  fertility-response claim.

## ACS Outcome Result From The Implemented Pipeline

The ACS housing-outcome screen does **not** currently support a simple positive
"loan-limit expansion raised young housing consumption" story.

For household heads ages `25-34`, the clustered metro/year FE diagnostic finds:

- ownership-rate interactions are small and not robustly positive after `2008`;
- mean rooms move negative after `2008`, with the `2008` coefficient about
  `-0.228` per one-log-point GSE exposure;
- the share in units with `6+` rooms also moves negative after `2008`, with the
  `2008` coefficient about `-0.053` per one-log-point GSE exposure.

These are screening coefficients, not final causal estimates. But the sign is
important. The HMDA first stage is visible, yet the broad ACS housing response
does not look like a clean expansion in young household housing quantity.

Interpretation:

- the crisis-period housing-market collapse may dominate the loan-limit
  relaxation in ACS aggregates;
- high-exposure places are high-cost boom-bust markets, not quasi-randomly
  treated regions;
- the policy may affect the composition of mortgages without increasing average
  housing space for young households;
- fertility regressions should not be run as a headline mechanism until this
  housing-response miss is explained or the design is tightened.

The natural refinement is to restrict the treatment to CBSAs with clear HMDA
GSE-sold-band takeup, estimate ACS outcomes on that narrower high-cost metro
sample, and add local HPI and unemployment controls. That refinement has now
been implemented.

## Refined HMDA-Takeup Outcome Result

The refined screen is:

\[
HighTakeup_c =
1\{\Delta Share^{GSEsold,newband}_{c,2008/09-2007}
\geq p75(\Delta Share^{GSEsold,newband}_{c,2008/09-2007}
  \mid Exposure_c>0)\}.
\]

The implemented cutoff is `0.0147`: the post-`2008`/`2009` average GSE-sold
new-band origination share minus its `2007` value. The HMDA file has `73`
positive-exposure metros and `19` high-takeup metros. The ACS regressions use
metro and year fixed effects, ACS household weights, local HPI relative to
`2007`, MSA unemployment, and metro-clustered standard errors.

This refinement does **not** rescue the positive housing-consumption mechanism.
Among household heads ages `25-34`, the high-takeup binary estimates within
positive-exposure metros are:

- `2008` owner rate: `0.004`, `p=0.718`;
- `2009` owner rate: `-0.017`, `p=0.275`;
- `2008` mean rooms: `-0.035`, `p=0.521`;
- `2009` mean rooms: `-0.173`, `p=0.021`;
- `2008` share in `6+` rooms: `-0.005`, `p=0.605`;
- `2009` share in `6+` rooms: `-0.036`, `p=0.009`.

The continuous-takeup specification is also mostly negative in the housing
margins, with a `2009` owner-rate coefficient of `-0.387` (`p=0.036`) and a
`2009` `6+` room-share coefficient of `-0.533` (`p<0.001`) per unit of the
takeup-delta measure. Because the treatment variable is a share change, this
continuous scale is large; the sign and instability are more informative than
the raw magnitude.

The refined natality screen should **not** be used as a fertility result. The
same high-takeup binary specification for ages `25-34` gives a positive but
imprecise `2008` birth-rate coefficient of `35.188` births per `1,000`
population (`p=0.306`), then small/noisy estimates afterward. The continuous
takeup measure has a `2008` coefficient of `739.037` (`p=0.043`) but no stable
post-pattern. More importantly, placebo age groups behave badly: the older
`35-44` group shows significant positive estimates in some post years. That is
not a credible first-time-buyer fertility response.

Recommendation: keep the FHA/GSE loan-limit work as a documented negative
screen / empirical appendix. The HMDA first stage is useful and publishable as
diagnostic evidence that the policy moved GSE-eligible mortgage composition, but
the current public ACS/natality aggregation should not be used as main causal
support for the paper's housing-space or fertility mechanism. To turn this into
a credible causal design, the next step would need finer outcome data or a much
tighter borrower/geography design, not another broad MSA-level fertility
regression.

## Sources

- HUD, 2008 FHA loan-limit announcement:
  <https://archives.hud.gov/news/2008/pr08-174.cfm>
- FHFA conforming loan-limit FAQ:
  <https://www.fhfa.gov/sites/default/files/2024-04/FHFA-CLL-FAQs.pdf>
- FHFA HARP fact sheet:
  <https://www.fhfa.gov/news/fact-sheet/home-affordable-refinance-program-harp>
- GAO first-time homebuyer credit testimony:
  <https://www.gao.gov/products/gao-10-166t>
- Bhutta, GSE affordable housing goals:
  <https://www.federalreserve.gov/econres/feds/gse-activity-and-mortgage-supply-in-lower-income-and-minority-neighborhoods-the-effect-of-the-affordable-housing-goals.htm>
- Jayaratne and Strahan, bank branch deregulation:
  <https://www.newyorkfed.org/research/staff_reports/research_papers/9513.html>
- Dettling and Kearney, house prices and birth rates:
  <https://www.nber.org/papers/w17485>
- Dettling and Kearney, modern mortgage and baby boom:
  <https://www.nber.org/papers/w33446>
- Lovenheim and Mumford, housing wealth and fertility:
  <https://econpapers.repec.org/RePEc:tpr:restat:v:95:y:2013:i:2:p:464-475>
- Abel and Fuster, HARP refinancing:
  <https://www.aeaweb.org/articles?id=10.1257/mac.20180116>
- Agarwal et al., HARP refinancing and spending:
  <https://academic.oup.com/restud/article/90/2/499/6619572>
- Bulman, Goodman, and Isen, financial resources and fertility:
  <https://www.nber.org/papers/w30743>
