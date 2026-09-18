# How Many Young US Households Are Down-Payment Constrained?

Evidence memo, September 18, 2026. Prepared to give the model's 4.5 percent
marginal constrained share an external benchmark.

> **Lead review (September 18, 13:30).** Checked against the papers: Haurin, Hendershott and Wachter (NBER w5630) is the NLSY panel of youth aged 20 to 33 for 1985 to 1990, 37 percent constrained, ownership probability 0.20 to 0.10 and 0.52 to 0.34 when constrained, all verified in the text. Kaplan, Mitman and Violante (NBER w23694) verified: "very few households are constrained in this way: rather than buying excessively small houses, they prefer to rent a house of the desired size", and the credit relaxation raises home ownership by about 3 percent (the working paper attributes the 3 percent to the productivity shock relaxing the payment-to-income limit and says the credit relaxation has "a similar size effect"). The open item in section (f) is resolved: the model's 4.5 percent is a stock object. It is the mass of households aged 26 to 38 with zero or one child at home, in the baseline stationary distribution, whose modal tenure-and-size choice differs between the solution with the down payment and the solution without it at the same pre-choice state (`code/model/sandbox/diagnostics_constrained.py`, "Binds share"). Every household re-chooses housing each period, so the stock comparison to the Haurin range is the right one. Caveat that stays: the model removes the down payment entirely, while every empirical benchmark compares two positive requirements, so 4.5 percent is an upper bound on the model's marginal share under any realistic relaxation.

## (a) The question

The model's credit friction is a single object. A household that buys must put
down \(1 - \phi\) of the purchase price in cash, with \(\phi = 0.80\), so the
required down payment is 20 percent of the price. There is no credit score, no
debt-to-income limit, and no low-down-payment government product. At the current
calibration, removing this requirement changes the housing choice of about
4.5 percent of family-forming households, defined as ages 26 to 38 with zero or
one child at home. The author wants to know whether 4.5 percent is a plausible
number, and what data object it should be compared against.

The 4.5 percent figure was supplied by the author. I did not find the experiment
that produces it anywhere in `code/model/`, so I take it as given and treat the
question as one of external benchmarking, not replication. One thing does need
checking on the model side before any comparison is final, and I flag it in
section (f): whether 4.5 percent is measured over the stock of households in the
relevant age and child cells, or over the flow of households actually making a
move in that period. The empirical literature is split between exactly those two
populations, and the numbers differ by roughly a factor of three.

## (b) Four different objects that all get called "the constrained share"

These are routinely conflated in the housing literature and in referee reports.
They are not the same number and they do not answer the same question.

**(i) Share who could not meet a 20 percent down payment today.** A pure
affordability accounting object. Compare a household's liquid and semi-liquid
assets to 0.20 times the price of the house it would want. No behavior is
involved. This is large, on the order of two-thirds of young renters and recent
movers, and it is the number the model's \(\phi\) most directly corresponds to
in levels.

**(ii) Share who could not meet the actual institutional minimum.** The
US minimum is 3.5 percent for a Federal Housing Administration loan and as low as
3 percent for some conventional products, not 20 percent. This share is much
smaller than (i). The model has no such product, so (ii) has no model
counterpart, and an empirical estimate of (ii) is not a valid target for this
model.

**(iii) Share whose housing or tenure choice actually changes when the constraint
is relaxed.** The marginal group. Most households who fail (i) do not become
owners when you relax the constraint, because they did not want to own at that
price, or because they would rather rent a unit of the size they want than buy a
small one. This is what the model's 4.5 percent measures, and it is the only one
of the four that is comparable to it.

**(iv) Share of actual buyers who received family help with the down payment.**
A conditional-on-purchase object, informative about how binding the constraint is
for people who cleared it, silent about the people who did not. It is evidence
that the constraint binds, not a measure of how many it binds for.

Objects (i) and (iii) are related by a behavioral response: (iii) equals (i)
multiplied by the probability that a constrained household flips when freed. The
literature below gives both pieces, so the product can be formed.

## (c) Literature on borrowing-constrained buyers

| Paper | Sample | "Constrained" definition | Share constrained (i) | Marginal response (iii) |
|---|---|---|---|---|
| Linneman and Wachter (1989) | SCF recent movers, 1981-83 | Desired house value above the maximum affordable given wealth at the maximum loan-to-value ratio; separately given income | 40% highly wealth constrained; 27% highly income constrained | Lower propensity to own; not verified numerically |
| Haurin, Hendershott and Wachter (1996/1997) | Panel of youth aged 20 to 33, 1985-90 | Same, but the household picks the loan-to-value ratio that minimizes the more binding of the two constraints, capped at 0.95 | 58% wealth constrained at a fixed 0.80 loan-to-value ratio; **37% highly constrained by either wealth or income** after optimal loan-to-value choice; 6% moderately | Ownership probability falls from 0.20 to about 0.10 for a low-wage couple, and from 0.52 to 0.34 for a higher-wage couple, so **10 to 18 percentage points** |
| Barakova, Bostic, Calem and Wachter (2003) | SCF 1989, 1995, 1998; recent movers aged 21 to 50 | Wealth: cannot meet 10 percent down on the preferred house. Income: debt payments above 38 percent of income. Credit: imputed score below 620 | 1998: **65.2% wealth, 29.5% income, 35.9% credit**. 1989: 67.0%, 49.3%, 21.1% | 1998 predicted ownership rises 0.29 to 0.48 when the wealth constraint alone is lifted (**+19 pp**), and to 0.56 when all three are lifted (**+27 pp**) |
| Acolin, Bricker, Calem and Wachter (2016) | SCF 2001, 2004-07, 2010-13; all households | Wealth, income and credit as in Barakova et al., combined into one indicator | Not reported separately in the published note | Being constrained cuts the likelihood of owning by 26% (2001), 23% (2004-07), 30% (2010-13). Aggregate ownership in 2010-13 is **2.3 pp** below the 2001 credit regime and **5.2 pp** below the 2004-07 regime |
| Fuster and Zafar (2021) | New York Fed Survey of Consumer Expectations special module, February 2014, 962 respondents (698 owners, 264 renters) | Behavioral: does stated willingness to pay rise when the required down payment falls from 20 percent to 5 percent | 59% would choose to put down less than 20 percent; 78% of renters would | **43% of all respondents and 58% of renters raise willingness to pay.** Average willingness to pay rises about 15 percent overall and about 40 percent for renters. Those who choose 20 percent or more show no change |
| Kaplan, Mitman and Violante (2020) | Calibrated lifecycle model with rental sector, long-term defaultable mortgages, maximum loan-to-value ratio 0.95 pre-boom | Model-internal | Not reported | A pure credit relaxation raises home ownership by about **3 percent**. They state that "very few households are constrained in this way" because a household that cannot buy the size it wants rents that size instead |
| Greenwald (2018) | Fannie Mae loan-level data, 2006 and 2014 | At the loan-to-value limit versus at the payment-to-income limit | Qualitative: a loan-to-value-constrained majority of borrowers, a payment-to-income-constrained minority | Not a share-of-households object |
| DeFusco, Johnson and Mondragon (2020) | CoreLogic loans 2010-15, jumbo versus conforming, around the 43 percent debt-to-income Qualified Mortgage threshold | Above the 43 percent debt-to-income cutoff | Not applicable | The rule **eliminated 15 percent** of the affected market and **shifted another 20 percent** below the threshold, so 35 percent of the affected segment responded |
| Bhutta and Ringo (2021) | Rate-lock and Home Mortgage Disclosure Act data, the January 2015 50 basis point cut in FHA mortgage insurance premiums | Households likely to use FHA (credit score below 680, loan-to-value above 80) | Not applicable | Purchase originations to this group rose about **14 percent**; up to 40 percent of that rise came from fewer denials; no response at all in the top income quartile |
| Engelhardt and Mayer (1994) | Chicago Title and Trust buyer surveys, 1976-82 and 1992, 18 metropolitan areas | Received family help toward the down payment | Object (iv): **about 20 percent** of first-time buyers got help from relatives; 4 percent financed the whole down payment that way; help averaged 50 percent of the down payment when received | Repeat buyers get only 2 percent of down-payment funds from relatives, which the authors read as evidence that gifts are targeted at constrained first-time buyers |

Two entries deserve emphasis because they pull in opposite directions.

Haurin, Hendershott and Wachter is the closest match to the model's population,
being young households only, and it is the only source that gives both pieces of
the product. Multiplying the 37 percent constrained share by the 10 to 18
percentage point ownership response gives a marginal share of **3.7 to 6.7
percent of young households**. That bracket contains 4.5 percent.

Kaplan, Mitman and Violante is the closest match to the model's structure. Their
reason for a small number is exactly the mechanism in this model: with a rental
market offering the same housing services, a household priced out of owning the
size it wants rents that size rather than buying a smaller unit, so the
down-payment requirement changes the tenure portfolio rather than housing
consumption. Their credit relaxation moves about 3 percent of households.

Barakova et al. is the outlier at 19 percentage points, and the reason is the
population. "Recent movers" are the flow of households transacting within roughly
two years, a group with a far higher baseline probability of buying than the
stock of households of the same age. If roughly a third of households aged 26 to
38 move within a two-year window, a 19 percentage point effect on movers is
about 6 percent of the stock. I did not verify that mobility rate against a
primary source, so treat the conversion as an order-of-magnitude adjustment
rather than a number.

## (d) Survey evidence

| Statistic | Value | Source and year | Sample |
|---|---|---|---|
| Families turned down for credit in the past 12 months | 10.1% (2022), 10.7% (2019) | Federal Reserve Bulletin, *Changes in U.S. Family Finances from 2019 to 2022*, Table 5 | All families, 2022 SCF |
| Did not apply for credit for fear of denial | 12.9% (2022), 12.7% (2019) | Same, Table 5 | All families |
| Either turned down or feared denial | 18.4% (2022), 18.4% (2019) | Same, Table 5 | All families |
| Same, by age of head | **Not published.** The Bulletin reports this series only economy-wide | Fed Bulletin | Would require SCF microdata |
| Renters who prefer a down payment below 20 percent | 78% | Fuster and Zafar, February 2014 module | 264 renters |
| Median chosen down payment | 9.1% for renters, 18.4% for owners | Same | Same |
| Renters raising willingness to pay when the requirement falls to 5 percent | 58% | Same | Same |
| Renters' mean stated probability of ever buying a home | 33.9% in the 2025 wave, down from 40.1% in 2024 | New York Fed SCE Housing Survey chart packet | 308 renters |
| Renters saying a mortgage would be very or somewhat difficult to get | 66.8% | Same, 2024 wave | Renters |
| First-time buyers using a gift or loan from family or friends for the down payment | 22% (19% gift, 3% loan) | National Association of Realtors, *2025 Profile of Home Buyers and Sellers*, Exhibit 5-4 | First-time buyers who made a down payment |
| First-time buyers naming "saving for the down payment" as the hardest step | 31% (11% among all buyers) | Same, Exhibit 3-7 | First-time buyers |
| Median down payment, first-time buyers | 10% (2025); 19% all buyers; 23% repeat buyers | Same, Exhibit 5-3 | Buyers, 2025 |
| "Mortgage-ready" renters aged 40 and under | 25% to 45% across 31 metropolitan areas, average 34%. Definition: no current mortgage, credit score at or above 620, debt-to-income at or below 25, no recent foreclosure, bankruptcy or serious delinquency | Urban Institute and Freddie Mac, *Barriers to Accessing Homeownership*, 2018 report on credit-bureau data as of September 2016 | Non-mortgage-holders aged 40 and under |
| Renters naming the down payment as a barrier | 68% | Zillow Housing Aspirations Report 2018, relayed by Urban Institute | Renters |
| Median net worth, renters and other non-owners | $10,400 | Fed Bulletin, 2022 SCF | All renters, not split by age |

Three gaps are worth stating plainly. The New York Fed's current annual housing
survey does not carry a perceived-required-down-payment question or a
main-obstacle question, so the only clean down-payment-sensitivity item in that
series is the 2014 module that became Fuster and Zafar. The Fannie Mae National
Housing Survey could not be opened. And there is no published tabulation of the
share of renters aged 25 to 40 whose liquid assets fall below 20 percent of the
local median house price. That specific object would have to be built from SCF
microdata. It is also, as section (f) argues, the moment worth building.

## (e) Reduced-form fertility papers and the share they move

| Paper | Design | Effect | Implied share moved |
|---|---|---|---|
| Dettling and Kearney (2014) | Vital Statistics births 1990-2007 by metropolitan area, house price index interacted with the local ownership rate | A $10,000 price rise raises owner fertility by 5 percent and cuts non-owner fertility by 2.4 percent | Net **+0.8 percent** of current-period fertility at the mean ownership rate. The 1997-2006 boom implies about a 9 percent rise in births |
| Lovenheim and Mumford (2013) | Panel Study of Income Dynamics, women 25 to 44, simulated metropolitan house price growth | A $100,000 housing wealth gain raises the annual probability of a birth by 0.82 to 0.89 percentage points, that is 16 to 18 percent | Against a 5.0 percent baseline annual birth probability, the 1999-2005 boom moved **0.43 to 0.64 percentage points** of owner women per year. No detectable effect on renters |
| Daysal, Lovenheim, Siersbæk and Wasser (2021) | Danish population registers, owner women 20 to 44 | 100,000 Danish kroner, about $12,000, raises the birth probability by 0.27 percentage points, that is 2.32 percent | Renter effect is -0.04 percent and statistically indistinguishable from zero |
| Hacamo (2021) | 2000 Census and 2005-18 American Community Survey, plus PSID. Triple difference from the January 2004 ruling letting nationally chartered banks bypass state antipredatory lending laws | Fully exposed households are 13 to 15 percentage points more likely to move into a new home with a mortgage, and **6 percentage points** more likely to buy with a mortgage and have a child. PSID replication gives 8.9 | The 6 percentage points is a credit-access marginal share on a young, first-time-buyer-eligible population |
| Kearney and Wilson (2018) | County shale production instrumenting local male earnings | An extra $1,000 per capita of simulated production raises births by 5.96 per 1,000 women aged 18 to 34; about 3 percent at peak intensity | The authors state explicitly that this is reduced form and not an income elasticity |
| Clark and Ferrer (2016/2019) | Canadian Survey of Labour and Income Dynamics, non-movers aged 18 to 40 | Owner effects are small and fragile; the renter fixed-effects estimate is negative and insignificant | Not a usable share |
| Atalay, Li and Whelan (2021) | Australian HILDA panel | A $100,000 price rise raises the probability of a child by 7.5 percent among owners | **Not verified.** Full text could not be opened |

The single most useful number here is Hacamo's 6 percentage points. It is the
only reduced-form fertility estimate whose treatment is a credit shock rather
than a wealth shock, and it is expressed as a share of households that flip a
joint buy-and-have-a-child outcome. Its population is young, exposed, and
first-time-buyer eligible, which is close to the model's family-forming cell.
That said, its treatment is a large, one-sided credit expansion in the early
2000s, not the removal of a 20 percent down payment, and the outcome is joint
rather than housing alone.

The wealth-shock papers are not measuring this model's mechanism. They identify
the effect of house price changes on people who already own. Lovenheim and
Mumford and Daysal et al. both find no renter effect at all, which in this
model's language means the down-payment channel does not show up as a fertility
response of renters to prices in their data. That is a caution about what the
model should be asked to reproduce, not evidence against the mechanism.

## (f) Bottom line and the moment to discipline it

The best defensible range for the marginal constrained share among young
family-forming households is **4 to 8 percent per period**, built from the one
source that is both about young households and reports the two pieces of the
product, Haurin, Hendershott and Wachter, whose 37 percent constrained share
times a 10 to 18 percentage point ownership response gives 3.7 to 6.7 percent,
and corroborated at the low end by Kaplan, Mitman and Violante's 3 percent and
at the high end by Barakova et al.'s 19 percentage points on movers, which
converts to roughly 6 percent on the stock. The model's 4.5 percent sits inside
that range and is, if anything, on the low side, because the model's experiment
removes the down payment entirely while every empirical benchmark listed above
compares two positive down payment regimes. Two structural features defend the
low number: the model has a rental sector offering the same housing services,
which is precisely why Kaplan, Mitman and Violante find few quantity-constrained
households, and the model's \(\phi\) is the only credit friction, so it should
map to the wealth constraint alone and not to the two-thirds-constrained figures
that bundle income and credit. The number that would not be defensible is a
comparison against the 65 percent wealth-constrained share or the 68 percent of
renters citing the down payment as a barrier, because those are object (i) and
(iv), not object (iii).

Before any of this is reported, confirm whether the model's 4.5 percent is
computed over the stock of households in the age and child cells or over the flow
of households that actually adjust housing that period. If it is a stock object
it is comparable to the Haurin range and is fine. If it is a flow object, it
should be compared to Barakova's mover-based 19 percentage points instead, and
4.5 percent would then be low by a factor of four.

**Proposed disciplining moment.** Use a level moment, not an elasticity. Build
from SCF microdata the share of renters aged 26 to 38 whose liquid and
semi-liquid assets, defined as in Barakova et al., fall below 20 percent of the
median owner-occupied unit price in their metropolitan area, and compute the same
object in the model's stationary distribution. This is object (i), it is
computable in both without any behavioral assumption, it maps to \(\phi = 0.80\)
exactly, and the Haurin and Barakova estimates of 58 to 65 percent give an
external anchor for whether the model's wealth distribution is right. Hold the
4.5 percent marginal share out as a validation object rather than a target, and
check it against 3.7 to 6.7 percent. If a second, behavioral check is wanted, the
Fuster and Zafar renter response is the cleanest one available: 58 percent of
renters raise willingness to pay when the required down payment falls from 20
percent to 5 percent, and renter willingness to pay rises about 40 percent. The
model can be asked the same question directly by re-solving at \(\phi = 0.95\).
That is a more informative check than the Acolin aggregate ownership counterfactual,
which mixes wealth, income and credit constraints that the model does not
separately contain.

## Not done

- Atalay, Li and Whelan (2021) could not be opened; the 7.5 percent figure is
  secondary-source only.
- Linneman and Wachter (1989), Zorn (1989), Duca and Rosenthal (1994) and
  Gyourko, Linneman and Wachter (1999) are paywalled. Their numbers here come
  from the Haurin, Hendershott and Wachter full text, which I did open.
- No age breakdown of the SCF credit-denial series exists in the published
  Bulletin; it would require the microdata.
- The Fannie Mae National Housing Survey site returned an access error on every
  attempt.
- The Joint Center for Housing Studies *America's Rental Housing* renter-savings
  tabulation was not retrieved.
- The SCF liquid-assets-versus-local-price share proposed above does not exist as
  a published statistic and would have to be constructed.

## References

Opened directly in this investigation unless noted.

- Acolin, Arthur, Jesse Bricker, Paul Calem and Susan Wachter. 2016. "Borrowing
  Constraints and Homeownership." *American Economic Review Papers and
  Proceedings* 106(5): 625-29. **Opened** (Wharton working paper 791 full text).
- Atalay, Kadir, Ang Li and Stephen Whelan. 2021. "Housing Wealth, Fertility
  Intentions and Fertility." *Journal of Housing Economics* 54. **Not verified.**
- Barakova, Irina, Raphael W. Bostic, Paul S. Calem and Susan M. Wachter. 2003.
  "Does Credit Quality Matter for Homeownership?" *Journal of Housing Economics*
  12(4): 318-36. **Opened**, including Tables 2 and 6.
- Bhutta, Neil and Daniel Ringo. 2021. "The Effect of Interest Rates on Home
  Buying: Evidence from a Shock to Mortgage Insurance Premiums." *Journal of
  Monetary Economics* 118: 195-211. **Opened** (FEDS 2017-086 working paper).
- Clark, Jeremy and Ana Ferrer. 2019. "The Effect of House Prices on Fertility:
  Evidence from Canada." *Economics: The Open-Access Journal* 13(1). **Opened**
  (Waterloo working paper).
- Daysal, N. Meltem, Michael F. Lovenheim, Nikolaj Siersbæk and David N. Wasser.
  2021. "Home Prices, Fertility, and Early-Life Health Outcomes." *Journal of
  Public Economics* 198. **Opened** (NBER working paper 27469).
- DeFusco, Anthony A., Stephanie Johnson and John Mondragon. 2020. "Regulating
  Household Leverage." *Review of Economic Studies* 87(2): 914-58. **Opened.**
- Dettling, Lisa J. and Melissa S. Kearney. 2014. "House Prices and Birth Rates:
  The Impact of the Real Estate Market on the Decision to Have a Baby." *Journal
  of Public Economics* 110: 82-100. **Opened** (NBER working paper 17485).
- Duca, John V. and Stuart S. Rosenthal. 1994. "Borrowing Constraints and Access
  to Owner-Occupied Housing." *Regional Science and Urban Economics* 24(3):
  301-22. **Not verified.**
- Engelhardt, Gary V. and Christopher J. Mayer. 1994. "Gifts for Home Purchase
  and Housing Market Behavior." *New England Economic Review* May/June: 47-58.
  **Opened.**
- Federal Reserve Board. 2023. "Changes in U.S. Family Finances from 2019 to
  2022: Evidence from the Survey of Consumer Finances." *Federal Reserve
  Bulletin*. **Opened**, Table 5.
- Fuster, Andreas and Basit Zafar. 2021. "The Sensitivity of Housing Demand to
  Financing Conditions: Evidence from a Survey." *American Economic Journal:
  Macroeconomics* 13(1): 231-65. **Opened** (New York Fed Staff Report 702).
- Greenwald, Daniel L. 2018. "The Mortgage Credit Channel of Macroeconomic
  Transmission." MIT Sloan working paper. **Opened**, though no numeric share at
  each constraint was located.
- Gyourko, Joseph, Peter Linneman and Susan Wachter. 1999. "Analyzing the
  Relationships among Race, Wealth, and Home Ownership in America." *Journal of
  Housing Economics* 8(2): 63-89. **Not verified.**
- Hacamo, Isaac. 2021. "The Babies of Mortgage Market Deregulation." *Review of
  Financial Studies* 34(2): 907-48. **Opened** (local copy,
  `/Users/tommasodesanto/Desktop/hhaa073.pdf`).
- Haurin, Donald R., Patric H. Hendershott and Susan M. Wachter. 1997.
  "Borrowing Constraints and the Tenure Choice of Young Households." *Journal of
  Housing Research* 8(2): 137-54. **Opened** (NBER working paper 5630, 1996).
- Kaplan, Greg, Kurt Mitman and Giovanni L. Violante. 2020. "The Housing Boom and
  Bust: Model Meets Evidence." *Journal of Political Economy* 128(9): 3285-345.
  **Opened** (NBER working paper 23694).
- Kearney, Melissa S. and Riley Wilson. 2018. "Male Earnings, Marriageable Men,
  and Nonmarital Fertility: Evidence from the Fracking Boom." *Review of
  Economics and Statistics* 100(4): 678-90. **Opened** (NBER working paper
  23408).
- Linneman, Peter and Susan Wachter. 1989. "The Impacts of Borrowing Constraints
  on Homeownership." *AREUEA Journal* 17(4): 389-402. **Not verified**; figures
  quoted from Haurin, Hendershott and Wachter.
- Lovenheim, Michael F. and Kevin J. Mumford. 2013. "Do Family Wealth Shocks
  Affect Fertility Choices? Evidence from the Housing Market." *Review of
  Economics and Statistics* 95(2): 464-75. **Opened.**
- National Association of Realtors. 2025. *Profile of Home Buyers and Sellers.*
  **Opened**, Exhibits 3-7, 5-3, 5-4.
- New York Fed. 2025. *SCE Housing Survey* chart packet. **Opened.**
- Urban Institute and Freddie Mac. 2018. *Barriers to Accessing Homeownership:
  Down Payment, Credit, and Affordability.* **Opened**; the 68 percent renter
  barrier figure originates in the 2018 Zillow Housing Aspirations Report.
- Zorn, Peter M. 1989. "Mobility-Tenure Decisions and Financial Credit: Do
  Mortgage Qualification Requirements Constrain Homeownership?" *AREUEA Journal*
  17(1): 1-16. **Not verified.**
