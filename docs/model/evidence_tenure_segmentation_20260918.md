# Is Family-Sized Housing Owner-Only? Evidence On The Tenure-Segmentation Assumption

Date: 2026-09-18. Status: evidence note. No model code was run and no model
or calibration change is proposed here.

> **Lead review (September 18, 13:30).** Two claims in this note are not verified by the lead and should be treated as the memo's, not the paper's, until checked. (1) The Kaplan, Mitman and Violante Table 6 renter shares by size class and the Section 5.3.1 "virtually indistinguishable" segmentation robustness are attributed to the published 2020 version; the NBER working paper (w23694) that I opened has the size grids {1.50, ..., 5.15} owner and {1.17, 1.50, 1.92} rental, as stated, but neither the table nor the robustness passage. (2) The statement that the cap is "the model's strongest fertility lever" (six to eight rooms: TFR +0.071, family-ownership gap 0.296 to 0.035) is taken from `output/model/fable_size_mapping_audit_20260701/HANDOFF_RESTART.md`, a July 1 audit on the earlier package at a coarse wealth grid whose mechanism numbers are recorded as unreliable. It has not been measured on the paper's model. The direct test was run on September 19: with the cap at eight rooms and nothing else changed, completed fertility is 1.876 against 1.872 at fixed child preference, childlessness and first-birth age are unchanged, and ownership falls from 0.46 to 0.31. On the paper's model the cap is an ownership lever, not a fertility lever, so the "strongest fertility lever" line below is withdrawn. The DUE finding (rental maximum equal to the owner maximum, 12.7, `DUE.txt` line 1211) is verified against the local text and does contradict the July cap note's attribution.

## (a) The question

The model gives renters and owners different size menus. A renter chooses
housing services \(h\) continuously on the interval \((0,\bar h^R]\) with
\(\bar h^R = 6\) rooms; an owner chooses from the discrete ladder
\(\mathcal H = \{2, 4, 6, 8, 9.5, 11\}\) rooms. Both objects are set in
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_housing_fertility/parameters.py`
(lines 122-125). The consequence is that the three largest rungs, 8, 9.5 and
11 rooms, are reachable only by buying. Because the model's child-space
requirement \(\bar h(n)\) rises with the number of children, a household that
wants a third or fourth child must cross the down-payment threshold
\((1-\phi)\) to get the space. That is the paper's central mechanism, so the
question is not cosmetic: is "large units are owner-only" a constraint on the
household's choice set, or is it a compact description of an allocation that
households would produce anyway through preferences, life-cycle wealth and
landlord economics? A cap that merely describes the observed allocation is not
a mechanism; a cap that binds is. The evidence below separates what is measured
from what is assumed.

## (b) What the paper's own figure shows

The slide is the frame "Big Units Are Owned --- Renting Is Not a Substitute" in
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/wedge2_hacamo_slides.tex`
(lines 78-100). Both panels come from the 2023 American Housing Survey national
public-use household file, version 1.1, weighted by `WEIGHT`, built by
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/ahs_supply_snapshot/analyze_ahs_family_unit_menu.R`.
The weighted read is 145.3 million housing units, 133.2 million occupied. Size
is the AHS bedroom count, binned 0-1, 2, 3, 4+. The left panel normalizes
within tenure; the right panel normalizes within size bin, with owner, renter
and vacant summing to 100 in each bin.

| Object (AHS 2023, national) | Value |
|---|---:|
| Share of owner-occupied units with 3+ bedrooms | 81.3% |
| Share of renter-occupied units with 3+ bedrooms | 30.2% |
| Owner share of 0-1 bedroom units | 12.8% |
| Owner share of 2 bedroom units | 39.4% |
| Owner share of 3 bedroom units | 75.4% |
| Owner share of 4+ bedroom units | 87.0% |
| Renter share of 4+ bedroom units | 9.2% |

Source file: `code/data/ahs_supply_snapshot/output_ahs_family_unit_menu_national/ahs_stock_by_bedroom_tenure.csv`.
The companion frame reports that 79% of households with children live in 3+
bedroom units against 59% of childless households, and that the parent share
within size bins rises 6, 21, 29, 41% across the four bins. I recomputed all of
these from the underlying CSVs and they reproduce exactly.

Two things the figure does not show. It is a cross-section of occupancy, not of
availability, so it cannot distinguish a choice-set restriction from sorting.
And it is stated in bedrooms, while the model is stated in rooms; the room
version of the same table is weaker, because 7+ room units are 89.0% owner,
8.0% renter and 3.0% vacant, and AHS does not publish a finer top bin.

## (c) What the data say

| Fact | Value | Source, year |
|---|---:|---|
| Renter share, studio | 86.5% | ACS 2024 1-year, B25042 |
| Renter share, 1 bedroom | 84.9% | ACS 2024 1-year, B25042 |
| Renter share, 2 bedrooms | 54.9% | ACS 2024 1-year, B25042 |
| Renter share, 3 bedrooms | 20.1% | ACS 2024 1-year, B25042 |
| Renter share, 4 bedrooms | 11.0% | ACS 2024 1-year, B25042 |
| Renter share, 5+ bedrooms | 8.3% | ACS 2024 1-year, B25042 |
| Renters living in 3+ bedroom units | 29.7% | ACS 2024 1-year, B25042 |
| Owners living in 3+ bedroom units | 80.5% | ACS 2024 1-year, B25042 |
| Renter share, 1-unit detached structures | 13.7% | ACS 2024 1-year, B25032 |
| Renter share, 5+ unit structures | 86.1% | ACS 2024 1-year, B25032 |
| 3+ bedroom share of units in 20+ unit buildings | 7.6% | AHS 2023 national |
| 3+ bedroom share of detached single-family units | 83.5% | AHS 2023 national |
| Renter-occupied 3+ bedroom units, count | 13.4 million | AHS 2023 national |
| Of those, share in single-family (detached + attached) | 68.7% | AHS 2023 national |
| Single-family rental stock | 14.9 million units, 2022; peak 16.1 million in 2016 | JCHS, *America's Rental Housing 2024* |
| Single-family built-as-rental completions | 26,000/yr in 2010 to 67,000/yr in 2022; 13% of new rental construction | JCHS 2024 |
| Renter share, 3 bedrooms, New York metro | 29.1% | ACS 2024 1-year, B25042 |
| Renter share, 3 bedrooms, Los Angeles metro | 28.2% | ACS 2024 1-year, B25042 |
| Renter share, 5+ bedrooms, New York metro | 8.9% | ACS 2024 1-year, B25042 |
| Prime-age childless renters with 6+ rooms | 13.8% | ACS/IPUMS 2012-2023, MMS metros |
| Prime-age childless renters with 7+ rooms | 6.0% | ACS/IPUMS 2012-2023, MMS metros |
| Owner minus renter mean rooms, prime-age childless | 2.42 rooms | ACS/IPUMS 2012-2023, MMS metros |
| Median gross rent by bedrooms (studio to 5+) | $1,321 / $1,301 / $1,490 / $1,677 / $2,069 / $2,182 | ACS 2024 1-year, B25031 |
| Implied rent per bedroom (1BR to 5+BR) | $1,301 / $745 / $559 / $517 / $436 | own arithmetic on B25031 |
| AHS mean rent per square foot by bedrooms | $2.09 / $1.49 / $1.18 / $1.14 | AHS 2023 national |
| Within-cell marginal rent premium above 5.5 rooms | $7.7 per room (se $15.6), statistically zero | ACS/IPUMS, 1,284,904 renter heads |
| Central 3+ bedroom rental price premium, metro-weighted | -0.02 log points | ACS/IPUMS 2012-2023 |

The ACS room and bedroom moments are in
`code/data/mms_center_periphery/output_intergen_one_market_targets/intergen_one_market_acs_housing_targets.csv`.
The within-cell rent regression is documented in
`output/model/fable_size_mapping_audit_20260701/SOFT_CAP_DESIGN_MEMO.md`, section 4:
rent on rooms with a kink, year by metro by location fixed effects, household
weights, standard errors clustered by metro. Below the knot the marginal room
is worth $85.1 (se $11.9); above it the extra premium is $7.7 (se $15.6) at a
5.5 room knot and negative at a 6 room knot. The renter room histogram declines
smoothly through the would-be kink (4 rooms 27.0%, 5 rooms 17.3%, 6 rooms
10.1%, 7 rooms 4.4%), so there is no bunching at the cap.

Three readings follow. First, large rental units exist in large absolute
numbers: 13.4 million occupied 3+ bedroom rentals nationally, which is not a
number compatible with a literal prohibition. Second, they are overwhelmingly
single-family houses rented out one at a time, not purpose-built family
apartments; only 7.6% of units in 20+ unit buildings have three or more
bedrooms. Third, the scarcity is a quantity phenomenon, not a price phenomenon:
rent per bedroom and rent per square foot fall monotonically with size, and the
within-market marginal price of family-sized space is statistically zero.

## (d) What the literature says

### Is the rental stock a different good?

The best direct evidence is Halket, Nesheim and Oswald (2020), which estimates
a selection model of how English dwellings are allocated across the
owner-occupied, private-rented and social sectors. Three of their results bear
on our assumption. The owner-occupied share of the Greater London stock rises
from 33.1% for dwellings under 50 square metres to 90.1% for dwellings over 100
square metres, with the private-rented share falling from 27.4% to 7.2%; the
owner share is 94.4% for detached houses and 20.7% for high-rise flats; and,
most usefully, unobserved quality in the rental sector declines with size, with
an average quality gap of roughly 22% between a 50 and a 100 square metre rental
property, a pattern absent on the owner side. That last finding is the strongest
available support for treating large rental and large owner units as different
goods rather than the same good under different tenure. It is estimated on
English data, so transferring it to the United States is an assumption.

Glaeser and Shapiro (2003) supply the American counterpart and locate the
segmentation on the structure-type axis rather than the size axis: 85.5% of
people in single-family detached houses own, 85.9% of people in multi-unit
buildings rent, and the explanation offered is the maintenance agency problem,
since the major maintenance decisions in a multi-unit building are building-
specific rather than apartment-specific and are therefore better made by a
single owner. This is the empirical form of the Henderson and Ioannides (1983)
rental externality: a tenant on a fixed-rent lease does not internalize the
effect of use on the structure, so landlords bear a moral-hazard cost, and the
cost is larger where the household's expected tenure is longer. Families are
precisely the long-tenure households. I was not able to open Henderson and
Ioannides directly, so that characterization is reported from secondary sources
and should be checked before it appears in the manuscript.

Two candidate sources turn out not to speak to the question. Sinai and Souleles
(2005) treat housing as a single good and argue for ownership as a hedge against
rent risk; there is no discussion of stock heterogeneity. Bachmann and Cooper
document gross flows between the renter and owner segments and say nothing about
the physical composition of either stock.

### How leading quantitative models handle the rental menu

| Paper | Rental size object | How disciplined |
|---|---|---|
| Kaplan, Mitman and Violante (2020) | partial segmentation on a discrete menu: owner sizes {1.50, 1.92, 2.46, 3.15, 4.03, 5.15}, rental sizes {1.17, 1.50, 1.92}; the smallest class cannot be owned and the largest four cannot be rented | inside the calibration, targeting an owner/renter average size ratio of 1.5 and an owner/renter earnings ratio of 2.1 |
| Greaney, Parkhomenko and Van Nieuwerburgh (2025), December 5 2025 draft | no size cap; owner sizes discrete on a six-point grid, rental size continuous, "maximum size is assumed to be the same as for owner-occupied houses: 12.7" | owner grid from ACS values and AHS size percentiles; no rental-bound row in the calibration table |
| Sommer and Sullivan (2018) | none; one linear technology converts a unit of housing stock into a unit of shelter services, so the two stocks are the same good and fully convertible, with household landlords | per-unit landlord management cost only |
| Landvoigt, Piazzesi and Schneider (2015) | no rental sector at all; continuous owner quality assigned in equilibrium | renting assumed away by appeal to the rarity of owner-to-renter transitions |
| Chambers, Garriga and Schlagenhauf (2009) | owner-side minimum size, the mirror object | estimated as part of the exactly identified problem |
| Halket and Pignatti (2015) | none; scarcity is an equilibrium search and screening outcome | Craigslist listings by zip and bedrooms, time on market, rent-to-price |
| Halket and Vasudev (2014) | none; the housing choice space is convex and does not depend on tenure. Renters end up in smaller units endogenously: average rental size 0.70 against average owned size 1.11 | tenure sorting from down payments, transaction costs and tax wedges |
| Couillard (2025) | none; household type is the full cross of age, tenure, family status and size | rent coefficients estimated with additive heterogeneity by dimension |

Kaplan, Mitman and Violante (2020) is the closest precedent and the one to
cite, but it comes with two warnings. Their Table 6 reports that 10% of renters
in the AHS occupy one of the top four size classes while the model puts 0% there
by construction, which is the same failure our cap produces. And their Section
5.3.1 re-runs the model with no segmentation, half segmentation and full
segmentation and reports results "virtually indistinguishable from baseline",
because in their setting renters are not constrained in housing services. Our
model is the opposite case: the cap is a first-order lever. A referee who knows
this literature will notice that we are leaning on a device its own originators
showed to be inessential in their application, without having shown it is
essential in ours for a reason other than assumption.

The reference urban model does not support the assumption. A July 2026 internal
memo recorded that it sets the rental cap at 1.32 times the owner grid minimum,
and the canonical \(\bar h^R = 6\) was defended partly as that rule applied to
our own sample
(`output/model/fable_size_mapping_audit_20260701/CAP_EXTERNAL_NOTE_20260706.md`).
The December 5, 2025 draft in the repository
(`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/DUE.txt`,
section 2.1.5 and the calibration table) says the opposite in plain text: "The
lower bound of the size of a rental unit is zero and the maximum size is
assumed to be the same as for owner-occupied houses: 12.7." The calibration
table lists the smallest and largest owner size and has no rental-bound row.
The model's rental and owner size sets do differ, but only in that one is
discrete and the other continuous, which is a computational convenience the
authors state as such. Either the earlier draft differed or the earlier memo
misread it. Until that is resolved, the current draft cannot be cited as a
precedent for tenure segmentation by size, and the external anchoring rule for
\(\bar h^R\) loses its stated pedigree.

Halket and Pignatti (2015) offer the microfoundation the hard cap is a reduced
form of. Houses are ex ante identical; rental scarcity arises because landlords
post lower rents but screen harder against long-duration households, so families
face low matching probabilities rather than high posted rents. Thin quantities
with flat prices is exactly the signature in our data. They describe the hard
cap as the shortcut they replace, naming Chambers, Garriga and Schlagenhauf
(2009), Fisher and Gervais (2011) and Amior and Halket (2014) as models that
impose large houses being supplied only on the owner-occupied market.

Finally, the conversion evidence cuts against a literal cap. Kaplan, Mitman and
Violante cite roughly 3 million single-family units converted from
owner-occupied to rental use between 2007 and 2011 against 4.5 million
owner-to-renter household transitions, and conclude that conversions account for
the bulk of those transitions. Lambie-Hanson, Li and Slonkosky document the
institutional side: institutional investors averaged 11.7% of single-family
purchases in their sample and all investor purchases reached 38.1% of
transactions, with the rise in institutional share explaining 58% of the
increase in real house price growth in affected counties. Gurun, Wu, Xiao and
Yang on the welfare effects of institutional landlords could not be opened and
is marked not verified.

## (e) Identification: what is established, what is conjectured

Established, in the sense that it is measured and reproducible:

1. The renter share falls steeply and monotonically in unit size, to roughly
   9 to 11% at four or more bedrooms, in both AHS 2023 and ACS 2024, nationally
   and in New York and Los Angeles.
2. Within year by metro by location cells, the marginal rent of a family-sized
   room is statistically zero, and there is no bunching of renter unit sizes at
   any candidate cap.
3. Rent per bedroom and rent per square foot fall with unit size. Large rentals
   are cheap per unit of space, not expensive.
4. Family-sized rentals are concentrated in single-family structures; large
   multifamily buildings essentially do not contain them.
5. In England, where the question has been studied with a selection model,
   large rental units are systematically lower in unobserved quality than small
   ones, by about 22% between 50 and 100 square metres, while owner-occupied
   quality does not vary with size (Halket, Nesheim and Oswald 2020).

Conjectured, in the sense that the paper cannot currently demonstrate it:

1. That the observed size distribution of renters reflects a choice-set
   restriction rather than sorting. Renters are younger, poorer and more mobile
   than owners; all three predict small units without any constraint.
2. That \(\bar h^R = 6\) is the right level. The current defense is that six
   rooms is the ninetieth percentile of renter unit sizes in the matched ACS
   sample. That is a description of the outcome, not identification: fitting a
   constraint to the upper tail of the distribution it is supposed to generate
   is circular, and the model then predicts zero renters above six rooms against
   6.0% of prime-age childless renters observed at seven or more.

What would identify a constraint rather than a preference outcome, and what
exists on each:

- **Rent per unit of size by size.** A binding quantity constraint with excess
  demand should show up as a premium on large rentals. It does not: the
  within-cell premium is zero and the gradient runs the other way. This is
  evidence *against* a price-rationed cap and *for* a non-price rationing story,
  but it is also consistent with no constraint at all.
- **Unobserved quality by size and tenure.** If large rentals are a genuinely
  worse good, a household wanting family space really does have to buy to get
  it, and the cap is a defensible reduced form. Halket, Nesheim and Oswald
  (2020) find exactly this, but on English data with an English social-housing
  sector. No American replication has been located. This is the single most
  valuable piece of evidence the paper could add.
- **Vacancy and time on market by size.** Under quantity rationing, large
  rentals should clear faster and sit vacant less. Halket and Pignatti measured
  time on market by zip and bedrooms with Craigslist listings and found
  screening against long-duration households. The equivalent tabulation for the
  current sample does not exist: I could not find a published Census Housing
  Vacancy Survey table by bedroom count, and the AHS Table Creator cross-tab has
  not been run here.
- **Conversion frictions.** The cleanest test is whether the stock moves between
  tenures at the family-sized margin. It does, at scale: roughly 3 million
  single-family units were converted from owner-occupied to rental use between
  2007 and 2011, and institutional investors reached 11.7% of single-family
  purchases in the Lambie-Hanson, Li and Slonkosky sample. If large units were
  technologically owner-only, that conversion could not have happened. This is
  the strongest single argument against a literal cap and the best case for
  reading \(\bar h^R\) as a friction that binds in the short run rather than a
  technological restriction.
- **Condominium conversion law and multifamily zoning.** These would give
  policy variation in the supply of large non-owner units. No evidence has been
  assembled here.

What the paper can claim: that the family-sized stock is, as a matter of
measured fact, overwhelmingly owner-occupied; that this is a quantity fact and
not a price fact; and that a hard cap is a standard reduced-form device for
representing it, with Kaplan, Mitman and Violante (2020) as the closest
precedent and Halket and Pignatti (2015) as the microfoundation. What the paper
cannot claim: that the cap is identified, that it is a technological
restriction, or that the reference urban model imposes one.

## (f) Bottom line

The strongest defensible statement is that family-sized housing in the United
States is overwhelmingly owner-occupied, that the scarcity of large rentals
shows up in quantities and in unobserved quality rather than in posted prices,
and that the model's hard cap at six rooms is a reduced form for a screening or
availability friction of the Halket-Pignatti type, with the same status as the
partial segmentation in Kaplan, Mitman and Violante (2020). The cap level itself
is set to the ninetieth percentile of observed renter unit sizes, which is a
calibration convention rather than an identified parameter, and the model's
implication of zero renters above six rooms is false against a data share of
6.0% among prime-age childless renters. The weakest point, and the one a referee
will attack first, is that the cap is the model's strongest fertility lever,
since moving it from six to eight rooms raises the baseline total fertility rate
by 0.071 and collapses the family-ownership gap from 0.296 to 0.035, so a
first-order quantitative result rests on an unidentified assumption that the
originators of the device themselves showed to be inessential in their own
application. The second line of attack is conversion: roughly 3 million
single-family units moved from owner-occupied to rental use between 2007 and
2011, which rules out any technological reading of the cap. The honest framing
is to present the cap as a calibrated approximation to an availability and
quality friction, report the fit over a band such as
\(\bar h^R \in \{5.5, 6.0, 6.5\}\), state the counterfactual implication the cap
gets wrong, and let the mechanism rest on the down-payment threshold rather than
on segmentation alone.

## Not done

- Henderson and Ioannides (1983) could not be opened; JSTOR blocked access and
  no working-paper mirror was found. The rental-externality characterization
  above is from secondary sources, principally Glaeser and Shapiro (2003).
- Gurun, Wu, Xiao and Yang on institutional landlords and renter welfare could
  not be opened past the paywall. Its reported magnitudes are not verified.
- Glaeser and Gyourko's *Rethinking Federal Housing Policy* was not opened;
  Glaeser and Shapiro (2003) was opened instead and carries the stylized fact.
- No vacancy rate or time-on-market tabulation by bedroom count was obtained.
  Census Housing Vacancy Survey historical tables break out vacancy by
  units-in-structure and region, not bedrooms; AHS Table Creator was not run.
- The ACS 2024 figures were pulled through the Census Reporter API rather than
  data.census.gov or the Census API directly, and are 1-year rather than 5-year
  estimates. They should be re-pulled from the primary source before print.
- The internal discrepancy on the reference urban model's rental cap is
  reported, not resolved. Resolving it requires locating the earlier draft that
  the July 2026 memo read.
- No American replication of the Halket-Nesheim-Oswald quality-by-size-and-
  tenure result was located. This is the highest-value missing evidence.

## References

Each entry states whether the source was opened at first hand for this note,
opened by a delegated agent that reported the extracted text, or verified
earlier in the repository's own July 2026 memo.

- Amior, M. and Halket, J. (2014), "Do households use home-ownership to insure
  themselves? Evidence across US cities," *Quantitative Economics* 5(3).
  Not opened; named only as a citation inside Halket and Pignatti.
- Bachmann, R. and Cooper, D., "The Ins and Arounds in the U.S. Housing Market."
  Opened in full by the delegated agent; contains turnover flows only, nothing
  on stock composition.
- Chambers, M., Garriga, C. and Schlagenhauf, D. (2009), "Accounting for Changes
  in the Homeownership Rate," *International Economic Review* 50(3). Verified in
  the internal July 2026 memo against the St. Louis Fed working paper; not
  opened for this note.
- Couillard, B. (2025), "Build, Baby, Build: How Housing Shapes Fertility," job
  market paper, University of Toronto. Opened at first hand; local text extract
  at `docs/literature/couillard/BuildBabyBuild_Couillard_2025.txt`.
- Coven, J., Golder, S., Gupta, A. and Ndiaye, A. (2024), "Property Taxes and
  Housing Allocation Under Financial Constraints," CESifo WP 11203. Verified in
  the internal memo; not opened for this note.
- Glaeser, E. and Shapiro, J. (2003), "The Benefits of the Home Mortgage
  Interest Deduction," NBER WP 9284. Opened in full by the delegated agent.
- Greaney, B., Parkhomenko, A. and Van Nieuwerburgh, S. (2025), "Dynamic Urban
  Economics," draft of December 5, 2025. Opened at first hand; full local text
  at `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/DUE.txt`.
- Gurun, U., Wu, J., Xiao, S. C. and Yang, S. (2023), "Do Wall Street Landlords
  Undermine Renters' Welfare?", *Review of Financial Studies* 36(1). NOT
  VERIFIED; paywalled, abstract and search snippets only.
- Halket, J. and Pignatti Morano di Custoza, M. (2015), "Homeownership and the
  scarcity of rentals," *Journal of Monetary Economics* 76. Verified in the
  internal memo against the published PDF; not opened for this note.
- Halket, J. and Vasudev, S. (2014), "Saving up or settling down: Home ownership
  over the life cycle," *Review of Economic Dynamics* 17(2). Opened in full by
  the delegated agent at the UCL Discovery copy.
- Halket, J., Nesheim, L. and Oswald, F. (2020), "The Housing Stock, Housing
  Prices, and User Costs: The Roles of Location, Structure, and Unobserved
  Quality," *International Economic Review* 61(4). Opened in full by the
  delegated agent at the CEMMAP working paper CWP7315; the published Wiley
  version was paywalled.
- Henderson, J. V. and Ioannides, Y. M. (1983), "A Model of Housing Tenure
  Choice," *American Economic Review* 73(1), 98-113. NOT VERIFIED; JSTOR access
  blocked and no mirror found.
- Joint Center for Housing Studies (2024), *America's Rental Housing 2024*.
  Opened as PDF by the delegated agent.
- Kaplan, G., Mitman, K. and Violante, G. (2020), "The Housing Boom and Bust:
  Model Meets Evidence," *Journal of Political Economy* 128(9). Opened in full
  by the delegated agent at the authors' revision PDF; consistent with the
  internal memo's earlier check against NBER WP 23694.
- Lambie-Hanson, L., Li, W. and Slonkosky, M. (2019), "Institutional Investors
  and the U.S. Housing Recovery," Federal Reserve Bank of Philadelphia WP 19-45.
  Opened in full by the delegated agent.
- Landvoigt, T., Piazzesi, M. and Schneider, M. (2015), "The Housing Market(s)
  of San Diego," *American Economic Review* 105(4). Opened in full by the
  delegated agent at NBER WP 17723.
- Sinai, T. and Souleles, N. (2005), "Owner-Occupied Housing as a Hedge Against
  Rent Risk," *Quarterly Journal of Economics* 120(2). Opened in full by the
  delegated agent at the Philadelphia Fed working paper; contains nothing on
  stock heterogeneity.
- Sommer, K. and Sullivan, P. (2018), "Implications of US Tax Policy for House
  Prices, Rents, and Homeownership," *American Economic Review* 108(2), 241-274.
  Opened in full by the delegated agent at the 2013 working-paper version;
  consistent with the internal memo's check of the published section I.B.
- US Census Bureau, American Community Survey 2024 1-year estimates, tables
  B25031, B25032, B25042. Retrieved through the Census Reporter API, which is a
  pass-through for the Census B-table cells, not the Census API itself.
- US Census Bureau, American Housing Survey 2023 national public use file,
  version 1.1. Used directly through the repository build.
