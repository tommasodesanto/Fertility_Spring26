# Intergen Size-Tenure Supply And Pricing Audit

Date: 2026-06-23

Status: empirical/literature discussion memo. This does not change model logic,
solver logic, calibration targets, or production code.

## Question

The current one-market intergen model lets renters choose continuous housing
services up to \(h^R_{\max}=8\), while owners choose discrete rungs. In the
June 23 diagnostic point, renters are not mainly stuck at the rental cap, yet
the model puts too many childless renters into mid-sized units. The question is
whether this is a grid problem, a measurement problem, or a missing supply/menu
object.

The short answer is: the evidence points more to a missing size-by-tenure
product menu than to a simple cap or per-room price fix.

## Local Empirical Inputs

Lightweight summarizer:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
code/model/.venv/bin/python code/data/mms_center_periphery/summarize_intergen_size_tenure_supply_pricing.py
```

Derived output:

`code/data/mms_center_periphery/output_intergen_size_tenure_supply_pricing/`

The script reads existing precomputed ACS/MMS packet CSVs. It does not read the
raw 9.2GB IPUMS extract.

Source packets:

- `code/data/mms_center_periphery/output_family_size_supply/acs_family_size_supply_cells.csv`
- `code/data/mms_center_periphery/output_couillard_bedroom_supply/acs_bedroom_supply_cells.csv`
- `code/data/mms_center_periphery/output_intergen_quality_adjusted_room_targets/intergen_quality_adjusted_room_structure_summary.csv`

These are descriptive ACS/IPUMS household-head facts for 2012--2023, matched to
the MMS geography, with middle PUMAs collapsed to center in the underlying
packets.

## Supply Menu Facts

All household heads, ACS/MMS room-bin packet:

| Tenure | Size bin | Share within tenure | Mean rooms | Mean rent | Mean rent per room |
|---|---:|---:|---:|---:|---:|
| Owner | 1--4 rooms | 0.129 | 3.603 |  |  |
| Owner | 5--6 rooms | 0.387 | 5.547 |  |  |
| Owner | 7+ rooms | 0.483 | 8.620 |  |  |
| Renter | 1--4 rooms | 0.630 | 3.108 | 1198 | 454 |
| Renter | 5--6 rooms | 0.278 | 5.374 | 1322 | 247 |
| Renter | 7+ rooms | 0.092 | 8.062 | 1623 | 206 |

All household heads, ACS/MMS bedroom packet:

| Tenure | Bedroom bin | Share within tenure | Mean bedrooms | Mean rooms | Median rent per bedroom |
|---|---:|---:|---:|---:|---:|
| Owner | 0--1 bedrooms | 0.005 | 1.000 | 1.357 |  |
| Owner | 2 bedrooms | 0.026 | 2.000 | 3.351 |  |
| Owner | 3+ bedrooms | 0.969 | 4.293 | 6.901 |  |
| Renter | 0--1 bedrooms | 0.067 | 1.000 | 1.124 | 1020 |
| Renter | 2 bedrooms | 0.293 | 2.000 | 3.000 | 504 |
| Renter | 3+ bedrooms | 0.640 | 3.550 | 5.059 | 364 |

Prime-age childless heads, ACS/MMS structure summary:

| Structure | Owner share within tenure | Renter share within tenure | Owner mean rooms | Renter mean rooms |
|---|---:|---:|---:|---:|
| 1-family detached | 0.735 | 0.179 | 6.719 | 5.505 |
| 1-family attached | 0.117 | 0.060 | 5.625 | 4.734 |
| Multifamily, 2+ units | 0.120 | 0.744 | 4.099 | 3.310 |
| Mobile home/trailer | 0.027 | 0.015 | 4.965 | 4.354 |

The empirical menu is strongly tenure-segmented. Owners are much more likely to
occupy large/detached units; renters are much more likely to occupy small and
multifamily units.

## Current Diagnostic Model Contrast

Current diagnostic point:

`output/model/intergen_room_distribution_current_best_20260623/`

Grid:

\[
H^{own}=\{2,4,6,8,10\},\qquad h^R_{\max}=8.
\]

Prime-age childless room-bin shares:

| Tenure | Bin | Data | Model | Model - data |
|---|---:|---:|---:|---:|
| Renter | \(\le 4\) rooms | 0.719 | 0.415 | -0.304 |
| Renter | 4--6 rooms | 0.221 | 0.465 | 0.244 |
| Renter | \(>6\) rooms | 0.060 | 0.120 | 0.059 |
| Owner | \(\le 4\) rooms | 0.201 | 0.234 | 0.034 |
| Owner | 4--6 rooms | 0.419 | 0.499 | 0.080 |
| Owner | \(>6\) rooms | 0.380 | 0.267 | -0.114 |

Other diagnostic facts:

- Renter cap share among prime-age childless renters: 0.005.
- Owner rung shares: 2 rooms 0.000, 4 rooms 0.234, 6 rooms 0.499, 8 rooms
  0.267, 10 rooms 0.000.
- Mean rooms: renter 4.561, owner 6.065, owner-renter gap 1.504.

This matters for the grid discussion. The current renter cap is barely binding,
so lowering \(h^R_{\max}\) is a blunt way to move the model. The model's
problem is mostly interior renter demand: too little mass in small rentals and
too much in mid-sized rentals. The owner grid also already includes a 2-room
rung, but that rung is unused, so extending ownership to zero is unlikely to
solve the size distribution by itself.

## Pricing Interpretation

The rent evidence is not "large rentals have a high per-room price." In the
ACS/MMS room packet, total monthly rent rises from about 1198 for 1--4 room
rentals to 1623 for 7+ room rentals, but mean rent per room falls from about
454 to 206. The bedroom packet has the same pattern: median rent per bedroom is
1020 for 0--1 bedroom rentals, 504 for 2-bedroom rentals, and 364 for 3+
bedroom rentals.

So a convex rental price schedule by physical rooms is not the clean empirical
fix. A larger rental is more expensive in total dollars, but not more expensive
per room in the pooled ACS/MMS descriptives. This likely reflects location,
structure, vintage, and selection: large rentals are often different products,
not simply more of the same homogeneous housing service.

The model implication is that one aggregate housing-services price \(q\) is too
coarse for this margin. A better object is a size-by-tenure product menu or
supply composition:

\[
\mathcal H^R = \{small\ rental,\ large\ rental\}, \qquad
\mathcal H^O = \{starter\ owner,\ large/quality\ owner\},
\]

with either distinct prices, availability/capacity, or taste/quality terms. The
key point is that this should be introduced as an economically interpretable
market/menu object, not as arbitrary grid surgery.

## Literature Check

The literature supports taking tenure-product segmentation seriously.

- Census ACS documentation for 2023 housing availability emphasizes that
  bedrooms should be read jointly with structure type and points users to
  tables for rooms, bedrooms, structure type, and occupancy status. It reports
  that single-unit structures with three or more bedrooms are 54.3 percent of
  all housing units in 2023, while common multiunit products are much smaller:
  [Census ACS-61, 2025](https://www2.census.gov/library/publications/2025/demo/acs-61.pdf).
- Coulson and Fisher's "Structure and Tenure" states the basic empirical fact
  directly: single-family units overwhelmingly owner-occupied, multifamily
  units overwhelmingly rented; their evidence highlights unit size and quality
  as strong predictors of tenure and building ownership:
  [Coulson and Fisher, 2012](https://www.smeal.psu.edu/bires/documents/Structure%20and%20Tenure3-21-12.pdf).
- Henderson and Ioannides is the classic tenure-choice reference; a useful HUD
  review summarizes the mechanism as a principal-agent/maintenance problem:
  owner-occupants internalize maintenance incentives more directly than
  renters and landlords:
  [Green, HUD User, 2000](https://www.huduser.gov/portal/publications/pdf/brd/05Green.pdf).
- Glaeser argues that rental housing in the United States is fundamentally
  linked to multifamily structures, so policies affecting tenure also affect
  density and structure composition:
  [Glaeser, 2011](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1914468).
- Sommer, Sullivan, and Verbrugge model house prices and rents with explicit
  owned and rental property markets. That is closer to the object we may need
  than a single homogeneous housing-services market:
  [Sommer, Sullivan, and Verbrugge, 2013](https://ideas.repec.org/a/eee/moneco/v60y2013i7p854-870.html).
- Sinai and Souleles show that homeownership provides a hedge against future
  rent risk. This is not a size-menu paper, but it reinforces that ownership is
  not only a financing technology for the same flow service; tenure changes the
  intertemporal risk object:
  [Sinai and Souleles, 2005](https://academic.oup.com/qje/article-abstract/120/2/763/1933972).
- The 2026 Joint Center rental report flags modern rental-market context:
  rental cost burdens remain high, multifamily construction is cooling, and the
  aging rental stock needs investment. This supports treating rental supply as
  its own constrained object:
  [JCHS, America's Rental Housing 2026](https://www.jchs.harvard.edu/americas-rental-housing-2026).

## What This Means For The Model

I would not make the next move "lower \(h^R_{\max}\)" as a production fix. It
may improve the scalar room gap, but it treats a product-composition fact as a
pure upper-bound fact. Since the current cap is barely binding, this would move
policy functions through an artificial constraint rather than through observed
supply/menu economics.

I also would not use the structure-adjusted room gap as the hard target by
itself. The previous quality-adjustment audit showed that conditioning on
structure collapses the owner-renter gap, but structure is exactly the product
margin the current model lacks. Controlling it away without a replacement
moment would erase identification of the tenure-product block.

The cleaner interpretation is:

1. Keep the raw or location-adjusted owner-renter size gap as a diagnostic
   pressure point.
2. Add explicit empirical targets for product composition if the model gains a
   product/menu state, e.g. renter share in small units, owner share in large
   units, or detached/single-family shares by tenure.
3. Consider a two-product rental/owner menu before another long calibration:
   small rental, family rental, starter owner, large owner. This preserves the
   idea that agents choose housing, but lets them choose from empirically
   different markets.
4. If staying with one aggregate housing-services market for now, do not expect
   the search to fully resolve renter size distribution, owner-renter room gap,
   young ownership, and old retention jointly. The model is missing a dimension
   that the data say is first-order.

## Identification-Preserving Target Implication

A possible future target block, conditional on adding product/menu states, is:

| Parameter block | Current pressure | Candidate replacement/addition |
|---|---|---|
| Housing need/Stone-Geary block \((\bar h_0,\bar h_n,\bar h_{\mathrm{jump}})\) | mean rooms and birth housing responses | renter small-unit share, owner large-unit share, parent-childless room gap |
| Tenure utility/entry block \((\chi,b_{\mathrm{entry}},\phi,\psi^{PTI})\) | ownership rates and family ownership gaps | ownership by age plus renter/owner product shares |
| Supply/menu block | aggregate \(H^D=H^S\) only | distinct rental/owner product quantities or fixed empirical menu shares |

This would not reduce identification. It would shift the room-margin discipline
from a single mean room gap toward distributional/product moments that match
the observed market.

## Bottom Line

Your concern is right in the narrow sense that "agents should choose" housing
sizes, not be forced into arbitrary bins. But the data suggest the current
choice set is too abstract: it lets renters buy interior amounts of homogeneous
housing services that, empirically, are tied to a different stock of structures,
locations, and quality. The most defensible model fix is not arbitrary grid
manipulation; it is a small size-by-tenure supply/menu extension.
