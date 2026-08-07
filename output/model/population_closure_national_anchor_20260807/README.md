# National outside-origin entry anchor

## Result and recommendation

The headline national anchor is **3.714%**: among ACS household heads aged 18--22, the `HHWT`-weighted share who are foreign-born and whose reported year of immigration is within four calendar years of the survey year. This is the closest available national analogue to an outside entrant in the model; it is far below the current quota-closure floors (14.26% and 14.32%). Thus neither arm is feasible under this anchor. The main caveat is conceptual: immigration arrival is not necessarily the same date as forming a separate household, and the ACS year-only `YRIMMIG` measure (which is partly imputed) cannot identify that transition. It nevertheless measures national foreign inflow rather than domestic relocation, which nets to zero nationally.

## Data and implementation

Source: [`extract27.dta`](../../../code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta), the on-disk ACS/IPUMS extract used by the MMS workflow. The annual ACS records cover 2012--2023. The extract carries `YEAR`, `SAMPLE`, `STATEFIP`, `GQ`, `PERNUM`, `RELATE`, `AGE`, `HHWT`, `PERWT`, `BPL`, `BPLD`, `CITIZEN`, and `YRIMMIG`. This calculation uses `CITIZEN` rather than birthplace: IPUMS codes 2--5 are foreign-born, while code 1 (born abroad of American parents) is excluded.

`stage_acs_blocks.py` reads non-overlapping fixed record blocks directly from the 9.2-GB Stata file and writes only the small selected samples to `staged_blocks/`; `summarize_staged_entry.py` builds the final CSV. The block staging is purely a resource-safe way to stream the existing extract. The blocks span source records 20,000,000--59,046,775; all retained annual records are exactly 2012--2023. Within household units (`GQ` 1, 2, or 5), `PERNUM == 1` and `RELATE == 1` agree exactly: 14,829,596 observations satisfy both and none satisfies only one.

## Exact definitions

The national geographic sample is the 50 states plus DC (`STATEFIP` 1--56), excluding Puerto Rico and group quarters for the head measures. A household head has `GQ` in {1,2,5}, `PERNUM == 1`, `RELATE == 1`, and positive `HHWT`. Head measures are weighted by `HHWT`. The person analogue has age 18--22 and positive `PERWT`, with no head or group-quarters restriction, and is weighted by `PERWT`.

For a respondent observed in year \(t\), “within \(k\) years” means \(0 \leq t-\mathrm{YRIMMIG} \leq k\), where \(k\in\{4,8\}\). The stock measure is simply foreign-born, regardless of `YRIMMIG`. The final CSV, [`national_outside_entry_measures.csv`](national_outside_entry_measures.csv), contains each measure both pooled and by year, including unweighted sample size and weighted numerator and denominator.

## Pooled results and quota feasibility

The renewal multiplier is \((1-s)/s\). The quota requirement is \(s_{out}\geq 1-B_0/E_0\), equal to 14.26% for the floor arm and 14.32% for the tilt arm. Every candidate below therefore leaves **both arms infeasible**.

| Measure | Value | Unweighted count | Floor / tilt feasible? | \((1-s)/s\) |
|---|---:|---:|---|---:|
| Heads 18--22, foreign-born, arrival \(\leq4\) years (headline) | 3.714% | 214,919 | No / No | 25.92 |
| Heads 18--22, foreign-born, arrival \(\leq8\) years | 5.361% | 214,919 | No / No | 17.65 |
| Heads 18--25, foreign-born, arrival \(\leq4\) years | 3.465% | 574,005 | No / No | 27.86 |
| Heads 18--25, foreign-born, arrival \(\leq8\) years | 5.558% | 574,005 | No / No | 16.99 |
| Heads 18--22, foreign-born stock (upper bound) | 9.036% | 214,919 | No / No | 10.07 |
| Persons 18--22, foreign-born, arrival \(\leq4\) years | 2.885% | 2,455,577 | No / No | 33.66 |

No annual head cell has fewer than 5,000 unweighted heads; the minimum is 12,805 (2020, ages 18--22). All reported shares are strictly between zero and one.

## Verification

`summarize_staged_entry.py` was run twice on the same staged inputs; the two `national_outside_entry_measures.csv` files were byte-identical (`cmp -s`).

## Lead interpretation (appended after review, 2026-08-07 night)

The measurement is accepted as correct for its stated object, and the object
itself is the finding. Every entry-age-window definition (2.9--9.9 percent)
sits far below the 14.3 percent feasibility floor, yet the United States is
approximately stationary under below-replacement fertility. The
reconciliation: real immigration sustains the population mostly at ages
above the model's stylized entry window (arrivals in the mid-20s through
40s), so the share of age-18--22 household entrants who are recent
immigrants is a different -- and much smaller -- object than the total
immigration contribution to cohort replacement. A back-of-envelope
age-aggregated anchor (net international migration per period relative to
cohort-equivalent inflow, published Census magnitudes, TO BE VERIFIED) is on
the order of 0.17--0.22, which brackets the current across-CBSA 0.169 and
keeps both certified arms feasible.

Decision this puts to the author (with the advisor):

1. Anchor the fixed inflow to the age-aggregated object (net international
   migration in entrant-household equivalents). Keeps the stylized age-18
   entry as an explicit compression of arrival ages; the 0.169 provisional
   value is approximately right and the quota results stand as computed.
2. Take the entry-age flow literally (3.7 percent). Then the closure as
   specified is infeasible at every certified arm, and honesty requires
   either multi-age entry (an age-structured immigrant inflow vector) or a
   fertility block that clears a much higher native-replacement bar.

Option 1 is the recommended reading for the current paper; option 2 is the
model extension under which the entry-age measurement becomes the right
anchor. Do not average the two concepts.
