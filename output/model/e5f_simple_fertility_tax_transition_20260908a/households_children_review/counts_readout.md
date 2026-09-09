# Households with children: count and population audit

Here a child household means at least one child in the model's dependent/at-home proxy state. It is not a literal count of households with children under 18: child departure is geometric, with probability 4/18 each four-year step.

A birth changes an existing household's family state. It does not create an additional household decision unit immediately. The total number of households and the number with dependent children are different stocks.

All counts below are per 100 households in the common 2023 economy. They are normalized model counts, not US household totals. Each row compares the rebated 1% baseline with the rebated 2% reform at the same date.

## All household ages

| Year | Total HH, baseline | Child HH, baseline | Child HH, reform | Child-HH reform % | Child-HH share, baseline | Child-HH share, reform |
|---|---:|---:|---:|---:|---:|---:|
| 2023 | 100.0000 | 36.4149 | 36.4488 | 0.0929 | 36.4149 | 36.4488 |
| 2043 | 93.5753 | 25.6914 | 25.8189 | 0.4960 | 27.4554 | 27.5853 |
| 2063 | 74.0860 | 19.5240 | 19.6991 | 0.8972 | 26.3531 | 26.5366 |

## Young households: model nodes 26, 30, 34

| Year | Total HH, baseline | Child HH, baseline | Child HH, reform | Child-HH reform % | Child-HH share, baseline | Child-HH share, reform |
|---|---:|---:|---:|---:|---:|---:|
| 2023 | 19.8786 | 9.4740 | 9.4928 | 0.1989 | 47.6591 | 47.7538 |
| 2043 | 15.2985 | 6.4727 | 6.5337 | 0.9421 | 42.3095 | 42.7081 |
| 2063 | 10.6908 | 4.5853 | 4.6515 | 1.4451 | 42.8897 | 43.2084 |

## The household-state partition

Zero-parity households have no previous modeled birth. Ever-parent households include both households with dependent children and positive-parity households with no child currently at home. The latter are labeled empty-nest states here. Accordingly H = zero parity + dependent-child households + empty-nest households, and ever-parent households = dependent-child households + empty-nest households. All identities and lifecycle totals pass.

| Year | Zero parity, baseline | Ever parent, baseline | Empty nest, baseline | Ever parent, reform | Empty nest, reform |
|---|---:|---:|---:|---:|---:|
| 2023 | 24.6013 | 75.3987 | 38.9838 | 75.4274 | 38.9786 |
| 2043 | 28.7084 | 64.8669 | 39.1755 | 65.0219 | 39.2030 |
| 2063 | 25.9736 | 48.1124 | 28.5884 | 48.3818 | 28.6827 |

## Interpretation

In the rebated 1% baseline, dependent-child household counts fall 46.38% from 2023 to 2063; total household counts fall 25.91%. At 2063 the reform nevertheless has 0.897% more dependent-child households than its matched baseline. A positive reform effect therefore does not imply growth relative to 2023.

On impact, total household mass is unchanged, but dependent-child household mass increases 0.09292%. The change in ever-parent mass is 0.0003235978; the change in empty-nest mass is -0.0000585111. Their difference exactly gives the change in dependent-child household mass.

A first birth maps (n,m)=(0,0) to (1,1): the same household moves from zero parity into both ever-parent and dependent-child categories. A later birth maps (n,m) to (n+1,m+1). It creates no new household and no new ever-parent household. If m was already positive, it adds no dependent-child household either; if m=0 at positive parity, it moves an existing empty-nest household back into the dependent-child category. Child maturation can move households in the opposite direction.

Total household mass changes through entrants and terminal-age exits in the forward law. Current births feed the maintained births/2.1 queue and can affect new-household entry twenty years later. This is distinct from counting children as resident persons today.

The saved tables cannot recover the number of children inside these households, detailed parity counts, spouses/non-head adults, or resident-population totals. That requires richer state/person accounting. The young age map and the maintained household-formation/child-maturation rules remain limitations.

Verification: 94 input files checked against stage/graph manifests, 750 checks passed; all 44 case/date/age-group records are in counts.csv. counts_audit.json records hashes, definitions and limitations.

## Worked impact example: 100,000 initial households

These are net reform-minus-baseline differences, not identified individual switchers. Using raw explicit births avoids mixing household categories with the extra top-code adjustment.

| Component | Net additional count |
|---|---:|
| First births / new ever-parent households | 28.6561 |
| Later births returning a household to the dependent proxy category | 5.1814 |
| Additional births in already-dependent proxy households, residual | 8.8905 |
| All raw explicit births | 42.7280 |
| Additional households with a dependent proxy child | 33.8375 |
| Additional total household decision units | 0.0000 |

The extra dependent-household count equals new ever-parent households plus returns from zero dependent state. The remaining births increase the number of modeled children inside households already in that category. The independently summed first-birth series equals the increase in ever-parent household mass within numerical precision. This impact identity does not describe later-date cumulative cohort flows.

## Entry handoff and the baseline trajectory

The field entry_flow_E always reports the youngest pre-fertility age-cell stock, not an independently measured annual entry flow. In 2023 that cell inherits the historical age-distribution bridge. From 2027 onward it equals the cohort inserted by the previous date's renewal queue. Its increase from 0.0143371220 to 0.0638508304 is therefore a 4.4535-fold change across that handoff. Every later youngest-cell stock matches the previous date's entrant_flow_next within numerical precision.

This is not monotonic contraction from the starting date. Total baseline households rise 1.2301% between 2023 and 2027, peak in 2027, then decline at every subsequent date; they first fall below their 2023 level in 2035. The first increase exactly equals the inserted cohort 0.0638508304 minus reported exits 0.0499597386.

Dependent-proxy household counts already fall 6.3870% between 2023 and 2027, and their share falls 2.7401 percentage points. Their count and share both peak in 2023 and decline at every reported date. A growing total household count in 2027 therefore coexists with fewer households in the dependent-proxy state.

| Year | Total baseline HH per 100 initial | Dependent-proxy HH per 100 initial | Dependent-proxy share % | Reform effect on share, pp |
|---|---:|---:|---:|---:|
| 2023 | 100.0000 | 36.4149 | 36.4149 | 0.0338 |
| 2027 | 101.2301 | 34.0891 | 33.6749 | 0.0614 |
| 2031 | 100.9124 | 31.8843 | 31.5960 | 0.0842 |
| 2035 | 99.2017 | 29.7102 | 29.9493 | 0.1031 |
| 2039 | 96.6285 | 27.6028 | 28.5660 | 0.1187 |
| 2043 | 93.5753 | 25.6914 | 27.4554 | 0.1299 |
| 2047 | 89.8955 | 24.1003 | 26.8093 | 0.1420 |
| 2051 | 86.0467 | 22.9584 | 26.6813 | 0.1529 |
| 2055 | 82.2818 | 21.8116 | 26.5084 | 0.1632 |
| 2059 | 78.1219 | 20.6489 | 26.4316 | 0.1736 |
| 2063 | 74.0860 | 19.5240 | 26.3531 | 0.1835 |

The reform raises the dependent-proxy share at every reported date, by 0.03384 pp initially and 0.18348 pp in 2063. These policy differences are separate from the much larger shared baseline change. The entry discontinuity identifies the inherited empirical-to-queue handoff; it does not by itself establish a coding error, newly created households at childbirth, or a policy-induced migration response.
