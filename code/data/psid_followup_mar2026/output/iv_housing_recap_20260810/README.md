# Fertility and Housing Estimates: August 10 Recap

## Numbers to keep straight

The first table is the main first-birth event-study evidence. The second table
contains the preferred corrected five-year fertility-IV triangulation. The
third table records the two exploratory robustness designs run afterward.
These are not identical estimands, so they should not be averaged together.

## Main first-birth event studies at event time +3

| Outcome | Coefficient (SE) | Interpretation | Status |
|---|---:|---|---|
| Rooms | 0.805 (0.167) | Additional rooms relative to event time -2 | Household-ID-FE source, but under the active 0/99 room-code hold |
| Moved for size | 0.020 (0.006) | 2.0 percentage-point increase in the interview-level probability | Household and year FE; current reproducible TWFE estimate |
| Homeownership | 0.113 (0.009) | 11.3 percentage-point increase in ownership | Legacy household-ID-FE output; current source provenance is less clean |
| Moved from renting to owning | Not reported | Historical transition was constructed with an invalid calendar lag | Withdrawn |

## Corrected five-year endpoint 2SLS: preferred IV triangulation

| Outcome | Twins: second-child IV | Same sex: third-child IV |
|---|---:|---:|
| Rooms change | 0.578 (0.852), p=0.497 | -1.229 (1.053), p=0.243 |
| Moved for size | 0.092 (0.252), p=0.716 | 0.104 (0.347), p=0.765 |
| Ever owner, baseline renters | 0.098 (0.303), p=0.747 | 0.275 (0.326), p=0.398 |
| Moved to ownership, baseline renters | 0.051 (0.282), p=0.856 | 0.335 (0.249), p=0.179 |

These remain the cleanest IV magnitudes because the treatment and five-year
housing windows are aligned and the design uses a nearby pre-birth baseline.
None is individually significant. The most favorable same-sex reduced form is
a 4.3 percentage-point increase in moving to ownership (SE 3.0 points,
p=0.146).

## Exploratory IV robustness

| Design and outcome | Twins 2SLS (SE), p | Same-sex 2SLS (SE), p |
|---|---:|---:|
| Full-panel rooms | 0.706 (0.424), p=0.096 | -3.114 (2.353), p=0.186 |
| Full-panel moved for size | -0.033 (0.104), p=0.751 | -0.048 (0.355), p=0.892 |
| Full-panel ownership | 0.238 (0.251), p=0.344 | 0.210 (0.661), p=0.751 |
| Full-panel moved to ownership | -0.022 (0.067), p=0.747 | 0.189 (0.298), p=0.526 |
| Post-only rooms mean, exact age/race FE | 0.776 (0.828), p=0.349 | -1.895 (1.823), p=0.299 |
| Post-only moved for size, exact age/race FE | 0.124 (0.220), p=0.571 | 0.171 (0.468), p=0.714 |
| Post-only ownership mean, exact age/race FE | 0.330 (0.179), p=0.065 | 0.081 (0.372), p=0.827 |
| Post-only moved to ownership, exact age/race FE | -0.094 (0.089), p=0.287 | 0.165 (0.224), p=0.460 |

The conventional quadratic-age post-only twins ownership specification gives
0.392 (SE 0.178, p=0.028). Exact maternal-age indicators move its p-value to
0.065, so it is a positive signal rather than robust conventional significance.
The full-panel twins rooms result is likewise suggestive, but its weighted
pre-period joint test has p=0.082.

## Short scientific summary

1. The main first-birth event studies show sizable increases in rooms,
   space-related moving, and ownership.
2. The twins IV repeatedly gives a positive rooms magnitude of roughly
   0.6--0.8 rooms, but precision depends strongly on the construction.
3. Ownership is positive in several twins specifications; the post-only result
   is closest to significance but is sensitive to maternal-age controls.
4. Same-sex results support ownership-related adjustment more than rooms. The
   rooms sign is consistently negative in the corrected IV exercises, and the
   exclusion restriction is questionable because sex composition can directly
   affect room sharing.
5. Treat the IVs as triangulation. The main event studies remain the primary
   empirical objects, with the rooms estimate still under its coding hold.
