# One-Market Intergen Candidate Replacement Targets

Date: 2026-06-18

## Status

These are extracted candidate values, not final SMM targets. They are intended
to make the replacement-moment discussion concrete. Before any value enters the
objective, reauditing is required for source file, sample, household unit,
weights, formula, room/bedroom convention, wealth definition, timing, and
model-object mapping.

The recommended baseline convention is: ACS/MMS household heads for
housing/tenure/rooms/renter costs; PSID for lifecycle wealth, old-age
completed-child wealth gaps, and birth-related room increments; nonhousing or
liquid wealth as the primary old-age bequest object, with total net worth as a
robustness check.

## Newly Extracted ACS/MMS Housing Targets

Script:
`code/data/mms_center_periphery/build_intergen_one_market_housing_targets.R`

Output:
`code/data/mms_center_periphery/output_intergen_one_market_targets/intergen_one_market_acs_housing_targets.csv`

Sample: ACS/IPUMS `extract27`, household heads, 2012--2023, matched MMS metro
sample, middle collapsed to center, ages 30--55.

| Candidate moment | Value |
|---|---:|
| Childless owner mean rooms | 6.224 |
| Childless renter mean rooms | 3.805 |
| Childless owner minus renter mean rooms | 2.419 |
| Childless owner share with rooms \(\ge 5\) | 0.799 |
| Childless renter share with rooms \(\ge 5\) | 0.281 |
| Childless owner share with rooms \(\ge 6\) | 0.596 |
| Childless renter share with rooms \(\ge 6\) | 0.138 |
| Childless owner share with rooms \(\ge 7\) | 0.380 |
| Childless renter share with rooms \(\ge 7\) | 0.060 |
| Childless renter mean rent-to-income, \(12 rent / income\) | 0.412 |
| Childless renter median rent-to-income, \(12 rent / income\) | 0.240 |
| Parent owner minus renter mean rooms | 2.400 |

Interpretation: the discrete owner median-room target should be replaced by a
smooth room-distribution object such as the owner share with rooms \(\ge 6\),
the owner-renter mean room gap, or a small vector of tenure-specific room-bin
shares. The rent-to-income mean is high because it is sensitive to low-income
tails; the median is close to the old aggregate cost-share target but measures
only renters.

## Newly Extracted PSID Old-Age Wealth Targets

Script:
`code/data/psid_followup_mar2026/build_intergen_oldage_wealth_targets.R`

Output:
`code/data/psid_followup_mar2026/output/intergen_oldage_wealth_targets/intergen_oldage_wealth_targets.csv`

Sample: PSID `PSIDSHELF_MOBILITY`, 1984--2019, ages 65--75, completed children
from `RELCHIREP`.

| Candidate moment | Value |
|---|---:|
| Old homeownership rate | 0.864 |
| Parent old homeownership rate | 0.872 |
| Childless old homeownership rate | 0.789 |
| Parent minus childless old ownership gap | 0.083 |
| Old nonhousing net worth / income, mean | 6.419 |
| Old nonhousing net worth / income, median | 2.230 |
| Parent old nonhousing net worth / income, mean | 6.510 |
| Childless old nonhousing net worth / income, mean | 5.503 |
| Parent minus childless old nonhousing net worth / income gap | 1.007 |
| Old total net worth / income, mean | 9.672 |
| Old total net worth / income, median | 5.264 |
| Parent minus childless old total net worth / income gap | 1.284 |
| Old liquid savings / income, mean | 0.991 |
| Old liquid savings / income, median | 0.195 |
| Parent minus childless old liquid savings / income gap | -0.105 |

Interpretation: the ownership numbers replicate the broad old-age ownership
object, but the useful replacement for the bequest block is the wealth object.
The mean nonhousing-wealth gap by parent status is positive in ratio form but
negative in dollars, so this target needs special auditing before use. The
candidate primary object is the old nonhousing net-worth-to-income ratio, not
total net worth.

## Existing PSID Housing-Increment Packets

Existing files:

- `code/data/psid_followup_mar2026/output/no_location_family_space_packet_20260523/psid_birth_rooms_event_summary.csv`
- `code/data/psid_followup_mar2026/output/rooms_second_birth_quick_v1/rooms_second_birth_quick_moment_v1.csv`
- `code/data/psid_followup_mar2026/output/rooms_second_birth_quick_variants_v1/rooms_second_birth_quick_variants_v1.csv`

Currently visible candidate values:

| Candidate moment | Value |
|---|---:|
| First-birth room response, post-3 coefficient | 0.664 |
| First-birth room response, post-5 coefficient | 0.843 |
| Second-birth event-study room response, post-3 coefficient | 0.704 |
| Second-birth event-study room response, post-5 coefficient | 0.724 |
| Second-birth quick mean room change, post-3 | 0.488 |
| Second-birth quick mean room change, post-5 | 0.525 |

Interpretation: this block should not be blindly inserted into the objective.
The first-birth target is relatively stable. The second-child housing target
depends on event-study versus quick-window construction, horizon, and treatment
definition. This is a data-target audit item, not a new model-search item.

## Existing PSID Young Wealth Targets

Existing file:
`code/data/psid_followup_mar2026/output/entry_wealth_v1/entry_wealth_candidate_targets_focus_v1.csv`

Selected candidate values:

| Candidate moment | Mean | Median |
|---|---:|---:|
| Young childless, liquid nonhousing net worth / income | 0.870 | 0.141 |
| Young childless renters, liquid nonhousing net worth / income | 0.179 | 0.062 |
| Young pre-birth, liquid nonhousing net worth / income | 0.824 | 0.166 |
| Young pre-birth renters, liquid nonhousing net worth / income | 0.205 | 0.080 |

Interpretation: these are plausible replacements/supplements for the
entry-wealth and savings block, but the exact moment depends on whether the
model's young state is all childless young households, pre-birth households, or
renters near the ownership-access margin.

## Open Audit Items

1. Decide whether the production ACS housing sample is MMS-metro household
   heads, all national household heads, or another no-location population.
2. Choose room-bin shares or mean room gaps as the replacement for owner median
   rooms; do not use a discrete median as the local identifying moment.
3. Decide whether renter rent-to-income alone is an adequate \((\alpha,\chi)\)
   target or whether we also need an owner imputed user-cost object.
4. Reconcile old-age wealth definitions with the model estate object:
   nonhousing net worth, liquid savings, and total net worth give different
   parent-childless patterns.
5. Audit the second-child room-increment target across existing PSID variants
   before treating \(\bar h_n\) as identified by a single number.
