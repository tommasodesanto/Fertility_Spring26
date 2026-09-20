# Rental room-size support diagnostic

This is descriptive evidence for the literal six-room rental cap. It uses the
same 2005–06 ACS household-head stream, active 42-metro sample, HHWT weights,
and household filters as the housing target replay. `ROOMS` codes 1–10 are
kept, with code 10+ represented by capped room bin 9. Tenure codes 1 and 2 are
owner and renter; any other tenure code is retained as `unknown` and excluded
from owner/renter shares. The primary sample is all matched-metro structures;
the ownership sample is age 30–55 with `UNITSSTR` in 3–10. Child rows split
zero versus one-or-more current children using `NCHILD`.

In the primary active-42 sample, renters are 6.94% of valid-tenure household
weight among homes above six rooms, while 6.02% of renter weight is in those
homes. The unweighted figures are 19,425 of 301,008 renter records, or
6.45%. The renter share within the capped 7–9 room bins falls from 9.66% at
seven rooms to 4.31% at nine rooms. Among households with one or more current
children, 7.31% of valid-tenure weight in >6-room homes is renter weight and
9.58% of renter weight is >6 rooms; for zero-current-child households the
corresponding figures are 6.44% and 3.88%. Unknown-tenure records are zero in
this extract.

The ownership 30–55 DUE sample gives 7.20% of >6-room valid-tenure weight as
renter weight and 7.23% of renter weight in >6-room homes (12,224 of 155,432
renter records unweighted). These are descriptive support statistics only;
they do not identify a rental price premium and do not change targets, weights,
or the model specification.

## Reproduction

The full pass, including the unchanged four-target replay, is regenerated with:

```bash
PY=/Users/tommasodesanto/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3
BASE=$(mktemp -d /tmp/rental_support_full.XXXXXX)
$PY code/empirical/housing/build_initial_housing_profile_diagnostic.py \
  --output "$BASE/base" \
  --rental-support-output "$BASE/rental_size_support"
```

The full replay completed in 48.64 seconds over 24 chunks. Its four active-42
target scalars match the prior receipt exactly (absolute tolerance (10^{-10})).
The source is `extract27.dta`; the canonical source SHA-256 is recorded in
`rental_size_metadata.json` as
`edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e`.
The post-fix builder source SHA-256 is
`05d90429c2a98169009aa6124afd91b0ee89c1c081b543177eba6dddb2f882ab`.
The saved summary values are unchanged by these compatibility and naming fixes.
