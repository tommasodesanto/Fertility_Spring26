# Active ACS Room-Target Receipt

Generated: `2026-08-17 01:49:39 EDT`
Bootstrap: `B=1000`, seed `20260705`, 42 `MET2013` clusters resampled with replacement.

## Result

- `prime30_55_parent_3plus_minus_1to2_mean_rooms`: target `0.36769955881000`, reproduced `0.367699558813823`, metro-bootstrap SE `0.069002912`; active declared synthetic SE `0.018384978` (weight `2958.514988`).
- `aggregate_mean_occupied_rooms_18_85`: target `5.779970481941968`, reproduced `5.779970481941968`, metro-bootstrap SE `0.113293163`; active declared synthetic SE `0.288998524` (weight `11.973159`).
- Same-draw covariance: `0.003818430`; correlation: `0.488442801`.

These metro-bootstrap SEs are a project-design uncertainty measure, not an official ACS replicate-weight variance estimate. The extract does not carry ACS replicate weights. They are reported for audit and are not substituted for the active calibration's declared synthetic five-percent SEs.

## Exact empirical objects

The source is the pooled 2012--2023 IPUMS ACS `extract27.dta`. The base sample keeps identifiable `MET2013` observations in households (`GQ` 1 or 2), the first person in the household (`PERNUM=1`), positive `HHWT`, owner or renter tenure, and positive literal `ROOMS`. It joins 2012--2021 observations to the 2010-PUMA lookup and 2022--2023 observations to the 2020-PUMA lookup on state, PUMA, and CBSA. Matched center and periphery PUMAs are retained and the lookup's middle category is collapsed to center.

Mean rooms additionally keeps ages 18--85. The family-size contrast keeps ages 30--55 with at least one own child in the household and a youngest own child younger than 18; it compares `NCHILD>=3` with `NCHILD` 1--2. `NCHILD` is current own children in the household, not completed parity. Both estimators use `HHWT`.

The pinned geography build summary records a 0.30 core-tract population target, 0.50 center-PUMA cutoff, 0.10 periphery-PUMA cutoff, 51 retained cities, and 1,361 PUMA-city rows in each period lookup. The realized target sample contains 42 distinct `MET2013` codes.

## Why the older moment_se.csv is NA

The committed `code/data/moment_standard_errors/output/moment_se.csv` was last regenerated with the PSID-only source branch. Its README and timing log contain no ACS load or bootstrap entry, and its reproduction log says `ACS source not run`. Thus the NA is a run-selection artifact, not a failed ACS point-estimate gate. The older 14-moment harness also does not contain the active aggregate mean-rooms row.

## Provenance limitation

This receipt independently pins the raw extract, lookup files, relevant builders, and the analysis cache, and reproduces the target formulas from the pinned cache. The legacy cache predates a source-hash manifest, so its relationship to the raw extract is reconstructed rather than cryptographically recorded at cache creation. A future raw rebuild should write these hashes when the cache is created; changing any pinned input makes this script fail before emitting a receipt.

## Files

- `target_receipt.csv`: point estimates, samples, weights, and bootstrap SEs.
- `bootstrap_covariance.csv` and `bootstrap_correlation.csv`: same-draw two-moment matrices.
- `provenance.csv`: SHA-256, byte size, and modification time for every input and builder.
