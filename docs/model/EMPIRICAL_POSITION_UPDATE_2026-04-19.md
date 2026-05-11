# Empirical Position Update

Updated interpretation after the ACS renter-only elasticity runs.

## Keep

- `P2a` should be treated as an **ACS descriptive fact**, not a PSID headline.
- The right statement is:
  `Among renters, parent households have lower income elasticity of housing quantity/services.`
- Main empirical support:
  ACS renters, `log rooms` on `log income`:
  - childless slope = `0.129`
  - parent slope = `0.102`
  - parent interaction = `-0.027`
  - source: `output/acs_renter_elasticity_v1/acs_renter_elasticity_summary_v1.csv`
- Recommended theory-slide wording:
  `Parents are less income-elastic in housing.`
  Clarification in notes/text:
  `In the data, this shows up most cleanly among renters using ACS rooms as the housing-quantity proxy.`

## Drop

- Do **not** keep:
  `Parents are less responsive to local rent.`
- ACS renter `rooms` on local rent-per-room gives:
  - childless slope = `-0.298`
  - parent slope = `-0.301`
  - interaction = `-0.003`, `p = 0.226`
  - source: `output/acs_renter_price_elasticity_v1/acs_renter_price_elasticity_summary_v1.csv`
- So `P2b` is not supported in the current data.

## Separate Keep From Lemma 2'

- Keep the **cross-metro** prediction as a separate ACS fact from `lemma2prime_v2`.
- Implemented result:
  - center-residence gap on log center/periphery rent gap: slope = `6.23`, `p = 0.237`
  - within-CBSA center-destination gap on log center/periphery rent gap: slope = `23.53`, `p = 0.00047`
  - source: `output/acs_income_crossmetro_v1/acs_cross_metro_regressions_v1.csv`
- Interpretation:
  the mover-destination version is clearly supportive; the stock center-residence version is positive but imprecise.

## Do Not Claim

- Do not claim:
  - `First-birth relocation rate decreasing in pre-birth wealth`
  - `Flight moves come from low/middle-wealth households`
- Current PSID wealth-by-quintile event studies are not monotone.
- ACS income-proxy flight results are mixed and not clean enough for a slide claim.

## Recommended Slide Language

- `Parents are less income-elastic in housing quantity among renters (ACS rooms).`
- `We do not currently show a clean parent-vs-childless difference in local-rent responsiveness.`
- `We do not currently have clean evidence that family flight is concentrated in low-/middle-resource households.`
- `Metros with larger center-periphery rent gaps exhibit stronger parent sorting on the mover destination margin.`
