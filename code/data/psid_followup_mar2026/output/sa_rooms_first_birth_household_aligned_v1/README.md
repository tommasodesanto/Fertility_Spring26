# Corrected PSID first-birth rooms target

The paper-facing target is the change in occupied rooms from event time -1 to
event time +3 around a first birth. This four-year contrast is the closest data
analogue to one model period. The corrected estimate is
`0.720246262382` with clustered standard error
`0.0852600513386` (inverse-variance weight
`137.5652749`).

The Sun--Abraham event study uses 49,457
woman-year observations for 4,112 women. Every
observation is unique at the household-year level, and the primary sample
excludes physical dwellings containing multiple PSID family units. The sample
contains current women age 18 or older who are the reference person or
spouse/partner, have observed age and positive PSID longitudinal weight, and
live in a valid single-family-unit dwelling. Women whose full relationship
history confirms zero children are the comparison group. The regression
includes person and survey-year fixed effects, age and education controls,
applies `IW` as a probability weight, and clusters by person ID.

`ACTUALROOMS_` is shifted forward by one observed interview within individual
before the vintage-specific non-room codes are removed. This repairs the
one-wave extraction error in the prior target. The earlier unweighted,
last-treated-control estimate is withdrawn.

The treatment date is the first biological child across all twenty recorded
child relationships. Known parents without a biological-birth year and people
with unknown child-count histories are excluded from the control group.
Event time -2 is the unique omitted period; event time -6 is explicitly
estimated. The reported standard error uses the full covariance between the
+3 and -1 coefficients. Earlier leads reject a flat prepath, so the complete
event-study profile remains a required diagnostic. The contrast may include
adjustment associated with a later birth and is not a mechanically pure 0-to-1
space requirement. Its use in calibration is therefore a provisional mapping
to the model's one-period birth-versus-no-birth branch, not a causal claim.

Contract ID: `psid_first_birth_rooms_household_aligned_sa_v1_20260817`. Source, code, log, and output hashes
are recorded in `metadata.json`.
