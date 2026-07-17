# Intergen Quality-Adjusted Size Target Audit

Date: 2026-06-23

Status: empirical target diagnostic. This does not change model logic, solver
logic, or the active SMM target set.

## Question

The current one-market intergen calibration targets raw ACS rooms for childless
prime-age owners and renters:

\[
\bar h^{data}_{owner}=6.224,\qquad
\bar h^{data}_{renter}=3.805,\qquad
\Delta h^{data}_{owner-renter}=2.419.
\]

The model has no geography and no explicit housing-quality or structure-type
dimension. The audit asks whether the raw owner-renter room gap is mostly a
geography/quality artifact, and whether a quality-adjusted size object would be
a more appropriate hard target.

The measurement premise is standard. BLS treats shelter as the service provided
by a unit and quality-adjusts rent/OER for changes in unit services and aging:
<https://www.bls.gov/cpi/factsheets/owners-equivalent-rent-and-rent.htm>. BEA's
ACS-based housing-services work imputes owner rental equivalence from tenant
rents using structure type, rooms, bedrooms, age of structure, and geography:
<https://www.bea.gov/sites/default/files/2019-11/improving-measures-of-national-and-regional-housing-services-us-accounts.pdf>.
OECD's hedonic housing-price guidance emphasizes that location is a key
determinant of housing value:
<https://www.oecd.org/content/dam/oecd/en/publications/reports/2011/02/hedonic-price-indexes-for-housing_g17a1f54/5kghzxpt6g6f-en.pdf>.
RAND's hedonic housing-services work separates location, space, interior
quality, and exterior quality:
<https://www.rand.org/pubs/reports/R2450.html>.

This supports auditing quality-adjusted size. It does not by itself justify
dropping the raw room target.

## Construction

Script:
`code/data/mms_center_periphery/audit_intergen_quality_adjusted_room_targets.R`

Outputs:
`code/data/mms_center_periphery/output_intergen_quality_adjusted_room_targets/`

Sample matches the existing ACS room-target builder:
ACS/IPUMS `extract27`, household heads, 2012--2023, matched MMS metro sample,
middle collapsed to center, age 30--55, childless, owner or renter.

For each control set \(X\), the script forms exact cells \(c\), keeps common
support cells with both owners and renters, computes tenure-specific weighted
mean rooms within each cell, and averages the within-cell owner-renter gap:

\[
G_X =
\sum_{c\in\mathcal C_X}
\omega_c
\left(\bar h_{owner,c}-\bar h_{renter,c}\right),
\qquad
\omega_c =
\frac{W_{owner,c}+W_{renter,c}}
{\sum_{d\in\mathcal C_X}(W_{owner,d}+W_{renter,d})}.
\]

The model comparison values are from the current June 23 diagnostic best
14-moment/scalar point, not a calibrated benchmark:

\[
h^{model}_{owner}=6.053,\qquad
h^{model}_{renter}=4.560,\qquad
\Delta h^{model}_{owner-renter}=1.493.
\]

## Main Results

`M-T` means model minus target.

| Target object | Target owner | Model owner | M-T owner | Target renter | Model renter | M-T renter | Target gap | Model gap | M-T gap | Support |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Raw current target | 6.224 | 6.053 | -0.171 | 3.805 | 4.560 | 0.755 | 2.419 | 1.493 | -0.926 | 1.000 |
| Same year, metro, MMS location, age | 6.056 | 6.053 | -0.003 | 3.903 | 4.560 | 0.657 | 2.152 | 1.493 | -0.660 | 1.000 |
| Add structure type | 5.498 | 6.053 | 0.555 | 4.625 | 4.560 | -0.065 | 0.872 | 1.493 | 0.620 | 0.903 |
| Add structure type and vintage | 5.666 | 6.053 | 0.387 | 4.873 | 4.560 | -0.313 | 0.793 | 1.493 | 0.700 | 0.616 |
| Add bedrooms too | 5.361 | 6.053 | 0.692 | 5.102 | 4.560 | -0.542 | 0.259 | 1.493 | 1.234 | 0.424 |

Location composition is real but not large enough to resolve the calibration
tension. In the childless age-30--55 sample, `58.4%` of renter weight is in MMS
center locations, versus `41.3%` of owner weight. Holding year, metro, MMS
location, and age fixed lowers the owner-renter room-gap target from `2.419` to
`2.152`, but the current diagnostic model gap is still only `1.493`.

Structure composition is much larger. Once units are compared within the same
structure type, the gap falls to `0.872`; adding building vintage gives `0.793`.
But this is not automatically a better hard target. Structure type is not just
nuisance quality. It is plausibly part of the tenure product: owners are much
more likely to occupy detached single-family units, while renters are much more
likely to occupy multifamily units. Controlling it away may remove exactly the
housing ladder margin that \(\chi\), \(\bar h_0\), \(\bar h_{\mathrm{jump}}\),
and the owner-rung menu are supposed to discipline.

The bedroom-controlled specification is diagnostic only. Bedrooms are a
physical size margin, so conditioning on bedrooms over-controls the object we
are trying to target.

## Interpretation

The audit does not say "the room-gap target was wrong." It says:

1. Geography/year/age adjustment is defensible and modest. If the no-location
   model needs a cleaner hard target, `2.152` is a reasonable candidate
   replacement for the raw owner-renter room gap `2.419`.
2. This does not solve the current fit. Under the location-adjusted target, the
   model still understates the owner-renter room gap by `0.660` rooms, and the
   miss remains mostly on renters: model renters are `4.560` rooms versus a
   location-adjusted target of `3.903`.
3. Structure adjustment changes the sign of the fit, but it is too aggressive
   to use as a hard target without adding a replacement moment for the same
   economic block. Otherwise the calibration would stop disciplining the
   owner-product/tenure-ladder distinction.
4. The empirical result supports the user's concern about missing geography and
   housing quality, but the main economic message is more specific: the current
   scalar room target conflates physical size with structure/quality, and the
   model has no separate starter-owner versus large-owner/quality ladder.

## Identification Implications

If the raw room-gap moment is replaced by the location/year/age-adjusted gap,
identification is mostly preserved: it remains an owner-renter housing-services
separation moment for the housing-need and tenure blocks.

If the structure-adjusted gap is used instead, identification is not preserved
unless another hard moment replaces the removed structure/product information.
Candidate replacements include:

- owner and renter shares in detached/single-family units, if the model gains a
  structure or quality state;
- owner share with rooms \(\ge 6\) and renter share with rooms \(\ge 6\), kept
  as raw distributional objects;
- an explicit starter-owner versus large-owner rung share, if the owner ladder
  is redesigned;
- an external restriction on the owner-rung/quality menu.

Without one of these, using the structure-adjusted gap would make the room block
easier by erasing the product margin the model currently fails to represent.

## Recommended Next Step

Do not launch a new global search yet solely on the structure-adjusted gap.
First rescore the existing June 23 frontier under two diagnostic alternatives:

1. raw room-gap target `2.419`;
2. location/year/age-adjusted room-gap target `2.152`.

Keep the other hard targets fixed and report target, model, and gap side by
side. If the joint frontier still cannot deliver young ownership, old exit, and
room separation under the `2.152` gap, the next productive step is a model/menu
change: separate starter ownership from large/quality ownership, rather than
declaring the size moment solved.
