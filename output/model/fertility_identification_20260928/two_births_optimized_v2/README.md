# Corrected two-birth experiment

Author authorized additional work on this experiment before overnight planning.
The reference is **2007 stationary reference — block0506, September 28 verified export**.
Version 1 and its failed checkpoint remain immutable. This version corrects only
inner probability arithmetic: shifted exponentials, identical inclusive values,
identical infeasibility mask and 1e-12 normalization tolerance. The earlier
saved-array audit found the error only at 22 unreached menus; the correction
passed its regression and eight existing solver tests (Torch 18757604).

Economic experiment remains unchanged: after success, one optional further
attempt using the existing later-birth Gumbel scale/inclusive value and an
independent draw from the same age-specific conception probability. At most two
births per four-year cell, three children in the state space. Linear within-cell
age projection remains. These assumptions are experimental, not adopted.

All ten reference calibration coordinates remain fixed. Child benefit is
adjusted to the completed-fertility target 2.1 with the original renewal check.
The numerical starting guess uses the authenticated failed point's already
normalized child benefit. Income, entry, other fiscal/housing primitives,
targets, weights, bounds and scientific gates remain unchanged. No shared
source changes, transition, search or promotion.

One smoke: 16 synthetic tests plus a full flag-off reference solve, 20-minute
internal cap. Following lead review and source-pin verification, one corrected
normalized evaluation: max8 stationary solves,
20-minute supervisor/22-minute Slurm cap. The native normalization does not
independently enforce the supplied 1050-second deadline; the supervisor owns
the actual time limit. Expected 2 minutes for smoke and
3–8 minutes for evaluation using observed 131–239-second stationary solves;
queue waiting is additional. No automatic retry. Save heartbeat every15seconds,
latest/best completed-case summaries, full14fit/31parameter rows and17plots.

Smoke18758966 passed16tests plus the full flag-off reference replay in107seconds;
all12 arrays match exactly. Independent source review found no launch blocker.

## Verified normalized experiment

Torch **18759222** completed in 318 seconds, using one stationary evaluation
(242 seconds) at the previously normalized numerical starting guess. All native
scientific gates and the additional extra-choice cache audit pass. Maximum
probability-menu sum error is 2.220e-16; reached invalid-menu mass is zero.
The original 1e-12 probability tolerance is unchanged. Saved-table analysis
**18759306** verified all targets, weights, parameter restrictions, lifecycle
projections, source identities and 17 standard plot names/hashes without any
model solve. All 17 plots received visual review. The inherited upper-wealth
housing discontinuities and retirement profiles remain diagnostic caveats;
standard probability panels show the outer attempt only.

Allowing the second birth raises children among mothers by age 25 from 1.190
to 1.448, but motherhood falls from 45.011% to 35.656%. Their product, children
per woman, falls from 0.535 to 0.516 (target 0.810). Mean first-birth age rises
from 25.933 to 27.399 (target 25.976). Child benefit falls from the actual saved
reference value 0.136 to 0.090 to meet the completed-fertility target 2.100.
The ten calibration coordinates remain held at the reference estimates.
This is a verified experimental equilibrium, not a recalibration or an adopted
specification. It does not establish infeasibility after recalibration.

The symmetric arithmetic decomposition of the change in children by 25 is
+0.104 from children among mothers and −0.123 from the motherhood share. These
are accounting components, not causal effects. A separately authorized
one-equilibrium control holds the saved reference child benefit fixed to
separate the birth-rule change from the benefit adjustment, including equilibrium
responses. Torch **18760270**, 15-minute supervisor/17-minute allocation, one
worker. The control reports any completed-fertility and demographic renewal
misses; it cannot be treated as a replacement-stationary calibration candidate.
All other household, market, fiscal and probability gates remain in force.

### Identity

Reference checkpoint SHA256:
`b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d`.
Reference source manifest:
`07d84336a3112b251afe505908113d9c00585b91c34bd50f0dee108435db496d`.
Reference parameter table:
`4d229fe18ee9c43a8a709413a10690a845d0cf7028efe57925d4f4b1d7c8d036`.
Contract: `68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf`.
Objective: `04528ece6513e4da9436bb8b35b54cd7e9cecd3a6695d050a0a7aea5a71d5e60`.
Corrected experimental checkpoint:
`7e561325fe5c9ff285a9380e381009b1916aaaa160ebca981335da1d988ac6ba`.
Effective overlay manifest:
`2ea39afa424643100ba2c6882c245990e7957a65cf824d5f6913e21a4c0e6a2c`.
Effective solver:
`4e262ecd299f8852aba00a32cee40021890a97b3a34c381a0b375b9fa4a8ba7e`.

The actual reference child benefit is 0.1355551166583114. The historical value
0.14281100340255604 was a numerical starting guess, not this checkpoint's
parameter. The fixed-benefit control reads the authenticated checkpoint itself.
All checkpoints and generated full source copies stay on Torch.

### Full target fit: normalized two-birth experiment

Loss **394.225**, versus **19.581** under the original birth rule at the same
ten coordinates. Original objective and weights. Gaps are model minus target.
Full-precision side-by-side tables are in
`run_v1/saved_analysis/target_fit_comparison.csv` and
`run_v1/saved_analysis/parameters_comparison.csv`.

#### Targeted rows

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Childlessness, ages 40–44 | 0.198 | 0.242 | 0.043 | 35532.304 | 66.937 |
| Exactly one child among mothers, ages 40–44 | 0.214 | 0.190 | -0.024 | 26952.821 | 15.092 |
| Mean age at first birth | 25.976 | 27.399 | 1.423 | 139.828 | 283.043 |
| Wealth / annual earnings | 6.927 | 6.314 | -0.613 | 7.595 | 2.850 |
| Annual bequests / wealth | 0.007 | 0.007 | -2.244e-04 | 5165289.256 | 0.260 |
| Mean rooms | 5.729 | 5.825 | 0.096 | 128.021 | 1.173 |
| Ownership, ages 30–55 | 0.676 | 0.655 | -0.021 | 2339.362 | 1.032 |
| First-birth rooms response | 1.465 | 1.729 | 0.264 | 137.565 | 9.623 |
| Recent-parent ownership gap | 0.128 | 0.142 | 0.014 | 27055.823 | 5.625 |
| Children by age 25 (capped at 3) | 0.810 | 0.516 | -0.293 | 100.000 | 8.590 |

#### Normalization target

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Completed fertility (normalization) | 2.100 | 2.100 | 2.227e-05 | — | — |

#### Untargeted checks

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| First births at age 30+ | 0.249 | 0.317 | 0.067 | 0.000 | 0.000 |
| Older wealth/income, p90/p50 | 3.516 | 3.070 | -0.445 | 0.000 | 0.000 |
| Rooms gap: 3+ versus 1–2 resident children | 0.385 | 0.286 | -0.100 | 0.000 | 0.000 |

### All parameter estimates, original search bounds and external restrictions

No parameter coordinate is estimated in this experiment. The ten listed search
bounds describe the reference calibration; those ten values are held fixed here.
Both fertility choice scales retain the reference near-lower-bound flags.
Child benefit alone is adjusted to meet its normalization target.

| Parameter | Value | Lower | Upper | Near bound | Status in this experiment |
|---|---:|---:|---:|---|---|
| H0 | 6.294 | 0.200 | 80.000 | False | held at reference estimate during two-birth diagnostic |
| beta_annual | 0.963 | 0.940 | 0.990 | False | held at reference estimate during two-birth diagnostic |
| chi | 1.094 | 0.100 | 5.000 | False | held at reference estimate during two-birth diagnostic |
| first_birth_fixed_cost | 0.621 | 0.000 | 8.000 | False | held at reference estimate during two-birth diagnostic |
| kappa_fert | 0.176 | 0.020 | 50.000 | True | held at reference estimate during two-birth diagnostic |
| kappa_fert_continuation | 0.332 | 0.020 | 50.000 | True | held at reference estimate during two-birth diagnostic |
| theta0 | 0.125 | 0.000 | 8.000 | False | held at reference estimate during two-birth diagnostic |
| delta_alpha_jump | 0.135 | 0.000 | 0.250 | False | held at reference estimate during two-birth diagnostic |
| child_benefit_curvature | 0.061 | 0.000 | 0.800 | False | held at reference estimate during two-birth diagnostic |
| tenure_choice_kappa | 0.012 | 0.001 | 0.100 | False | held at reference estimate during two-birth diagnostic |
| psi_child | 0.090 | — | — | — | normalized to completed fertility 2.1 |
| child_benefit_CRRA_coefficient | 0.085 | — | — | — | derived from normalized one-child benefit |
| theta1 | 0.008 | — | — | — | fixed external restriction |
| sigma | 2.000 | — | — | — | fixed |
| alpha_cons | 0.733 | — | — | — | fixed CEX childless expenditure share |
| delta_alpha | 0.000 | — | — | — | fixed zero later-child loading |
| h_P | 0.000 | — | — | — | no housing floor |
| utility_reference_rent | 0.110 | — | — | — | fixed substantive utility normalization |
| q_annual | 0.020 | — | — | — | author-retained 2% annual real rate |
| financed_share | 0.800 | — | — | — | inherited credit contract |
| housing_supply_elasticity | 0.630 | — | — | — | fixed provisional external mapping |
| payroll_tax | 0.080 | — | — | — | derived from adopted pension ratio |
| pension_period | 0.918 | — | — | — | balanced PAYGO |
| annual_depreciation | 0.014 | — | — | — | adopted |
| period_depreciation | 0.055 | — | — | — | compounded |
| annual_property_tax | 0.011 | — | — | — | adopted |
| period_property_tax | 0.042 | — | — | — | linear period convention |
| selling_cost | 0.060 | — | — | — | retained |
| rental_cap | 6.000 | — | — | — | retained provisional |
| wealth_grid_nodes | 160.000 | — | — | — | retained exact grid |
| income_states | 15.000 | — | — | — | retained B15 |

### Diagnostics and measurement scope

The unchanged 17 standard PNGs are in `run_v1/worker/case/standard_diagnostics/`.
The supplemental fertility comparison is
`run_v1/saved_analysis/supplemental_lifecycle.png`; its three panels show children
per woman, motherhood, and children among mothers in identical five-year windows.
All observed counts are capped at three children. The CPS 2004/2006 profile is
cross-sectional; NCHS first births use 2003–2006 period counts. Neither is a
single cohort trajectory. The model approximates the 2007 stationary distribution;
this is not a 2007–2023 transition result. Both births use the retained linear
within-period age projection, with a common event-time proxy, not separately
modeled birth dates. That approximation remains unresolved.

To regenerate the supplemental packet from saved tables on Torch (no solves):
`python output/model/fertility_identification_20260928/two_births_optimized_v2/analyze_result.py --plot`.
The native evaluator generates the full standard packet as part of each
`run.sh evaluate NEW_FOLDER` invocation; this requires a new explicitly bounded
model run and the matching passed smoke receipt, not an automatic restart.

## Fixed-benefit control completed

Torch **18760270** passed its specified diagnostic checks in 272 seconds, one
197-second equilibrium solve. Saved comparison **18760896** passed in 3 seconds,
zero solves/imports/checkpoint reads: it authenticates both experiments against
the reference, same effective household/reporting source hashes, all 14 target
rows, all 31 parameter/restriction rows and all 17 standard PNG names/hashes.
An independent visual review inspected all 17 control plots; the lead additionally
checked fertility, ownership and market quantities. No new obvious visual failure;
existing extreme-wealth housing/ownership and retirement-profile caveats remain.
Market residual is 1.313e-6, within the unchanged native gate. Extra-probability
error is 2.220e-16 and reached invalid-menu mass is zero.

| Fertility measure | Data / target | One birth per period | Two births, fixed benefit | Two births, benefit normalized |
|---|---:|---:|---:|---:|
| Children by age 25 | 0.810 | 0.535 | 0.742 | 0.516 |
| Mothers by age 25 (%) | 45.725 | 45.011 | 49.186 | 35.656 |
| Children among mothers by age 25 | 1.770 | 1.190 | 1.509 | 1.448 |
| Mean age at first birth | 25.976 | 25.933 | 26.020 | 27.399 |
| Completed fertility (normalization) | 2.100 | 2.100 | 2.619 | 2.100 |

All ten calibration coordinates are the original reference estimates in all
three cases. At fixed child benefit, the experimental birth rule closes 75.444%
of the age-25 children gap, but completed fertility exceeds target by 0.519.
It also tracks average children much better through ages 30–34, then overshoots:
at ages 40–44, children per woman are 2.146 versus data 1.718. Reducing child
benefit to meet the completed-fertility target removes the early gain. In that
comparison the symmetric accounting components are −0.200 from motherhood and
−0.026 from children among mothers. The benefit change includes equilibrium
responses; these two accounting components are not separate causal estimates.

The fixed-benefit control has entry E=0.062 and birth-generated entry B=0.077;
(E−B)/E=−24.691%. The native demographic renewal check was evaluated and its
specific failure recorded. It was deliberately not required for this diagnostic.
This control is not a demographic steady state or a calibration candidate.
The completed-fertility target remains 2.1 in the tables, not a changed target.

This supports retaining the extra birth opportunity as a candidate mechanism,
but provides no evidence yet that all targets can be jointly fitted. A useful
next local diagnostic would vary the fixed cost of starting a family while
retaining completed-fertility normalization, then inspect motherhood, conditional
children and first-birth timing together. The original-model Jacobian does not
identify derivatives in this changed model. No further run is authorized by this
note, no overnight plan is launched, and no reference is promoted.

Control checkpoint SHA256:
`78c82a13fe5c414435de8f2fed3e62b7ef16b67301366d64f503f7d842dffb86`.
Control driver SHA256:
`11134c6c1f44b1c2746029c2adde3768416069fdb4e099b67399b87be262effd`.
Its child benefit is the exact saved reference value 0.1355551166583114.
The conditional diagnostic score is 611.340, excluding the unmet normalization
row; it cannot be ranked as an admissible calibration improvement.
Full-precision three-case evidence is in `fixed_benefit_v1/saved_comparison/`,
including the complete 31-row parameter comparison with original bounds and
restrictions. All ten originally estimated values and both near-bound flags are
identical to the normalized table above. Here child benefit and its derived CRRA
coefficient are held at the reference values; none is normalized or re-estimated.

### Control: targeted rows

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Childlessness, ages 40–44 | 0.198 | 0.110 | -0.088 | 35532.304 | 275.143 |
| Exactly one child among mothers, ages 40–44 | 0.214 | 0.132 | -0.082 | 26952.821 | 179.937 |
| Mean age at first birth | 25.976 | 26.020 | 0.044 | 139.828 | 0.271 |
| Wealth / annual earnings | 6.927 | 6.267 | -0.660 | 7.595 | 3.307 |
| Annual bequests / wealth | 0.007 | 0.007 | -2.101e-04 | 5165289.256 | 0.228 |
| Mean rooms | 5.729 | 5.894 | 0.165 | 128.021 | 3.467 |
| Ownership, ages 30–55 | 0.676 | 0.652 | -0.024 | 2339.362 | 1.403 |
| First-birth rooms response | 1.465 | 1.654 | 0.189 | 137.565 | 4.936 |
| Recent-parent ownership gap | 0.128 | 0.055 | -0.072 | 27055.823 | 142.196 |
| Children by age 25 (capped at 3) | 0.810 | 0.742 | -0.067 | 100.000 | 0.453 |

### Control: normalization target, deliberately unmet

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Completed fertility (normalization) | 2.100 | 2.619 | 0.519 | — | — |

### Control: untargeted checks

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| First births at age 30+ | 0.249 | 0.235 | -0.014 | 0.000 | 0.000 |
| Older wealth/income, p90/p50 | 3.516 | 3.060 | -0.456 | 0.000 | 0.000 |
| Rooms gap: 3+ versus 1–2 resident children | 0.385 | 0.238 | -0.147 | 0.000 | 0.000 |

Run saved comparisons on Torch with
`python output/model/fertility_identification_20260928/two_births_optimized_v2/analyze_control.py --plot`.
Native standard-packet function is
`evaluator.rt['audit'].standard_diagnostics(packet, output, validate_production_young=False)`
after `driver.load_runtime` and `driver.install_overlay`; use the saved packet on
Torch and a new output directory. This function renders the same 17 panels
without a new equilibrium solve. Do not import or render on the Mac.


The optional supplemental four-series plot was rendered and its twelve plotted
series checked on Torch **18761226** (5 seconds, zero solves). Lead visual review
passes. It is `fixed_benefit_v1/saved_comparison/supplemental_lifecycle_control.png`.
The earlier no-plot analysis receipt/source are preserved alongside it; the new
receipt records the rendering source and output hash. Both standard 17-plot sets
remain intact. The full project source remains unchanged.
