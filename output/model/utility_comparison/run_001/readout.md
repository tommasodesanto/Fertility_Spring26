# Overnight utility comparison: bounded failure readout

**The planned full eight-hour recalibration was not achieved.** All four
searches stopped in their initial populations; no differential-evolution
round ran. Both floor variants have verified selected results and complete
reports. The share variants' final repetitions and reports could not be
collected after SSH access became unavailable. This evidence does not identify
a preferred utility specification.

The floor variants raise housing needs with children. The share variants
instead change the consumption–housing utility shares with children. Each was
crossed with a linear or diminishing benefit from children currently at home.
These are experimental preference specifications. The common pension ratio
is author-adopted; its update is separate from these utility experiments.

At the same retained structural point, concavity worsens the floor model's
childlessness, exactly-one-child and recent-parent ownership fits. It slightly
improves the family-room gap and mean first-birth age, while leaving the large
mean-room and first-birth housing-response mismatches. This is a conditional
comparison of two saved points, not a comparison of successfully reoptimized
models. Each floor winner was selected from only three distinct successful
points.

## Complete verified floor fits

All thirteen moments are shown below. Completed fertility is the separate
normalization; the other twelve rows are weighted. Percentages are rescaled
only for display. The linked full-precision CSVs contain every target, model
value, gap, weight and loss contribution:
[floor linear](floor_linear/export/target_fit.csv) and
[floor concave](floor_concave/export/target_fit.csv).
The two share-arm fit tables remain **uncollected and unverified**; no values
or ranking are inferred from incomplete monitoring summaries.

| Moment | Target | Floor, linear benefit | Floor, concave benefit |
| --- | ---: | ---: | ---: |
| Completed fertility (normalization) | 2.100 | 2.100 | 2.100 |
| Childlessness, ages 40–44 (%) | 19.828 | 18.279 | 17.737 |
| Exactly one child among mothers, ages 40–44 (%) | 21.366 | 23.126 | 23.564 |
| Mean age at first birth (years) | 25.976 | 26.324 | 26.316 |
| First births at age 30 or older (%) | 24.928 | 24.600 | 24.597 |
| Wealth / annual gross earnings | 6.927 | 6.215 | 6.216 |
| Annual bequest flow / wealth (%) | 0.729 | 0.688 | 0.688 |
| Wealth/income p90 / median, ages 76–84 | 3.516 | 2.910 | 2.909 |
| Mean occupied rooms, capped at nine | 5.608 | 6.543 | 6.544 |
| Homeownership, ages 30–55 (%) | 67.626 | 76.181 | 76.161 |
| First-birth housing response (rooms) | 1.465 | 0.785 | 0.785 |
| Room gap: 3+ versus 1–2 children at home | 0.385 | 0.220 | 0.228 |
| Recent-parent ownership gap (percentage points) | 12.761 | 8.781 | 8.482 |

The floor weighted losses are 285.640 and 302.722. Their difference summarizes
all twelve weighted rows; it is not evidence of a calibrated winner.
Both completed-fertility normalizations pass the unchanged
\(5\times10^{-4}\) tolerance; rounding the table does not make their exact
normalization gaps zero.

## Every retained floor parameter and restriction

Both selected points are `initial_0001`, the inherited common structural
seed. All eight free coordinates are identical; the child-benefit coefficient
is separately derived to normalize fertility. These are selected saved
values, not estimates from a completed search. The two fertility-shock
parameters are flagged near the lower bound using 1% of the full physical
interval; this flag does not establish a bound optimum. Complete original
parameter tables are [linear](floor_linear/export/parameters.csv) and
[concave](floor_concave/export/parameters.csv).

| Parameter (saved identifier) | Linear | Concave | Lower | Upper | Near bound | Status |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| `H0` | 8.646 | 8.646 | 0.200 | 80.000 | No | Free coordinate; retained selected point |
| `beta_annual` | 0.962 | 0.962 | 0.940 | 0.990 | No | Free coordinate; retained selected point |
| `chi` | 1.127 | 1.127 | 0.100 | 5.000 | No | Free coordinate; retained selected point |
| `first_birth_fixed_cost` | 0.507 | 0.507 | 0.000 | 8.000 | No | Free coordinate; retained selected point |
| `kappa_fert` | 0.209 | 0.209 | 0.020 | 50.000 | Yes | Free coordinate; retained selected point |
| `kappa_fert_continuation` | 0.482 | 0.482 | 0.020 | 50.000 | Yes | Free coordinate; retained selected point |
| `theta0` | 0.088 | 0.088 | 0.000 | 8.000 | No | Free coordinate; retained selected point |
| `h_P` | 1.890 | 1.890 | 0.100 | 2.300 | No | Free coordinate; retained selected point |
| `theta1` | 0.008 | 0.008 | — | — | — | externally fixed B15 |
| `psi_child` | 0.134 | 0.144 | — | — | — | normalized to 2.1 |
| `payroll_tax` | 0.080 | 0.080 | — | — | — | derived from adopted pension ratio and baseline demographics |
| `pension_period` | 0.918 | 0.918 | — | — | — | endogenous balanced PAYGO |
| `housing_supply_elasticity` | 0.630 | 0.630 | — | — | — | retained external setting |
| `tenure_choice_kappa` | 0.005 | 0.005 | — | — | — | retained external setting |
| `alpha_cons` | 0.733 | 0.733 | — | — | — | retained external setting |
| `sigma` | 2.000 | 2.000 | — | — | — | retained external setting |
| `selling_cost` | 0.060 | 0.060 | — | — | — | retained external setting |
| `financed_share` | 0.800 | 0.800 | — | — | — | retained external setting |
| `annual_depreciation` | 0.014 | 0.014 | — | — | — | author adopted input |
| `period_depreciation` | 0.055 | 0.055 | — | — | — | four-year compounded |
| `annual_property_tax` | 0.011 | 0.011 | — | — | — | author adopted input |
| `period_property_tax` | 0.042 | 0.042 | — | — | — | four-year linear source convention |
| `income_process` | 15.000 | 15.000 | — | — | — | retained B15 persistent-state count |
| `entrant_conversion_factor` | 0.500 | 0.500 | — | — | — | legacy child-departure diagnostic; inactive in split-birth entry |
| `adult_entry_birth_to_household_conversion` | 0.476 | 0.476 | — | — | — | effective closed stationary birth conversion |
| `child_benefit_exponent` | 1.000 | 0.860 | — | — | — | fixed sensitivity; not estimated |
| `utility_reference_rent` | 0.110 | 0.110 | — | — | — | fixed experimental common normalization |
| `pension_to_gross_worker_earnings` | 0.229 | 0.229 | — | — | — | adopted CPS2007 target |

The legacy entrant-conversion diagnostic is inactive in split-birth entry;
the effective birth-to-household conversion is the separately recorded
\(1/2.1\). Twelve weighted rows versus eight floor or nine share free
coordinates is a count check, not a rank or identification result.

## What actually ran

The table below reports the last verified initial-bank counts at 02:14 EDT.
Attempted/success/rejection/failure/timeout columns exclude the separately
reused smoke seed. Both floor arms are terminal. Share initial trials had
finished and their two selected repetitions were active, but final share
receipts and Slurm states are unknown after access was lost.

| Arm | New initial attempts | Success | Named inadmissible | Raw failed | Raw timed out | Reused seed | Initial unrun | DE unrun |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Floor, linear | 10 | 2 | 7 | 1 | 0 | 1 | 29 | 120 |
| Floor, concave | 10 | 2 | 7 | 1 | 0 | 1 | 29 | 120 |
| Shares, linear | 39 | 19 | 17 | 1 | 2 | 1 | 0 | 120 |
| Shares, concave | 39 | 20 | 17 | 2 | 0 | 1 | 0 | 120 |

All four arms passed both exact sequential smokes. Each floor original also
matches **each** of its two final repetitions at zero absolute and relative
tolerance. Each share arm had two repetitions active at the last observation;
their outcomes are not counted as completed. The controller summary had not
yet refreshed its repeat counters, so the live heartbeat is the evidence for
those active repetitions. There were no incomplete initial trials at that
observation. Missing final share evidence is an unresolved collection state,
not zero failures or proof of repeatability.

The original maximum was 652 new objective attempts. It was a ceiling, not a
completed count. There are 538 unrun search slots in the last verified state,
including all 480 DE slots. Floor jobs exited with code 2 because the search
was incomplete, even though their repeated selected-point reports succeeded.
Raw state is retained in the [02:14 snapshot](status_20260926_0614.json), each
floor `complete.json`, and the [access receipt](collection_access.json).

## Exact stopping reasons and remaining uncertainty

Both floor arms stopped after `initial_0005` raised `InfeasibleThetaError` at
age 34. Dead-state population mass was about
\(1.194\times10^{-12}\), above the unchanged \(10^{-12}\) gate. The census
contains renters with zero liquid wealth, a child at home and negative budget
slack. Small mass does not excuse the violation. This failure can arise inside
an intermediate price or fertility-normalization solve; it does not establish
that the structural parameter vector is globally infeasible, or exclude a
model bug. The failed evaluations have no accepted equilibrium or scored
loss. Original records and censuses are retained under each floor arm.

Both share `initial_0012` objectives reached the 3,100-second cap. The saved
tracebacks identify `TimeoutError: objective wall budget exhausted`, wrapped
by Numba as `SystemError`. Raw records therefore say `failed`; the
[postmortem](share_deadline_postmortem.json) records the deadline cause
separately. Ten and eleven stationary evaluations had completed before the
interrupted evaluation. This is not evidence of a separate kernel defect.
Later raw share counts include two explicit linear-arm timeouts and a second
concave-arm `SystemError`; that second exception's originating cause has not
been verified from its traceback. It must not be relabeled from proximity to
the deadline alone.

SSH authentication failed on the normal configured connection beginning at
02:21 EDT. The local control socket was absent. No credential workaround,
source edit, retry, gate relaxation, clock extension or new production run
was performed. The cluster controllers do not depend on the laptop, but their
unobserved final state cannot be certified. Final share tables, both
original-versus-repeat comparisons, all seventeen figures, every report page
and terminal counts remain outstanding.

## Verified evidence and visual review

Each actual floor report contains thirteen targets, twenty-eight parameters,
seventeen standard figures and twenty-two pages:
[linear report](floor_linear/export/utility_floor_linear_review.pdf) and
[concave report](floor_concave/export/utility_floor_concave_review.pdf).
The earlier historical layout fixture is not used as actual-results QA.
All floor figures and pages were reviewed. No new report/table clipping was
found. Inherited dense legends overlap curves on pages 14–17, income-state
labels overlap on page 18, and policy legends cross tick labels on pages 19
and 21. These limitations are recorded in the original
[linear visual receipt](floor_linear/export/visual_review_receipt.json) and
[concave visual receipt](floor_concave/export/visual_review_receipt.json).

The original selected case was independently compared with each repetition:
[linear verification](floor_linear/export/repetition_verification.json),
[concave verification](floor_concave/export/repetition_verification.json).
The scientific receipts preserve the passed market, fiscal, demographic,
household-budget, purchase-accounting, occupied-state value and probability
checks. Exact repetition is not grid convergence and does not validate
underwater-owner transitions. Final local readout checks verified all target
row arithmetic, the summed losses, the twenty-eight unique parameters, PDF
hashes and retained metadata hashes; no model was imported or solved.

The immutable launch contract is
`c3fa4d5b7925a54a747c72511538b8182029e1493e4d02e5ef703a30fcecc28a`;
scientific source commit is `d5dbf04d68ff000e1a1b8e66994cdffd31a01580`;
source inventory is
`237904131d159f775c7ae89d1bbf1e8d1dd70c79ad012c80a36d948658d6f9c6`.
The common target/weight/measurement fingerprint is
`8fad0155053df30fc0ccf968c733fcd86e916278ad37cb4b1be62462274557d6`.
The full [target provenance](../launch_v1/target_provenance.json) and arm
receipts retain their complete objective fingerprints and restrictions.

## Interpretation and next decision

The fixed-reference-rent normalization and curvature 0.860 remain
experimental. The former makes material expenditure compensation follow the
equivalence scale at one fixed benchmark rent; childless utility is unchanged.
The pension-to-gross-worker-earnings ratio 0.229 is author-adopted. Earnings,
entry wealth and income-rank coupling, timing, existing transfers, mortality
and targets inherit the frozen reference. No pending estate or mortality
proposal was inserted into this experiment.

The first-birth housing observer remains an unmatched stationary proxy for
the PSID panel estimator. The model's all-positive-estates observer remains
unmatched to the child-directed SCF target. These distinctions limit fit
interpretation. The separate borrowing review leaves underwater-owner
transition validation/fix for later; this stationary experiment cannot close
that issue. The main research task combines estate and mortality evidence.

The optional saved birth-versus-wait value-gap diagnostic is unavailable:
separate action values and eligibility masks were not retained. The bounded
[availability review](supplemental_birth_wait/availability.json) used the
saved-source schema; no new solve or probability inversion was performed and
no zero-shock equilibrium is claimed.

No utility adoption is supported by this incomplete comparison. The immediate
missing work is retrieval and verification of existing share outputs after
normal Torch access returns. A separate
[eight-hour recovery option](../../../../docs/model/utility_four_arm_recovery_v1.md)
allows at most forty new points per arm with time reserved for verification
and reporting. It is a review proposal, not a launched or authorized run.
The earlier 42-hour calculation is only an unapproved cost illustration.
