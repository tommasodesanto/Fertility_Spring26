# Small positive-credit diagnostic, Phase A and closed GE interface

This isolated experiment retains the **2007 stationary reference — block0506,
September 28 verified export** inputs. The only economic change is the
corrected renter rule with constant unsecured limit (D>0), currently selected
at the smallest one-cent mesh value satisfying the native entrant budget bound.
It is an experimental credit choice, not an adopted baseline or recalibration.
Entry wealth/income pairing, every occupied entrant, earnings, preferences
(including \(\psi_{child}=0.1355551166583114\)), mortgage \(\phi\), fiscal
inputs, owner raw-sale repayment, nonnegative death estates, and physical
housing supply remain pinned.

The bound checks all positive cells of the conditional entrant support and
income grid at \(q/q_{ref}\in\{0.85,1,1.15\}\). The native non-wedge exhaustive
saving objective requires positive surplus after \(c_{bar}+r h_{bar}\); its
upper candidate uses a \(10^{-6}\) search buffer. `c_min=.04` and the extra
`0.01` housing in `kernels.py` are reporting fallbacks, not economic minimum
spending in the active exact-allocation branch. Authenticated local input
inspection gives the unrounded entrant requirement
\(D_{entry}=0.13404978853236008\). The operational one-cent rule
\(\lceil(D_{entry}+10^{-6})/0.01\rceil 0.01\) selects \(D=0.14\), with
0.0059502 units of slack at the limiting cell. This is **only a necessary
entrant condition**, not a lifetime repayment certificate. A native Bellman
and KFE solve at the selected cap must pass, including the inherited entrant
dead-mass gate; a failure stops without automatically increasing the cap.
The fixed-price robustness domain is a diagnostic price guard, not proof of
GE existence outside that interval.

`driver.py smoke` authenticates the pinned bundle and traverses the Phase A
loop with a mock lifecycle call; it must report zero lifecycle evaluations.
`driver.py full` authenticates the frozen observer, runs the selected-cap
native fixed-price case, then calls `phase_b_ge.run_phase_b` under the same
budget. The controller permits at most three Phase A lifecycle evaluations
within 900 seconds and ten total in 2400 seconds. Each case has a 300-second
deadline. The current Phase A design makes only one attempt. Progress,
completed-case summary, deadlines, failure, and pinned input identity are
written during the run. Saved solution arrays are compressed `.npz`; no large
pickle is written. The native fiscal, renewal, target-fit and standard-plot
gates belong to Phase B and must pass before a result is described as closed.

The verified 18-file engine is copied into `source/small_credit_lab/` with no
shared core edits. Only this source and the tiny pinned manifest are staged;
the 120-MB authenticated bundle already exists on Torch at
`/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs`.
`launch_torch.sh` binds the frozen project at the exact Mac root read-only as
in the prior reviewed runtime validation, the staged source read-only, and
results/cache separately. It requests one CPU and 24 GiB for 40 minutes, with
a 2400-second deadline from launcher entry. Stage a source SHA manifest before
the run. The user explicitly authorized immediate `full` launch without the
smoke run; the bundle and frozen observer still authenticate at runtime. Job **18869900** was submitted once and completed successfully in **7m47s**.

Local bundle-only no-solve check: `local_smoke/completed.json` reports the
Phase A mock PASS with \(D=0.14\). The expanded frozen-observer smoke cannot
authenticate the altered local working tree, as documented for the reviewed
refactor; the Torch smoke must use the frozen read-only project bind. The completed production run authenticated inputs and passed the native
Bellman/KFE and closed-GE gates, including the exact selected repeat. The source includes no automatic retry or cap
increase.

## Verified result — September 30

Job 18869900 completed with six lifecycle evaluations, without a production
smoke as explicitly requested. The experimental changes relative to the frozen
reference are removal of the renter age taper, constant unsecured allowance
D=0.14, and full repayment on sale into renting. Entry and all calibrated
preferences, including child benefit psi, remain fixed.

At prescribed reference prices, completed cohort fertility is 2.1061050576
versus reference 2.0999983368. This is a recomputed cohort outcome, not an
immediate response or transition. In the closed stationary GE, prices and
mapped rents rise 0.54127%, household population rises 0.66723%, and completed
fertility is 2.0999980348. Replacement emerges through price adjustment, not
psi normalization. Household population is not a person headcount.

Renewal residual is -9.35796e-7 (absolute tolerance 1e-6); population-scaled
housing residual is 8.88e-16 and actual PAYGO residual 1.548e-13. Entry
censoring, renter floor violations, purchase violations and negative estates
have zero mass. Native queue raw L1 is 5.81554e-8; renewal-adjusted L1 is
2.76263e-14. Finite-grid and estate-creditor accounting caveats remain.

Selected and repeated tables and all 17 standard PNGs match exactly. See
[collection verification](collected_v1/collection_verification.json),
[full 14-row fit table](collected_v1/full/phase_b_ge/selected_root/target_fit.csv),
[full 31-row parameter table](collected_v1/full/phase_b_ge/selected_root/parameters.csv).
Checkpoints and numerical arrays remain on Torch. This diagnostic does not
replace the frozen reference; no transition was computed. The numerical budget
is closed.
