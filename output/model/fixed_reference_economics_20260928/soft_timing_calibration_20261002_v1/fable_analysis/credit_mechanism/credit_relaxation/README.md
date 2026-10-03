# Fixed-price credit relaxation at the revised-timing chain 13 point

This is an **experimental partial-equilibrium credit counterfactual**, not a
recalibration, general-equilibrium result, or dated transition. The only changed
economic primitive is the uniform financed share \(\phi\), from 0.8 to 1.0.
The price \(p=0.7760569760205563\), income and entrant law, payroll tax,
housing supply coefficient, preferences, and all other loaded model fields are
held fixed. The production engine retains its redundant purchase screen.
`prepared.json` records every production Python source hash, input snapshot
hashes, the pinned chain-13 array hash, the planned two sequential one-core
solves, and the 600-second per-solve and 1500-second outer deadlines.

The [matched native-node overlay](matched_native_policy_overlay.png) compares
first-birth attempts and housing/tenure policy schedules for ages 22–25,
inherited renters, income states 4 and 6, at identical beginning net financial
wealth \(b\). Solid lines condition on no children; dashed lines condition on
one child at home. These family-state schedules are not forced-birth transitions.
The financial-credit experiment itself is a causal change **within the fixed-
price model**. `overlay_data.json` holds exact weighted means. We use the same
\(\phi=0.8\) pre-fertility childless native-node weights in both arms, so those
means compare policies at fixed state composition. The two solved distributions
are saved separately and are not mixed into those means.

| Income state, fixed baseline weights | \(\phi=0.8\) | \(\phi=1.0\) |
|---|---:|---:|
| 4: childless ownership probability | 0.0377 | 0.2951 |
| 4: one-child ownership probability | 0.000494 | 0.000507 |
| 4: childless first-birth attempt probability | 0.005615 | 0.005439 |
| 6: childless ownership probability | 0.4366 | 0.4393 |
| 6: one-child ownership probability | 0.4345 | 0.4381 |
| 6: childless first-birth attempt probability | 0.737862 | 0.738171 |

`owner_room_decomposition.json` shows why state 4's large ownership response
does not carry through to one-child families. The two-room owner probability
rises from 0.033672 to 0.289265, accounting for 99.3% of the total childless
ownership gain. The two-room rung is unavailable in the one-child state because
the executed parent housing floor is 2.593760 rooms. One-child families can
still rent up to six rooms. This is evidence for a **rung mismatch** at these
states, not evidence that every larger-family housing need forces ownership.

The baseline \(\phi=0.8\) fixed-price solve was compared against the pinned
cluster `selected_repeat` arrays with a strict global absolute tolerance of
\(10^{-10}\). That gate failed on numerical differences, chiefly \(1.91\times
10^{-6}\) in an infeasible large-negative value; it was not relaxed. The saved
price bits matched exactly. `historical_pin_discrepancy.json` records every
core-array maximum and location, and `failure.json` preserves the stopped
attempt. A subsequent authenticated same-host reference, already verified
against the deployed production engine, matched the fresh baseline exactly on
all ten core arrays. `samehost_control.json` pins that reference, its verification
receipts, and the fresh baseline hashes. The continuation then ran **only**
\(\phi=1.0\); no baseline rerun was needed.

`completed.json` records both solve times (about five seconds each), the exact
same-host control, executed primitive input difference of only `phi`, mass
checks and propagation to `shared.phi_choice` and `shared.phi_state`. Each arm
has saved arrays, executed input snapshot, and the established 17 standard
diagnostic figures under `phi_080/` or `phi_100/`. Housing-market figures in
those folders are fixed-price diagnostics and do not claim market clearing.

Regenerate the overlay without solving:

```sh
output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python \
  output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/credit_mechanism/credit_relaxation/plot_overlay.py
```
