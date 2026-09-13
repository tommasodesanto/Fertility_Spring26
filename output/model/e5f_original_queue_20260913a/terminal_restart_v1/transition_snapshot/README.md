# 100-period transition: third completed mapping

Collected September 13, 2026 at 19:23 UTC from job 17699174. This is a
snapshot of an unconverged 100-period (400-year) path with a verified stationary
endpoint. All 100 dated house prices match evaluation 3 in latest_completed.json.
The full root score is 8.48344 against 0.0002 required. Maximum housing, pension,
and rebate relative imbalances are 7.4765%, 0.7082%, and 4.2417% respectively.

The final saved row is the period starting at year 396. Household mass there is
0.388248 versus terminal mass 0.349454 (11.10% above). The final carried full
distribution and entry-queue gaps are not persisted until this root exits;
these population comparisons do not substitute for those checks.

The six-panel graph uses native saved quantities. The birth panel displays
birth counts per four-year period, not period TFR: the running driver only
persists its age-standardized fertility observations on root exit. The same
initial and terminal equilibria as the ten-period experiment supply references.

Regenerate without a model solve:

```sh
python code/model/tools/build_e5f_stationary_shock_figures.py --case-dir output/model/e5f_original_queue_20260913a/terminal_restart_v1/transition_snapshot
```

rent_identity_check.json verifies every rent in both the 10- and 100-period
snapshots against carrying costs plus expected capital losses. Errors are
below 6e-17. In the short path the terminal price jump contributes 0.168014
of final rent 0.232492; in the long path it contributes 0.010411 of 0.048772.
Both contain boundary effects, and neither is an accepted equilibrium path.
The native economic solver and active cluster jobs were not changed.
