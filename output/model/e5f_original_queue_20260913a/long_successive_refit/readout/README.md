# Requested transition and policy readout — September 14

`big_run_transition.png` and `.pdf` show the latest completed mapping frozen
for this readout: stage 0, trial 1, round 1, mapping 7. This is a 104-date
forecast after the first candidate 2007 surprise, with preference 0.1289153146
held constant. It is **not the completed successive-surprise fitted history**.
No model solve or slide edit was performed.

The four panels preserve the main transition graph's objects: period fertility,
household heads, aggregate housing services, and house prices. Orange ends at
the last saved decision date, 2419. Red references are the independently verified
steady state for this candidate, placed at the 2423 boundary; the plotted path
is not joined to or forced through those references. The initial stationary
reference uses the same specification and exactly matches initial age masses.
Terminal fertility is stationary adjusted birth flow divided by entrant flow;
the closure has equal masses across fertile ages and no fertile-age mortality.

Largest housing, pension, and rebate gaps are 0.003773%, 0.001464%, and
0.003348%. First-period model fertility is 1.953246 versus target 1.974875,
outside the 0.005 fit tolerance. No accepted stage receipt existed at collection.
Terminal household mass remains 1.742% away on the saved next-state check;
distribution relative L1 distance is 9.440%. Horizon robustness remains open.
The displayed final decision-date household index is 74.484, versus 73.168 at
the stationary endpoint. This timing differs from the next-state distance check.

Regenerate the graph and its verified table without solving:

```sh
/Users/tommasodesanto/miniconda3/bin/python code/model/tools/build_e5f_big_run_readout.py
```

`frozen_big_run.json` contains native arrays, best-coordinate receipts, terminal
verification, and source hashes. `verification.json` records exact coordinate
matching, artist-array checks, birth-flow identities, and selected years.

## Policy comparison

`policy_comparison.csv` and `.json` compare each 2% arm with its own 1% baseline.
Both taxes rebate all receipts equally. The slide comparator is the frozen
one-shock packet in `inherited_2023_tax/transition_readout/frozen_policy_paths.json`.
The four-announced-shock comparator is original mapping 7 of
`announced_original_queue_20260913c/output/run`, whose rows hash matches the
recovered 2023 state receipt exactly. The 2% policy is round 2 mapping 7 of
`fourshock_tax_20260914b`; its best coordinates match native arrays exactly.
These histories are distinct from the unfinished big successive-surprise refit.

At 2063, the one-shock slide has births +1.031010% and household stock +0.113584%;
four announced shocks have births +1.001056% and household stock +0.111223%.
Both start with the same household mass as their respective baseline and have
slightly lower births on impact. The four-shock 2% path remains provisional:
its maximum housing imbalance is approximately 0.2423% and rebate imbalance
0.2417%. Small differences between these policy effects should not be treated
as precisely resolved while numerical convergence remains outstanding.

National illustrations use 128,086,487 householder heads and 334,914,896 people
in 2023, holding headship fixed. Additional household stock in 2063 is the model
policy-minus-baseline mass divided by common initial baseline mass, multiplied
by the national householder count: 125,673 for the slide and 124,455 for four
shocks. Corresponding people equivalents are 328,603 and 325,420. These are
illustrative national scalings, not separately modeled resident-population
forecasts and not cumulative household formations. The 1.805% previously
reported increase refers to the **stationary endpoint**, not 2063. It uses
the full final preference decline, not this first-candidate graph's preference.

`fourshock_policy_source.json` preserves the bounded extraction, matched source
hashes, and 2023–2063 arrays. The national householder input is
`code/data/Spatial_aggregate_withmicrodata/output/national_householder_housing_path/national_householder_housing_path.csv`;
the population input is ACS 2023 table B01003. The geographic normalization
limitations already recorded in the decision ledger remain applicable.
