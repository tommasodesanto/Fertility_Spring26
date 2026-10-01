# Housing decomposition of verified local chain 11 postcheck

Original-weight loss 114.13604148137114. Experimental candidate; not adopted. q=0.689763713921058.
Saved g is the realized post-housing distribution. Renter housing uses only g[:,0] and hR[:,0]; owners use physical room rungs [2,4,6,8,10]. Current children m=min(n,cs,3) under pinned independent_count mode. Age intervals use uniform four-year-cell overlap.

| Group | Mass share | Ownership | Renter rooms | Owner rooms | Renter cap % | Aggregate room contribution |
|---|---:|---:|---:|---:|---:|---:|
|All ages 18-85|100.00%|59.29%|4.4494|7.3096|28.91%|6.1451|
|All ages; current children 0|61.56%|57.33%|3.8825|6.3243|16.40%|3.2519|
|All ages; current children 1|28.02%|59.91%|5.4761|8.7152|51.66%|2.0784|
|All ages; current children 2|9.18%|68.57%|5.4932|8.8490|51.76%|0.7158|
|All ages; current children 3|1.23%|73.52%|5.4999|8.9482|50.33%|0.0990|
|Ages 18-29|18.52%|25.09%|4.1015|8.0973|24.87%|0.9453|
|Ages 18-29; never had children|9.76%|8.54%|3.2626|6.4257|4.24%|0.3448|
|Ages 30-45|24.69%|53.90%|4.6388|8.4061|31.49%|1.6469|
|Ages 46-65|30.87%|66.42%|4.6403|7.7198|31.13%|2.0636|
|Ages 66-85|25.92%|80.36%|4.5852|6.0294|29.61%|1.4893|
|Ages 25-34; never had children|4.86%|9.95%|3.2827|6.2325|4.74%|0.1738|

The existing matched first-birth observer gives control 5.975093924961874 rooms and treated 6.955108849086157 rooms, a 0.9800149241242826 room response one four-year period later. These include renters and owners. Compact observer JSON retains branch levels but not cohort distributions; its mean near six does not itself establish a renter cap.

Complete target fit and all 31 parameters remain in /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/chain_11/postcheck_results/selected_postcheck/phase_b_ge/selected_root.

## Matched-cohort decomposition completed

The existing dated first-birth observer was rerun with authenticated saved parameters, shared arrays and policies, with model solves forbidden. A Python function-return profiler captured its realized current distributions; the original observer formula was preserved. The six saved scalar results—control rooms, treated rooms, origin mass, destination mass, response and continuation births—reproduce bit for bit. The independent tenure sum differs from the control observer by one floating-point unit (8.88e-16).

| Destination branch | Renters | Owners | Renter rooms | Owner rooms | Renters at six-room cap | Total rooms |
|---|---:|---:|---:|---:|---:|---:|
| Held-childless control | 53.57% | 46.43% | 4.9365 | 7.1732 | 32.01% | 5.9751 |
| First-birth treated | 49.67% | 50.33% | 5.4770 | 8.4137 | 55.91% | 6.9551 |

The near-six control average includes substantial ownership. It does not mean control renters all hit their cap: their conditional mean is 4.937 and 32.0% bind. Cap binding is more common among treated renters (55.9%). Aggregate young-childless renters above use a different risk set and cannot substitute for this matched comparison.

The first-birth response of 0.980015 rooms equals the change in renter room contribution, +0.075978, plus the change in owner contribution, +0.904037. This is an accounting decomposition combining changes in tenure fractions and conditional room choices; it is not a causal attribution of the mechanism.

These branches share the origin successful-first-birth risk-set weights, including income and wealth, and advance one four-year period with existing policies. Treated households may have another child at the destination; controls are held childless. Origin mass is 0.05117423014844; treated continuation births are 0.01885554564488.

The complete matched table is matched_birth_housing.csv; source pins, exact reproduction, zero solve count and memory evidence are in matched_birth_verification.json. Reproduce using code/model/.venv/bin/python followed by the absolute path of matched_cohort.py. No GE, Bellman, lifecycle or search call is permitted by that script.

Origin-to-destination tenure switches remain unavailable: branch destination marginals do not preserve the required joint tags. Purchase affordability is also not an output of the existing observer; its exact state-cash-rule measurement would be a separate extension. Neither limitation blocks the realized tenure and cap result. The earlier matched_birth_saved_summary.csv records the compact observer's limited original fields and is superseded for tenure/cap values by matched_birth_housing.csv.

This remains local chain 11 at loss114.136041, not the newer verified Torch point at loss109. No calibration point was adopted.
