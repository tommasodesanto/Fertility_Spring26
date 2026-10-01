# Fixed-price expenditure-share and A(m) comparison

Experimental author-authorized matched three-arm comparison; no adoption or recalibration. The reference is the provisional with-A `nonnegative_mean_120x9_s1/037_r2_gn_1.0` point, whose exact Torch parameter/target hashes match `utility_calibration_round1_v1/deployment/attempt2/comparison_requested/with_A_{parameters,target_fit}.csv`. The authenticated reference price is **0.8198539089139171**, from `reference/with_A_closure.json`; reference psi is **0.1355551166583114**, share jump **0.1303383207216736**, and h_P is zero.

All arms keep the old shared parameter vector, benefit curvature and equivalence scale, earnings, nonnegative mean-preserving five-node entrant law, 120×9 grid, corrected zero unsecured borrowing limit, 2% annual rate, H0, fiscal objects, housing menus, and original 14-target/weight contract. Economic changes are experimental: arm 1 retains the old variable first-child shares and A(m); arm 2 turns off A(m) with those same shares; arm 3 also turns off the first-child share jump, holding shares at their childless level. All three have h_P=0. These are the existing native utility switches, with no model-source change.

Comparison 1 versus 2 isolates A(m) conditional on the old variable shares. Comparison 2 versus 3 isolates those shares with A(m) disabled. This is not a full factorial or unique additive attribution: A(m)=1 at constant shares and the share/compensation interaction remains relevant.

The price and psi are fixed. No fertility normalization or GE root is used. Actual fertility/TFR, birth renewal and housing-market residuals are reported; equilibrium and demographic renewal are not imposed. Household, purchase accounting, no entrant projection, PAYGO and metadata gates are retained. Each arm must provide full 14-target fit and all 31 parameter/restriction rows, the stable 17 diagnostic PNGs, actual shared utility arrays and fingerprints, and supplemental parent housing/rent profiles.

One Torch job is bounded to one CPU, 24 GiB, one thread, 1200 seconds total, three lifecycle calls, and 600 seconds per cell including reporting. The exact zero-lifecycle initializer and mocked three-cell controller run before the native loop. A failure stops the packet, with no retry or broader numerical search. Source packaging reuses the reviewed round-three and repaired prescribed-price runtime inventories, including the full small_credit_lab package. The first table is saved after every completed cell.

Lead reviewed the driver and source bindings. The exact Torch zero-lifecycle initializer validated all three actual 31-parameter mappings and utility arrays; the mocked three-cell controller passed. All 212 frozen source pins were checked again before submission. Job **18927323** was submitted October 1 at **01:10:30 New York** and completed all three cells with three lifecycle calls total, no fatal failure, and no retry. All arms are prescribed-price comparisons; their closures report that neither the market nor renewal root is imposed. Each arm has 14 target rows, 31 parameter/restriction rows, and 17 standard diagnostic PNGs locally under `collected/{with_A,no_A,constant_alpha}/`. Arrays are excluded from this compact backup. The compact completion receipt is `collected/monitor_status_20261001T0123NY.json`. The two earlier zero-lifecycle initializer failures are retained under deployment/initializer_attempt1 and initializer_attempt2. The repaired passed receipts are in deployment/initializer_passed, and deployment/submission_receipt.json records the numerical budget. No existing calibration/search, credit, transition or manuscript artifact is modified.

Reference Torch receipt: `/scratch/td2248/projects/entry_calibration_multistart_v1/results/nonnegative_mean_120x9_s1/run/037_r2_gn_1.0/phase_b_ge/selected_root/`.


Completed prescribed-price outputs:

| Arm | Original-weight loss | Mean rooms | Ownership 30–55 | First-birth rooms | Early fertility | Renewal residual | Market residual |
|---|---:|---:|---:|---:|---:|---:|---:|
| with_A | 18.128926 | 5.754419 | 0.663621 | 1.607464 | 0.524656 | 1.12e-7 | 0.038821 |
| no_A | 4123.190441 | 5.943430 | 0.665337 | 1.315158 | 0.983731 | 0.402624 | 0.007250 |
| constant_alpha | 1082.983025 | 4.964653 | 0.621678 | 0.147705 | 0.629279 | 0.123928 | 0.170738 |

Each arm’s `target_fit.csv`, `parameters.csv`, `closure.json`, and 17 standard PNGs are retained under `collected/<arm>/`. These outputs are prescribed-price results, not market- or renewal-clearing roots.
