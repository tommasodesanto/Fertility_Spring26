Assessment: the original 64 points cover \((\kappa,\lambda)\) broadly, including 18 low-\(\kappa\) cells, but do **not** cover the hypothesized low-\(\kappa\), *small absolute* 2007–23 child-preference decline region. The saved stopped run has no completed initial-population cases, so it does not establish that low \(\kappa\) was actually eliminated; it does establish that the controller would not preserve that region if its proposals were rejected.

1. Physical coverage

- Domain: \(\kappa\in[0.005,10]\), \(\lambda\in[0.02,1]\), and \(\Delta\psi_{2007\to2023}\in[-1.5,0.2]\), with log, log, and asinh transforms respectively ([adapter, line 36](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_overnight_case.py:36)).
- Population construction is one exact anchor, 47 explicit grid points, and 16 seed-near perturbations. The grid is
  \[
  \kappa\in\{.01,.03,.1,.3,1,2,4,8\},\quad
  \lambda\in\{.02,.05,.2,.5,.8,1\},
  \]
  omitting \((2,.8)\) because the anchor occupies it ([search controller, lines 131–156](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_nested_long_search.py:131)).
- Thus 18 points have \(\kappa\le .1\), each paired with every listed \(\lambda\). The low-scale corner \((.01,.02)\) is explicitly case 2.
- However, the realized preference-change support of all 64 saved centers is only
  \[
  \Delta\psi\in[-0.5092,-0.1578].
  \]
  For \(\kappa\le .1\), it is \([-0.4364,-0.2231]\); in particular case 2 is \((\kappa,\lambda,\Delta\psi)=(.01,.02,-.2452)\). There is no low-\(\kappa\) point with a decline near zero or, say, \(-.05\) to \(-.10\). This is a real joint-coverage gap, not a transform artifact: the repaired profile defines the coordinate as the change from the normalized 2007 intercept ([transition driver, lines 109–114 and 1880–1893](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_transition_calibration.py:109)).

2. What the evidence establishes

- The valid anchor establishes one complete, unoptimized feasible history at \((2,.8,-.3287)\), with loss \(485.2815\), normalized old completed fertility \(2.1000255\), and no evidence about low-\(\kappa\) fit. Its parameter table confirms \(\kappa=2,\lambda=.8\) are free coordinates ([anchor summary](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_full_20260906a/support_repair_e/smoke/smoke_anchor/task_001/summary.json)).
- The repaired replay is **not** a low-scale test: case 26 has \(\kappa=1,\lambda=.02,\Delta\psi=-.3766\). It normalized old fertility to \(2.10000344\), then rejected at the 2015-origin first-birth branch. It neither reached nor identifies the 2019–23 target ([replay assessment](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_full_20260906a/support_repair_e/replay_assessment.json)).
- The original search stopped on case 26’s then-unclassified stationary exception; its ledger contains only the two anchors and two all-coordinate smoke probes, not any completed initial member. Therefore no evidence shows whether the 18 low-\(\kappa\) proposals would have survived ([stopped ledger](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_full_20260906a/wide32/failed_search/all_cases.csv)).

3. Can DE recover discarded low scale?

Not reliably. The repaired controller can now classify a dated zero-support failure and continue ([lines 89–100](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_nested_long_search.py:89)); but after the initial stage every rejected slot is replaced by a cyclic draw from the globally ranked valid survivors, not by a nearby or same-\(\kappa\) survivor ([lines 495–505](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_nested_long_search.py:495)). There is no diversity floor or low-scale re-injection.

DE’s affine difference mutation can algebraically move outside the surviving coordinate range before clipping, so recovery is possible in principle ([lines 167–174](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_nested_long_search.py:167)). But with at most eight generations and selection only on complete losses, it provides no deterministic recovery of a low-\(\kappa\), small-\(|\Delta\psi|\) region once all its explicit initial members are rejected. The controller therefore risks eliminating that *joint* region; it does not prove that low \(\kappa\) itself is infeasible or cannot fit.

4. Conditional small amendment

If the pending small-scale canary and full smoke clear, a minimal unchanged-budget design would **replace, not add**, six existing low-\(\kappa\) grid rows:
\[
\kappa\in\{.01,.03\},\quad
\lambda\in\{.02,.2,.8\},\quad
\Delta\psi=-.05.
\]
This keeps 64 histories, the eleven-dimensional domain, targets, weights, and gates unchanged, while testing the missing joint direction. It would require changing only the explicit initial-population construction, then regenerating the centers and immutable contract because the controller and contract hash those artifacts ([builder, lines 52–70](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/build_e5f_joint_nested_long_contract.py:52)). Required checks: exact physical-coordinate table from regenerated centers; domain/transform validation; unchanged-anchor reproduction; completed smoke and canary review; and confirmation that any rejection remains recorded without synthetic loss.

I would not adopt that amendment yet: the evidence is sufficient to identify a joint-coverage deficiency, but not to infer that low \(\kappa\) needs a smaller decline or that its existing grid points will fail.