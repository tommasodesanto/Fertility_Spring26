## Blockers before the no-shock smoke

1. **The source-identity gate omits executable transition dependencies.** `SOURCE_NAMES` does not pin the native perfect-foresight driver or the estate/DUE audit, although this engine calls them through `pf` and `estate` during each mapping. A source drift could therefore pass `check_sources()`.

   - [run_e5f_preference_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_preference_transition.py:29)
   - [run_e5f_preference_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_preference_transition.py:293)
   - [run_e5f_preference_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_preference_transition.py:328)

   Add the direct native driver and audit modules (or replace this list with one complete, verified executable-source manifest) before launch.

2. **There is no test of the new engine’s actual no-shock path.** The only new four-shock test labels itself “Synthetic,” uses a one-period mocked evaluator, and never imports or exercises `run_e5f_preference_transition.py`. The older native PF tests exercise the legacy 2023/static-elastic path, not this driver’s 2007 override, fixed stock, endpoint scaling, split queues, DUE audit, cache, plan gate, or checkpoint/progress behavior.

   - [test_e5f_four_shock_acceleration.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/test_e5f_four_shock_acceleration.py:1)
   - [test_e5f_four_shock_acceleration.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/test_e5f_four_shock_acceleration.py:14)
   - [test_run_e5f_perfect_foresight_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/test_run_e5f_perfect_foresight_transition.py:149)

   Add a Torch-only, pinned no-shock integration smoke: all \(\psi_t\) equal saved \(\psi\), fixed physical stock, saved endpoint, both queues, and cached versus uncached equivalence. Do not describe the current synthetic acceleration tests as complete engine coverage.

3. **The shocked-run disable is self-authorizing.** `--execute`, a matching hash of the supplied plan, and `execution_enabled: true` are sufficient. The plan itself can set that Boolean and provide levels/endpoint; there is no separately pinned approval receipt and no current “no-shock only” condition. The unresolved draft fields prevent an accidental run today, but not an unauthorized shocked run once populated.

   - [run_e5f_preference_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_preference_transition.py:176)
   - [run_e5f_preference_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_preference_transition.py:461)

   Require either a separately pinned, explicit shocked-transition authorization or enforce equality to saved \(\psi\) for the present no-shock-only smoke.

The reviewed timing is otherwise consistent: all four announced levels enter the initial backward path, rows are explicitly relabeled from 2007, fixed stock is constant in the market residual, and both adjusted and raw entry queues are retained and terminal-checked.