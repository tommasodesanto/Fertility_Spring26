Established chronology and intent:

- The strict gate was introduced in commit `4f8572f14` (2026-09-02, “Stabilize positive-kappa transition accounting”), not with the September 7 joint-nested experiment. The diff adds a helper explicitly described as “exact mass for a pure redistribution after a fail-closed gate.” It rejects nonfinite/nonpositive mass or a relative discrepancy above \(5\times10^{-9}\), then rescales only a passing distribution. [run_e5f_transition_calibration.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_transition_calibration.py:572)

- That same commit inserted the gate after advancing each matched treated/control first-birth branch, and subsequently after destination feasibility, fertility, and current-choice transformations. The intent is accounting preservation for a matched one-period DID cohort—not permissive numerical repair. [begin branch](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_transition_calibration.py:1006), [advancement gate](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_transition_calibration.py:1085), [destination gates](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_transition_calibration.py:1154)

- The pre-existing “small mass pruning” is separate and older: the active optimized Markov transition skips location-origin, tenure, and tenure-destination blocks below \(10^{-15}\). Git blame dates these lines to the optimized-port commit `e31e2f9a` (2026-07-19); history traces the same rule back through the July 18 seq fork and June 4 Markov-income implementation. [solver.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:5209)

- This creates the current mechanism hypothesis: an entire subcohort of mass \(\sim3.55\times10^{-9}\) is not itself below \(10^{-15}\), but its post-location/tenure fragments can be. Those fragments are discarded, so the later branch-level relative gate correctly sees a \(1.745\times10^{-5}\) discrepancy. That causal link is an inference from the code, not yet an instrumented proof.

Current primary evidence:

- Job `17093420` is failed at initial population after four completed histories; task 018 reports
  actual `3.5540670288250579e-09`, expected `3.5541290338151895e-09`, relative gap `1.745e-05`, tolerance `5e-09`. [search_state.json](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_full_20260906a/wide32_support_g/search/search_state.json:7)
  This is distinct from the earlier canary timeout and replay-26’s zero-support rejection.

- The current status/memory predates that 06:09 UTC failure: it states no numerical gate was relaxed and documents the earlier 2015 zero-support rejection as a typed candidate rejection, not a calibration result. [CALIBRATION_STATUS.md](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/CALIBRATION_STATUS.md:1), [daily note](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/memory/daily/2026-09-07.md:55)

Scale-invariance finding:

- The gate’s rescaling \(g\mapsto g\,M/\sum g\) is scale-invariant for a nonzero branch that passes it: conditional distributions and conditional housing means are unchanged. But it has no pre-existing protection against the absolute \(10^{-15}\) pruning threshold in the Markov operator; that threshold is inherently not scale-invariant for tiny subcohorts.
- I found no historical note establishing an authorized solution (e.g., relative pruning or branch renormalization before advancement). The available history establishes fail-closed intent, not permission to relax it.
Lead clarification: the four completed histories in the failed search state are imported smoke receipts. Zero new full search cases completed. The original September2 gate introduction was independently checked with git show4f8572f14.
