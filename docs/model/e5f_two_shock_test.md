# Two contemporaneous logistic differences: isolated computational test

Author authorized 7 September 2026. Branch `codex/two-contemporaneous-shocks`,
based on `aa2c71de`. Production and the prior nested-GEV experiment are preserved.

For tenure d and attempt a in {0,1}, the utility of a plan is
Q[d,a]+d*H+a*F. H and F are independent centered logistic differences. Housing
scale is tenure_choice_kappa; fertility uses kappa_fert for the first birth and
kappa_fert_continuation (falling back to kappa_fert) subsequently. No scale ordering.

The four deterministic plan values retain the prior joint experiment's menu:
tenure is committed across conception success/failure, while product size and
continuous decisions adjust within tenure afterward. Individual owner-product
shocks are absent. This test does not isolate shock timing from these retained
menu differences relative to production.

The operator uses an exact housing envelope plus a bounded fertility softplus
residual. Split quadrature resolves tenure switches and narrow fertility layers;
its embedded error estimate is checked against independent SciPy quadrature,
refinement, symmetry and the probability/value derivative identity. Only the
omitted-tail error has an analytical bound. A compiled explicit stack avoids
unsafe caching of recursive compiled calls observed in the first local attempt.

Scope: 16 feasibility masks under three scale pairs, separable limit, extreme
scale ratios, translation and derivatives, Bellman wiring, existing joint
population-accounting tests; then one exact saved-continuation reproduction and
two full lifecycle solves at inherited prices. The wealth and demographic grids,
parameters and numerical gates are unchanged. This is not market clearing or
calibration and does not produce an SMM loss.

Run budget: 1,200 seconds for the three-case lifecycle probe, single local core,
30-second heartbeat, fail on numerical/accounting errors. Prior cluster full
smoke took about 55 seconds; new operator benchmark was 1.66 seconds per 10,000
random four-plan states (maximum probability sum error 4.44e-16). The expected
local run is minutes, bounded below the 30-minute cluster-only threshold.
Cluster login was rejected during this session; no new cluster job submitted.

Python 3.12 with NumPy 2.2.6, SciPy 1.15.3 and Numba 0.61.2 is the isolated
local probe runtime. The old project runtime passed operator tests but cannot
read the NumPy-2 checkpoint. No production runtime was modified.

Reproduction from this branch:

```sh
PYTHONPATH=code/model:code/model/tools python -m unittest discover -s code/model/tools -p 'test_e5f_two_shock_choice.py' -v
PYTHONPATH=code/model:code/model/tools python -m unittest discover -s code/model/tools -p 'test_e5f_joint_nested_full.py' -v
PYTHONPATH=code/model:code/model/tools python -m unittest discover -s code/model/tools -p 'test_e5f_joint_nested_integration.py' -v
python code/model/tools/run_e5f_two_shock_probe.py --checkpoint /absolute/path/to/dated_state.pkl --output output/model/e5f_two_shock_test_20260907a --seconds 1200
```

The authoritative checkpoint is root output/model/e5f_overnight_independent_verification_20260905a/numerical_smoke/dated_state.pkl,
SHA256 bbe10a21a843facaf2bceed56e89281e00992d632cddab640ff0c572d3eb494f.
The probe writes source hashes, fixed prices, dimensions, scales, completed-case
summaries and policy arrays. Figures are suppressed under the author's explicit
request not to create illustrations without asking; this overrides the usual
standard diagnostic graph requirement for this test.

Independent Astra/max review confirmed the operator and mass-factorization
reuse; its flagged continuation-scale fallback, incompatible flags, nonfinite
result checks and thin-layer resolution were incorporated. It did not certify
numerical results or calibration.

## Local reference discipline

The pristine parent source differs from the saved cluster checkpoint only at
floating-point roundoff (maximum 3.553e-15 in V). It is therefore recomputed
in the same isolated runtime using `tools/check_e5f_two_shock_local_reference.py`.
The new code with both experimental flags off must reproduce all ten of those
pristine-parent arrays bit for bit. This check passed. Cross-platform differences
remain reported separately; no tolerance was relaxed to hide a code difference.
Pass its `reference_arrays.npz` with `--local-reference` to the probe.

The first full experimental call exposed an unnecessary dependency on the old
GEV lambda field in a production checkpoint. The additive branch no longer reads
that field; its Bellman wiring test deliberately omits it. All sixteen operator,
Bellman and mass-accounting tests pass in the final isolated runtime.
