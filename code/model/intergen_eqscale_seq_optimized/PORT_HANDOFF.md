# E1 optimized-port handoff

This directory is a physical copy of `intergen_eqscale_seq`; the three source
packages were not edited.  The port reimplements the accepted optimizations
from `intergen_housing_fertility_optimized` in the E1 extension architecture.

## Changed implementation objects

- `parameters.py`: direct-Brent and compiled-scatter defaults, guarded override
  registry, and income-transition rebuild/stationarity contracts.
- `solver.py`: price-cache-backed direct scalar search with directional
  bracketing and legacy-damped fallback; accepted Bellman payload reuse;
  safeguarded-Brent degenerate-denominator bisection guard; Fortran
  family-state decoder use; and the gated full-Markov-support
  `add_old_wealth_income_moments` retirement-statistics correction.
- `utils.py`, `state_layout.py`: explicit family-state layout contracts.
- `calibration.py`, `target_system.py`, `run_e1_chain.py`: atomic 15-row E1
  target system (including occupied rooms) and stable target fingerprint.
- `e1_profile.py`, `benchmark_e1.py`, `parity_panel_e1.py`,
  `promotion_worker_e1.py`, `promotion_collect_e1.py`, and
  `cluster/submit_e1_promotion.sh`: package-owned E1 benchmark/promotion
  entry points.
- `tests/test_optimized_port.py`: direct sequential payload-staleness,
  direct-versus-legacy roots, repeat-price cache, explicit default-nesting,
  and contract regression checks.  The inherited nesting test remains as-is:
  it exercises the production-compatible defaults in the original fork;
  this package's new direct/scatter defaults are explicitly switched off in
  the added nesting test to establish the equivalent compatibility envelope.

## Extension-specific payload rule

The retained payload is the reference's eleven objects plus a twelfth,
`fert2_probs`.  The sequential Bellman writes `P._fert2_probs`; the KFE reads
that field for second-birth splits.  `upgrade_fast_markov_solution` restores
the payload's twelfth member immediately before its full KFE.  Thus a rejected
last price cannot contaminate the accepted price's second-birth flows.  Audit:
`_fert2_probs` is the only Bellman-to-KFE `P._*` input; `_first_births_by_age`,
`_second_births_by_age`, and `_second_attempts_by_age` are KFE outputs.

The integrated regression starts direct Brent at `1.30 * E1_PRICE` with
`scalar_market_refine_max_expand=0`, records the fast Bellman prices, and
requires the accepted price to differ from the final evaluated price. Its
integrated second-birth flows and distribution then equal a fresh accepted-price
solve exactly. `p_init_override` is an explicit supported dynamic override.

### Golden-owner family-state correction

The parent fallback used `ceil((c+1)/n_child_states)` and a remainder based on
`n_child_states`; that is not the package's Fortran flattening
`c = parity + n_parity * child_state`. The optimized fork deliberately uses
`decode_flat_family_state(c, n_parity) = (c % n_parity, c // n_parity)` in both
golden-owner fallback sites. Thus it can differ from the parent only when the
Python golden-owner fallback is selected (rather than the compiled owner block)
and there is more than one parity/child state; there it corrects the borrowing
floor selected for a flattened family column. Fixed-price E1 parity is exact
because its active algebra uses the compiled path.

The fused age/KFE kernel was deliberately not ported: the reference measured
it slower.  The retirement-income old-wealth correction is ported but inactive
under E1 because `retirement_income_z_scale=0`; it replaces the income-collapsed
old-age mean/median family only when that scale is nonzero.

## Promotion tasks

`promotion_worker_e1.py` writes one atomic JSON checkpoint for each task.
`demand-map` compares both packages at 17 bracket-spanning fixed prices;
`root-case` runs ordinary and forced-wrong-direction optimized roots;
`strict-repeat` checks two strict roots bitwise; and `throughput` runs a
deterministic ten-candidate beta sweep in both packages. `--smoke` preserves
the task shape with a three-point/three-candidate compact packet.
`promotion_collect_e1.py` adjudicates all four checkpoints. The submit script
uses one CPU and selects the task through `PROMOTION_TASK`; `PROMOTION_TASK=collect`
runs the collector. The throughput packet reports both fixed-price and loose-GE
totals; the latter is the promotion-relevant throughput measure. The forced
root starts outside the bracket with expansion capped and fails if the recorded
bilateral-fallback flag is absent. No job is submitted by these tools.

## Gates and commands

```bash
code/model/.venv/bin/python3 -m py_compile $(find code/model/intergen_eqscale_seq_optimized -name '*.py')
PYTHONPATH=code/model code/model/.venv/bin/python3 -m pytest code/model/intergen_eqscale_seq_optimized/tests/ -q
PYTHONPATH=code/model NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  code/model/.venv/bin/python3 -m intergen_eqscale_seq_optimized.benchmark_e1 --package both --mode fixed-price
```

The measured fixed-price paired speedup is `1.21x`; value/policy/probability
parity is exact, distribution maximum difference is `4.3e-19`, and target-
moment maximum difference is `4.44e-16`. Strict repeats are bit-identical.
Loose-GE throughput remains pending the promotion battery. The port now
contains 93 discovered test functions (the lead's pre-adversarial run reported
89 passing tests); this sandbox's pytest process exits with status 139 before
collecting/reporting tests. No `sbatch` was issued and no E1 artifact was
changed.
