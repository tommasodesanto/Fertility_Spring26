**September 30 full-GE verification passed:** job 18880497 exactly matches the preserved 160×15, D=0.14 control in eight closure receipts, 87 arrays per saved bundle, all fit/parameter rows and 17 plots per final. Evidence: [verification packet](../../../../output/model/publication_refactor_20260929/single_market_verification_v1/README.md). Singleton axes remain; this is phase-1 equivalence, not completed location-axis removal or a speed certificate.

# Stationary single-market specialization — isolated phase 1

Reference: **2007 stationary reference — block0506, September 28 verified export**.
This package is an isolated experiment; it has not replaced `refactor_lab`, the
legacy model, or the transition runtime. No economic parameter or numerical
resolution is changed. Source pins are in `source_provenance.json` and the
minimal source comparison is `review.diff`.

## Supported contract

`contract.validate_contract` rejects unsupported structure before shared
construction, Bellman solution and stationary distribution propagation. The
supported model has one housing market, one Markov income process, sequential
birth choices, and independent counts of children at home. It keeps four
lifetime-child states and four children-at-home states. It excludes joint,
two-shock and fertility-nest experiments and permanent income levels. Numerical
parameters, wealth-grid size and income-node count are configurable. Transition
matrices must be finite, nonnegative and stochastic. The reference income grid
has 15 nodes; the guard also accepts valid nine-node processes.

The reviewed native exact-allocation contract is required. Compiled location
selection is required: the inherited absence of `use_loc_kernel` means its
reviewed `True` default, while an explicit `False` is rejected. Public primitive
helpers assume an already validated parameter object; use the guarded solver
entry points for model calculations.

## Actual changes

- Removed shared-clock branches from child-count dispatch and its selected
  shared, Bellman and distribution paths. Valid states retain
  `m = 0, ..., n` for `n = 0, ..., 3`.
- Removed joint-choice Bellman/distribution branches and the unused
  `joint_nested.py`, `two_shock_choice.py` and `fertility_nested.py` modules
  and imports. Removed the
  joint consumption reconstruction, which is unreachable under required exact
  allocation. Removed unreachable statements after selected returns.
- Removed nonsequential birth branches from the Bellman and distribution path.
- Replaced the compiled location kernel with a single-destination calculation.
  It does not enumerate destinations or read migration interpolation maps. It
  preserves the original shifted-value multiply, maximum, exponential and
  log-sum arithmetic to avoid changing floating-point results. The old Python
  location-choice fallback is removed.
- Preserved the separate `credit.py` binding and both explicit credit modes;
  corrected mode continues to require an explicit `D`. Credit formulas and
  borrowing/sale-solvency kernel formulas are unchanged.
- Copied the authenticated input loader and optional call counter so runtime
  imports are self-contained. The loader still accepts the original bundle
  schema and pins the original reference; this phase does not introduce an
  override/calibration interface.

## Remaining scope and limits

Singleton location axes remain in Bellman/KFE arrays, stats and reports.
Migration-map construction, location-indexed budget loops and ancillary
historical reporting helpers remain. Phase 1 removes the selected location
choice computation, **not all location architecture**. Removing array axes
requires coordinated kernel, KFE, statistic and observer-interface work.
No lifecycle replay, GE calculation, speed comparison, calibration, transition
or production promotion was run. Component equality is not a full solution
certificate. Internal population closure remains unchanged; an external
birth-renewal/population-scale workflow still owns its full closure.

The existing `run.py` accepts the inherited fixed-price/reference-GE interface.
Its inherited corrected-GE refusal remains, so use a separately reviewed full
closure driver for that experiment; this package does not certify that driver.

## Bounded verification

From repository root, using the existing Python 3.13 environment:

```sh
OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
PYTHONPATH=code/model output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python \
-m pytest -q code/model/experiments/stationary_single_market/tests/test_specialization.py
```

The suite checks authenticated input identity, source preservation, every valid
`n,m` child state, unsupported-structure rejection, configurable income-node
counts, exact shared precomputation, pure and compiled singleton-kernel equality
including dead-value boundaries, and credit bindings. Tests import the reviewed
`refactor_lab` only as an oracle. Runtime code imports neither it nor active
legacy model modules. See `verification_receipt.json` for the bounded result.

## Diff review map

1. `contract.py`: supported structure, no numerical parameter pins.
2. `engine/kernels.py`: singleton kernel, retained arithmetic sequence.
3. `engine/shared.py`, `parameters.py`: independent count and deleted dead paths.
4. `engine/household.py`: selected sequential Bellman and removed joint/fallback.
5. `engine/distribution.py`: same birth pools/flows, selected independent count.
6. `run.py`: package naming and guard; `credit.py`/`inputs.py` remain unchanged.
