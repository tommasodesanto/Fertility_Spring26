# Fixed unsecured-credit contract v1

Reference: `2007 stationary reference — block0506, September 28 verified export`.
The before manifest pins `parameters.py` at `66f86697c2c58ca3864305bf13dd2be71a008905b2beb573f1a4ebafabef5464`, `solver.py` at `b637a655a9344b63f4461ee0fa4796c04bd98188477c4e6ace2c48ae0fc8aec1`, and `kernels.py` at `639c9a21797dbc9f2a0e9a891f283c115353c2edfcb89c959a7fe9f32b86ca27`. `prepare_overlay.py` refuses either a hash mismatch or an existing destination.

`P.unsecured_credit_limit=None` keeps the legacy renter rollover-and-age-taper path. An explicit finite scalar \(D\geq0\) instead imposes renter saving \(b'\geq-D\), independently of current unsecured debt, age, income, and the existing taper arrays. At terminal age or with positive current death probability, the separate non-negative-estate restriction gives \(b'\geq\max\{-D,0\}=0\). Owners, purchases, legacy debt arrays, and mortality/estate rules are unchanged. `native_solvency_credit` and the scalar contract fail fast when combined; legacy/factored Bellman routes also fail fast for an explicit scalar, so the supported route is `solve_bellman_full_markov_income`.

An explicit scalar, including \(D=0\), turns on the owner-to-renter gate: raw liquidation wealth \(b+(1-\psi)pH\) must be non-negative before interpolation or price-grid clipping. The deterministic and logit tenure kernels and their Python fallbacks apply the same gate. The native renter saving kernel receives a separate scalar floor; no owner input array is changed.

`overrides.json` records the planned case \(D=0\); no positive magnitude is selected. The parent fixed-reference manifest SHA is `147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4`; the checkpoint SHA `b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d` is recorded identity only and was not downloaded or rehashed. Preferences, entry, fiscal objects, and prices are inherited unchanged; there is no run-controller or GE acceptance.

From the repository root, run the pure checks with:

```sh
NUMBA_DISABLE_JIT=1 code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/credit_no_taper_v1/fixed_credit_contract_v1/test_contract.py
```

They cover legacy `None` identity, zero and positive limits at several ages and
current balances, estate floors, invalid values, external overrides from old
parameter objects, direct rebuild validation, natural-credit conflict, sale
boundaries and an actual saving-kernel borrowing fixture. Owner kernels and call
sites pass AST comparisons. The separate lead checks below execute the
Python fallback and unsupported-route guard. No lifecycle, equilibrium or
checkpoint is read by these tests.

Known inherited-entry blocker: the earlier strict-\(D=0\) numerical run failed with two negative-wealth entrants. This packet deliberately does not change entry.

## Lead review and bounded compilation check

The lead reviewed every economic source change against the two constraints,
and independently executed the actual Python `_tenure_location_stage` fallback
from its source: deterministic and logit choices reject negative raw sale
balances and accept the exact zero boundary. The unsupported legacy Bellman
guard also executed successfully. This small check took 47 milliseconds and
performed zero lifecycle solves.

The first draft mistakenly combined the explicit floor with the old floor,
preventing positive borrowing from zero wealth. It was rejected before any model
run and is preserved in `overlay_initial_unaccepted/`. The corrected saving
kernel replaces the old renter floor when the scalar is supplied; its executable
fixture confirms a positive cap permits borrowing, while cap zero prohibits it.
Owner kernel and owner call sites pass exact AST comparisons against the reference.

Torch zero-lifecycle compiled smoke **18843355** was submitted once with a
five-minute limit, one CPU and 16 GiB. Source and test pins are checked before
execution. Remote sources are read-only; receipts go to `verification_v1/`.
This submission does not resume the expired numerical comparison budget.
