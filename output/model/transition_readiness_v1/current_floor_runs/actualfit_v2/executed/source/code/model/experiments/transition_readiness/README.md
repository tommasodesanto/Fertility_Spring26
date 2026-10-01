# Transition readiness preparation

This directory owns isolated preparation, not a model/reference change. Source
pins in `legacy_source_pins.json` are exactly the ten deployed September 29
budget-diagnostic-v2 sources. `pinned_tools/` preserves their original bytes.
The legacy manifest is September 28 block0506 with
psi_child=0.1355551166583114 and manifest SHA
147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4.

## Historical numerical diagnostic

`legacy_changed_psi.py` performs the experimental psi_child ×1.001 change as
one permanent 2007 surprise. It retains earnings, preferences other than psi,
entry, wealth grid, credit/repayment rules, fiscal objects, closed-population
birth-to-entry timing and housing stock. This inherits the historical repayment
limitation and cannot validate current-floor borrowing rules.

The endpoint uses the **unchanged** `NativeEstimator.endpoint`: at most 24
stationary root evaluations including fresh final reproduction, plus one native
one-step verification. Its demographic-renewal, fiscal, estate, household,
projection and native terminal gates are preserved. There is no manufactured
endpoint, rescaled fertility or policy fallback. At most three six-date mappings
follow: changed-psi input; exact dated pension update
`b_new = b_old * payroll_revenue / pension_outlays`; fresh reproduction if the
trial passes. House prices remain at the historical reference price for this
bounded diagnostic. Failing housing or terminal gates causes an uncertified
result and stops. No price-path estimation occurs.

The internal hard budget is 3550 seconds within the 3600-second job, with endpoint stage at most 2000 seconds and
each native mapping at most 1200 seconds, all also bounded by the global
deadline. Maximum conservative policy calls: 24 stationary + 2 for endpoint
one-step + 3 × 2 × 6 = 62. Actual calls can be smaller. The prior three six-date
maps took 1170.57 seconds, leaving 2379.43 seconds for endpoint/setup/plotting;
changed-endpoint timings are unknown, so completion is not guaranteed. One
thread, 64-GiB cache and 96-GiB job memory are the intended resource contract.
Heartbeat is written every 30 seconds; latest-completed, best-so-far, map
checkpoints and complete/failure receipts persist. Capture saves dates 0, 3 and
5, and the unchanged standard 17 plots are rendered for the trial. If failure
occurs before a completed trial, no fresh visual packet is claimed.

Zero-solve check:

```sh
python3 code/model/experiments/transition_readiness/legacy_changed_psi.py --preflight
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s code/model/experiments/transition_readiness -p test_readiness.py -v
```

Real diagnostic command, only inside the reviewed Torch snapshot/binds with
NUMBA_NUM_THREADS, OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS and
BLIS_NUM_THREADS all set to 1:

```sh
python code/model/experiments/transition_readiness/legacy_changed_psi.py \
  --output output/model/transition_readiness_v1/legacy_changed_psi_run
```

No 104/128-date map, fit, repeated unchanged smoke or model job was run during
preparation. The lead owns deployment and launch approval.

## Eventual selected-calibration adapter

`selected_adapter.py` is a callable authentication seam, not a stationary or
transition engine. Existing floor searches have computed provisional GEs; those
are not an authenticated author-selected, fresh-repeat checkpoint. The floor
driver is `output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/run_psi.py`,
which imports sibling `utility_floor_round2_v1/runner.py`. Its corrected credit,
entry law, grid and preferences must be pinned from the actually selected run.
The old block0506 loader cannot silently supply that state.

`authenticate(manifest_pin)` first checks an exact pinned manifest with schema
`selected_calibration_transition_adapter_v1`, author_selected=true and
provisional=false. Required pinned artifacts are source_manifest (nonempty
path/SHA pins), target_contract (fingerprint), native_ge_receipt and native_selected_repeat
receipts emitted by the actual floor runner, effective_parameters (complete
JSON), state_identity, checkpoint and selected_price_repeat. The state receipt
must bind checkpoint, source and parameter hashes, entry/grid identities and
native population/birth-entry-queue verification. The fresh-repeat receipt
must pass recorded gates and bind those same identities and target fingerprint.
The full economic contract covers preferences, earnings, entry, grid, credit,
sale/mortality repayment, fiscal, population, housing, geography, estate, targets
and timing. Each row has identity, classification and pinned evidence;
outstanding closure objects block production. This authentication does not
itself establish the mathematical validity of a gate; lead review remains needed.

`setup_authenticated(manifest_pin, engine_contract=..., setup=...)` invokes a
real engine setup callable only after those checks, and only if its complete
contract matches. Native setup must then deserialize the authenticated
checkpoint and verify its arrays against the state receipt; no checkpoint is
unpickled before authentication. Housing stock versus elastic-reference choice
must be explicit; the same birth-renewal/housing-population mathematics does not
automatically authenticate a different entry/grid/credit engine.

Concrete preflight and deliberately blocked production CLI:

```sh
python code/model/experiments/transition_readiness/selected_adapter.py \
  --manifest /absolute/path/to/author_selected_transition_manifest.json \
  --manifest-sha256 VERIFIED_SHA256
python code/model/experiments/transition_readiness/selected_adapter.py \
  --manifest /absolute/path/to/author_selected_transition_manifest.json \
  --manifest-sha256 VERIFIED_SHA256 --production
```

There is intentionally no synthetic selected manifest. Missing selected-price
repeat/authenticated checkpoint, unclosed housing/estate objects, dated-engine
validation or full-horizon comparison remain visible blockers. The production
CLI always refuses numerical execution; the callable seam is the concrete
integration point for a reviewed engine.

The exact `run()` control-flow smoke replaces only native endpoint/mapping/render
dependencies and verifies three actual loop calls, three checkpoints, heartbeat,
latest/best receipts and the production-blocked result. Its receipt is explicitly
zero-model-call orchestration evidence. The Torch launcher repeats this smoke and
source-import preflight before any real endpoint call. An independent SIGUSR1
watchdog enforces the absolute deadline despite unchanged inner SIGALRM guards.
