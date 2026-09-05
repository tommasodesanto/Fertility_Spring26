# E5F code architecture and correctness audit evidence

5 September 2026. Read-only independent review; no model runs, calibration or policy changes.

The report is `docs/model/e5f_full_code_correctness_efficiency_review_20260905.md` in the repository. The production source is the frozen fifty-file snapshot under `tmp/e5f_overnight_20260905a/frozen_production/code/model/`, checked against the selected task's complete source manifest. The current checkout adds a default-off young-ownership file and edits the driver; the report distinguishes this from the production snapshot.

## Contents

- `architecture_summary.json`, `footprint.csv`, `functions.csv`: current initialization imports, AST function sizes, line counts, current versus frozen comparisons, production and current manifest objects.
- `footprint_reduction.json`: interpretable footprint groups, inactive non-Markov dispatch functions, byte-identical copied modules. Import reachability is not execution coverage.
- `source_verification.json`: independently reconstructed frozen bundle fingerprint and all source comparisons.
- `side_channel_probe.json` / `.py`: the same policy bundle and inherited state produce different births when continuation probabilities on the parameter object change. It executes the actual fertility and evaluation functions with identity mocks for feasibility preparation and current-choice realization. It establishes a latent API defect, not an affected production outcome.
- `bracket_probe.json` / `.py`: a residual mock proves 18 repeated evaluations at a reached price bound. No household solve runs.
- `invariant_checks.json` / `pure_invariant_checks.py`: down-payment threshold, Bellman/forward adjoint, one birth per period, and twenty-year queue impulse checks.
- `saved_state_footprint.json` / `.py`: saved array sizes, zero invalid-family-state mass, all55 policy dates' fallback flags and hashes, and wrong verbose fertility formula example.
- `existing_test_receipt.json` / `existing_tests_direct.py`: fifteen existing pure tests invoked directly against frozen scientific modules. This is not a passing pytest/compiled suite.
- `verification_manifest.json`: artifact hashes and source/report verification.

## Reproduction

Commands run from:

`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`

Use system Python3.9.6/NumPy2.0.2 for `side_channel_probe.py`, `bracket_probe.py`, and `saved_state_footprint.py`, whose no-solve operations do not require the full parameter builder. Use:

```bash
NUMBA_DISABLE_JIT=1 code/model/.venv/bin/python output/model/e5f_code_architecture_audit_20260905a/pure_invariant_checks.py
NUMBA_DISABLE_JIT=1 code/model/.venv/bin/python output/model/e5f_code_architecture_audit_20260905a/existing_tests_direct.py
```

The second runtime is Python3.10.10/NumPy1.24.3/Numba0.58.1 with JIT disabled. The standalone footprint importer is `footprint.py`; it intentionally inspects the current checkout's initialization route and is therefore a dated inventory rather than a frozen execution trace.

The pytest CLI exited139 without output in three bounded attempts (ordinary setup, disabled JIT, then disabled automatic plugin loading). No pytest pass or compiled-kernel validation is claimed. Direct execution of fifteen selected test functions subsequently passed. The system Python cannot run the full parameter setup because it contains runtime union-type checks, and the bundled Python3.12 lacks matplotlib needed by the transition's imported diagnostic module. No dependency environment was changed.

The probes require the retained frozen snapshot and selected output packets. Their data reads are explicit. Their output files overwrite only their own named receipts. Successful probe source is retained here so reproducibility does not rely solely on scratch files.

## Successful standalone source receipts

Every script below completed successfully in its stated runtime. The checksums pin the reviewed source; generated receipts are separately hashed in `verification_manifest.json`.

| Script | SHA-256 |
|---|---|
| `footprint.py` | `0d3959f896b67705b0f48bc5ded86e90f5aa1bfe023fcd2f364a5bfc06aa257b` |
| `side_channel_probe.py` | `9df74b850a108d96f9cf9a6dde27a46f068aa8d5676e7f0167df70e8116e8aee` |
| `bracket_probe.py` | `f5a6575581e349e31b61546ecb3bf99bc1857536962c1e18847e1b20018e92c4` |
| `pure_invariant_checks.py` | `5db7ad6fb92fd70237377981028c7d5f98c9cd0854a7abd8ccaf0a1c132b1399` |
| `saved_state_footprint.py` | `c6fa87d67bee0ebed1a7ec62dd33aea8f780eee48258ef97f39ab00b6a2537ca` |
| `existing_tests_direct.py` | `a76d77532d14454eed0b7ecb0aca3f05053782bc4e033ea1abe7bbc678619228` |
