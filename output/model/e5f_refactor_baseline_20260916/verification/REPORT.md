# E5F before/after verification harness

The comparator is an exact, fail-closed check for explicit artifact files or
manifest-listed JSON, CSV, NPY, and NPZ bundles. Bundle manifests require the
supported schema `e5f_refactor_evidence_v1`, nonempty string source and input
identities, matching input fingerprints, safe existing paths, and SHA-256
claims that match the artifact bytes. Distinct reference and candidate source
identities are accepted only when explicitly supplied to the command. Bundle
comparison ignores only `provenance.source_identity`, direct artifact
`path`/`sha256` claims after those claims are verified, and the fixed
administrative allowlist inside top-level `run_metadata`. Artifact JSON is
otherwise exact, including timestamps, durations, booleans, numeric types, and
unknown metadata. Arrays require exact shape, dtype, and values; matching NaN
positions are allowed, while infinities must match exactly. NumPy loading uses
`allow_pickle=False`.

The unittest suite contains 12 focused methods with fresh temporary bundles for
each mutation. It covers positive distinct-source bundles with differing
manifest timings; source and input mismatches; array value, shape, dtype, and
tiny-number changes; discrete branch, boolean/integer, and numerical-gate
changes; exact JSON timestamp/duration/run metadata; unknown manifest and units
metadata; missing keys/files, stale hashes, path escapes, and self-checks;
object NPY/NPZ and corrupted NPZ; valid and changed-cell CSV; and structured
CLI errors for malformed JSON, invalid UTF-8, and missing or wrong-type paths.
It also checks that an artifact-only file named `manifest.json` receives exact
JSON treatment and that the reserved bundle artifact name is rejected.

Validation run:

```text
RUNTIME=/Users/tommasodesanto/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3
$RUNTIME -m py_compile code/model/tools/verify_e5f_refactor.py code/model/tools/test_verify_e5f_refactor.py
PYTHONPATH=code/model/tools $RUNTIME code/model/tools/test_verify_e5f_refactor.py -v
Ran 12 tests in 0.245s — OK
```

The harness certifies artifact or artifact-bundle equivalence only. It does not
replay the native E5F model, establish source-history equivalence, or certify a
model solve, compiled replay, cluster job, or production calibration.
