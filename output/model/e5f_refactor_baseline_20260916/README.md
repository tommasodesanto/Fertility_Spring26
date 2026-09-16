# E5F refactor baseline package

This package records the six-stage behavior-preserving refactor workflow. Stage 1 inventory is complete, but freezing the baseline remains pending restoration of the exact certified recovery source/input bundle. Stage 2’s reviewed dependency-light focused harness passes 12 tests; native cross-version replay remains pending restored inputs. Stage 3 Kimi profiling review is complete (provider-reported completed-step cost 0.415824), while native benchmarks remain deferred.

Read [baseline/REPORT.md](baseline/REPORT.md) and [baseline/manifest.json](baseline/manifest.json) for source hashes, availability, and dirty-tree handling. The verification harness is documented in [verification/REPORT.md](verification/REPORT.md) and implemented at `code/model/tools/verify_e5f_refactor.py`. The profiling assessment is [profiling/REPORT.md](profiling/REPORT.md); raw Kimi streams remain evidence only and are not part of the source package.

Once two explicit bundle directories exist, run the comparison harness from the repository root:

```bash
python3 code/model/tools/verify_e5f_refactor.py <reference-bundle> <candidate-bundle> --output-dir output/model/e5f_refactor_baseline_20260916/verification/run --reference-source-id <reference-id> --candidate-source-id <candidate-id>
```

The harness requires each bundle's manifest and artifact hashes, permits separately pinned source identities, and requires matching input/scientific contracts. No production source changes are authorized in stages 1–3.
