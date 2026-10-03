# Audited correction after Fable's same-session revision

**October 3, 2026, about 07:10–07:13 New York.** The same-session Fable 5.1 revision exited successfully at 07:09:45; its `exit_receipt.json` records exit code zero. The independent [audit receipt](../lead_review/AUDIT_RECEIPT.md) then identified a closing-screen arithmetic error and terminal-age overstatement. This note records the narrow follow-up performed within `fable_analysis/`; no model or equilibrium solve was run.

The executed soft kernel uses `dp_choice=(dp_arr-income_for_purchase)/R` and compares beginning liquid wealth `b` with that threshold (`code/model/experiments/purchase_timing_sandbox/source/refactor_lab/engine/household.py` around lines 894–901; `engine/kernels.py` around lines 315–329). The screen is therefore $Rb+y\ge(1-\phi)Q$, or equivalently $b+y/R\ge(1-\phi)Q/R$. `analysis/native_grid_analysis.py` had used $b+y/R\ge(1-\phi)Q$. I corrected only this screen, reran the **saved-array** script, and refreshed its four figures and JSON. The command was:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MPLBACKEND=Agg code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/analysis/native_grid_analysis.py
```

The script completed with exit code zero. It emitted an existing zero-mass wealth-tercile division warning; the calculation of the closing-screen pass shares was unaffected. Comparing the new JSON with the preserved pre-revision JSON shows **exactly four changed scalar fields**, all `share_passing_closing_screen`. The original arm now has 0.971/0.951/0.992 at ages 18/22/26; the alternative arm has 0.971/0.952/0.933. The memo's formula, table, range, and cross-arm sentence were corrected. No other JSON result changed.

The memo now limits the under-2% owner at-floor statement to ages 38–78 and flags age 82 separately (30.5% original, 70.9% alternative). It says “at the floor” because policy equality within $10^{-6}$ is not a counterfactual test of whether the floor changes choices. The historical hard/quarter language now distinguishes closing resources $A=b+y$ from a beginning-wealth-only cash test. Figure F1 now labels the 0.676 value as conditional on held first-birth flows; F2 identifies CPS completed interview age 25 as $[25,26)$ with the model's 0.875 interpolation. The memo adds the extensive-margin accounting scenario (0.349 early first-birth share would imply any-birth share about 0.504 versus CPS 0.457; holding the latter requires early first-birth share at least about 0.403 with certain second births). These are conditional calculations, not impossibility claims.

The original launch snapshot manifest was not altered. [Saved-array source hashes](saved_array_source_receipt.json) separately pin the two exact-repeat NPZs used by the analysis. The pre-revision memo, code and JSON remain in this folder under `pre_revision_manifest.json`.
