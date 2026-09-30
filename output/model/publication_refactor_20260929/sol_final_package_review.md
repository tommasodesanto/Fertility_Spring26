# Final package wiring/provenance review (read-only)

**Disposition:** No numerical-source mismatch found. I hashed all 18 active `code/model/refactor_lab/engine/*.py` files against the immutable `output/model/publication_refactor_20260929/indexed_src_gridfix/refactor_lab/engine/` stage; every file is byte-identical, including promoted `kernels.py` SHA-256 `fe7d43af…`. I recomputed every local code/transform hash named in `engine/promotion_receipt.json` and its four evidence-file hashes; all match. The 27-component-test JUnit has zero failures/errors, and saved driver smoke summaries show local/Torch pass, fail-fast, timeout and GE budget routing. These are wiring checks; full Torch GE/reporting job 18851943 was still running at review time.

**Correct the provenance description before final presentation.** `README.md:71-79` says engine bodies are “the executed original bytes” and “the only numerical change” is indexed saving. `engine/materialize_receipt.json` instead names the separately reviewed *credit-overlay* source for `parameters.py`, `solver.py` and `kernels.py` (overlay SHAs `d7b4d23c…`, `cefc1627…`, `379d179a…`), then the indexed transformation changes one kernel body. State this two-stage provenance explicitly: the reference mode preserves the frozen borrowing behavior; the package also contains the separately labelled corrected-credit rule; indexed saving is the sole **performance transformation after overlay materialization**, not the only numerical difference from the original source tree. The receipt chain itself is coherent; this is a wording correction, not a reason to rerun solves.

**Clarify the regeneration claim.** `engine/promotion_receipt.json` says `regenerated_from_live_engine_identical=true`, echoed in `opus_report.md:8-9`. The current live engine is already indexed, while `verification/make_indexed_stage.py:56-63` expressly requires the *untransformed* scalar definition and no helper collisions; invoking it on the current live engine would stop. If regeneration was performed from a preserved pre-promotion/materialized stage, name and hash that input; otherwise remove/rename this Boolean. The active and tested-stage hashes, diff, saving source and generator hashes remain verified.

**Tighten the normal/verification boundary or soften its claim.** `README.md:69` and `opus_report.md:27,40` say a normal run imports no verification module unless `--count-calls` is selected. `run.py:187-190` imports `.verification.callcount` on **every equilibrium run**, even when the option is off. The solve still uses only `refactor_lab.engine` numerical modules and does not import the old model, so this is a layout/independence correction: move the import behind the flag (while keeping the exception handling valid), or describe the unconditional lightweight import accurately. Fixed-price normal execution does not have this issue.

I left the known README shell-command corrections to Luna as assigned. Existing historical receipts naming old verification paths are evidence for their tested stages, not current-path wiring evidence.

**Launch wrapper follow-up (Luna owns fix):** `verification/verify_torch.sh:11` locates `verify.sh` using `dirname "$0"`. Slurm `sbatch` runs its spooled copy, so this resolution fails from the spool directory despite the ordinary Torch dry smoke. Require an explicit staged `LAB_SRC` path and execute `$LAB_SRC/refactor_lab/verification/verify.sh`, then smoke the wrapper from a relocated script path. The already running GE batch used explicit staged paths and is unaffected.

## Closure after Luna's corrections

- **Wrapper resolved.** Current `verification/verify_torch.sh:10-17` requires `LAB_SRC` and executes its staged driver. The saved `pass10_wrapper_check` copied wrapper matches current bytes; the fake spool case forwards `MODE=torch` and all arguments, and the missing-variable case fails explicitly. This was a no-solve check.
- **Regeneration metadata resolved.** `engine/promotion_receipt.json` now says `regenerated_before_promotion_identical`, names `scalar_src_gridfix/refactor_lab/engine`, and pins its input kernel SHA. I recomputed that preserved kernel hash: `19dceb70…`, matching the receipt and transform-before hash. No current-live regeneration is claimed.
- **Normal import boundary described accurately in the README and `verification/__init__.py`.** They now state that `run.py` imports the lightweight `verification.callcount` in equilibrium, with instrumentation inactive unless requested. The numerical engine remains independent of the old package. `opus_report.md:27,46` still contains stale sentences implying the import occurs only with the flag; line 40 of the same report is correct.
- **Provenance improved but first sentence remains overbroad.** `README.md:75-83` now names the reviewed three-file credit overlay and the subsequent indexed performance transformation. Its opening sentence at line 73 still says the engine bodies are the executed original bytes, which is false without “except the reviewed overlay and indexed kernel” or equivalent qualification. Revise that sentence and scope `promotion_receipt.json`'s “no other numerical edit” to the *post-overlay performance promotion*.
- **New one-core command typo:** `README.md:46,110` sets `NUMBA_THREADS=1`, which Numba does not use. Replace with `NUMBA_NUM_THREADS=1` to make the displayed manual command genuinely one-core; `verification/verify.sh:110` already uses the correct variable.

These are documentation/launch-command corrections; none alters the tested engine bytes. Source identity is accepted for the isolated test package. The separate full Torch GE/reporting result is still pending.

## Lead closure after final review

The remaining documentation findings are corrected: README uses
`NUMBA_NUM_THREADS`, states the reviewed credit overlay explicitly, and
describes the unconditional lightweight counter import accurately. The Opus
report now makes the same distinction. The relocated-wrapper fixture passes
and rejects missing `LAB_SRC`. No numerical source changed. The subsequent
GE18851943 comparison passed87native arrays, identical14/31tables and17actual
PNG file hashes; see `ge_pair_18851943/`. Final18engine.py files remain
byte-identical to the tested indexed stage.
