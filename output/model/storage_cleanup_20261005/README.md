# First storage cleanup — October 5, 2026

Completed the author-authorized fast cleanup of one historical Torch experiment:
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/nightpair_20260925_v1`.

Removed **680 nonbest `initial_state.pkl.gz` files**, releasing **110,037,324,288 allocated bytes (110.037 GB)**. The directory decreased from **117.940 GB to 7.903 GB**. A 4,096-byte scratch write and fsync passed; the temporary probe was removed. `myquota` has not reflected the deletion in its displayed account total; no updated account-wide usage is asserted.

Kept all 40 workers’ best cases, including both named selected reference cases; all selected exports, smoke/repeat runs, source, inputs and runtime; and every trial’s receipts, parameter/fit records, tables, plots and logs. All 720 point receipts share one objective fingerprint and report verified experimental points. No model, calibration acceptance gate, active transition evidence or pipeline was changed, and no model job was launched.

The lead independently derived the 40 worker minima from the saved receipt extraction and checked that the manifest is exactly their complement. Every deletion had a unique, single-link regular file and unchanged metadata guards. All 40 retained checkpoint metadata records matched before and after. This was a metadata check, not a new content-hash or scientific audit of retained arrays.

[Concise completion receipt](cleanup_summary.json) pins the local detailed evidence: proposal, exact path/stat manifest, applied script, actual receipt extraction, streamed deletion receipts, retained-checkpoint metadata and size/quota checks. Large model arrays remain in their original locations. **The cleanup has been applied once; do not rerun the execution command.**

Proposed future calibration retention, not implemented here: compact parameters, loss, model moments, solve status and provenance for every guess; full arrays only for current and previous best per chain plus selected final/reference cases; lightweight optimizer restart checkpoints; publish a new best atomically before retiring the previous saved version. Keep the existing scientific gates. Broader data organization and pipeline work remain a separate discussion.
