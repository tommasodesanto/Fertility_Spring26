# Normalized calibration v2 deployment

The v2 search restarts from the six best postcheck-verified v1 centers, each with its exact center and three deterministic nearby starts. No v1 final simplex was serialized, so this is an explicit restart rather than a simplex resume. Verified seed provenance is in `verified_seeds.json`; the best candidate's full ROOT/REPEAT verification and 17-plot packet is under `best_v1_candidate/`.

The frozen 120-by-9 model, original target weights, free `psi`, internally derived `H0`, and population normalization `N=1` are unchanged. Each of 24 chains has one CPU, 24 GiB, one thread, a 10,800-second actual-start budget, at most 150 objective calls, a 900-second final reserve, and the 700-second native guard. There is no automatic retry.

Staging archive SHA-256: `fa05236a51802d68fdf49b684dd1b09a99d24edf499fc4f4249e527d5275acd2` (240 packaged source files; all 234 source pins included and verified). The staged package is at `/scratch/td2248/projects/normalized_calibration_v2`.

Submitted 2026-10-01 at 19:48 EDT under `torch_pr_570_general`: incumbent native smoke gate job `18989877`; 24-chain array job `18989878`, dependency `afterok:18989877`. The submission receipt is `submission_receipt.json`. At 19:50:21 EDT, the smoke gate was RUNNING on `cs609` (elapsed 18 seconds) and the array remained PENDING on its `afterok` dependency. The gate output file existed but was still empty; no smoke result or production chain had started at this snapshot.
