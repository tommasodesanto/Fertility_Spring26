# GE pair 18851943 collection

Collected from `/scratch/td2248/projects/publication_refactor_20260929/results/ge_pair_18851943` after Slurm reported completion with `rc=0`, phase `all_passed`, and 627 s recorded by `completed.json`. The scheduler receipt records 631 s elapsed, one allocated CPU, 12G requested, and batch MaxRSS 4,962,448K (~4.73 GiB).

The compact package includes the run plan/progress/summary/step and cache receipts; strict array comparison and full per-array rows; lab GE receipt and original-engine solve receipt; both certificate receipts, gates, parameter-identity receipts, observer JSON and CSV tables; diagnostic summaries; and only the lab certificate's 17 standard PNGs (1,754,914 bytes total). It excludes solution NPZs, acceptance pickle, Numba caches and raw streams/logs.

**Mechanical matches verified locally:** strict comparison passed with `strict_paths=true`, 87 array paths compared and 87 marked exact, no missing or extra lab arrays; these 87 arrays do not include the reconstructed pre-fertility distribution. Each GE certificate separately checks the stationary reconstruction and operator; the earlier fixed-price 113-path comparison includes that distribution. Fit tables have 14 data rows and identical CSV bytes (SHA-256 `5f5866ac59829bb8085b53af8c8f2a62324a097eb1c9136a6ba75f06796af1a8`); parameter tables have 31 data rows and identical CSV bytes (SHA-256 `9ecb45e649c86328861e047e5e6394b543303cc81dbc97d3bd1bd9c289c9c269`). The reported GE price is exactly equal across both solve receipts and both certificate receipts: `0.789870850567871`. Both receipts report 17 plot hashes and the maps match; SHA-256 of each collected lab PNG also matches the lab receipt's per-file hash. The receipt field `plots_hash_checked` is false on both sides; therefore this records equality of the reported maps and lab-file/receipt consistency, not an independent frozen-reference plot check.

Initial parameter identities match in both engines; final parameter identity receipts report no unexplained fields, and old-vs-lab final identity is true after excluding only `eq_iter` and the run-specific inherited-distribution evidence directory. The certificates share status `ge_gates_passed_renewal_consistent_not_equivalence`. Four probability-array nonfinite counters are zero on each side; consult the saved gates and observer JSON for the full diagnostic details.

Engine solve-stage times were 204.870 s (lab) and 275.951 s (old engine). Driver step times: lab_ge 219.584s, lab_ge_certify 66.526s, old_ge 326.777s, compare_lab_old 2.816s. Wrapper wall times include imports, certificate/reporting, and serialization and are kept separate in `verify/steps.tsv`.

No economics interpretation is made here; the root agent owns review of the saved gates and results. Machine-readable checks are in `matched_reporting.json`.

## Lead closeout

The lead additionally hashed the actual original GE PNG files remotely and
compared them with the collected lab PNGs and both receipts: all 17 match.
See `plot_file_comparison.json`. This is old-versus-new GE equality at the
common solved price; `plots_hash_checked=false` only denotes that the oracle
did not compare GE pictures with the different fixed-price checkpoint.
Native market residual 2.2985e-6 and pension residual 3.17e-13 pass; no occupied
budget, purchase, principal-floor or feasibility-projection mass violations
were found. The estate ledger retains its original provisional valuation and
does not certify counterparty settlement. No economics was changed.
