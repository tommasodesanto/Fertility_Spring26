# SameSex2 negative-rooms diagnosis — plan (code/tests ready, no compute run yet)

Question: does national SameSex2 rooms RF = -0.0376 [-0.0476,-0.0276] pp
(job 18274624) reflect a coding/link/order problem, a direct sex-
composition/room-sharing channel, or a robust-but-nonidentified
association? Not assumed answered either way; no specification search for
a positive sign.

## New code (this pass, local tests only, production files untouched)
- `code/empirical/acs/kleven_pseudo/samesex_diagnosis_lib.R`
- `code/empirical/acs/kleven_pseudo/test_samesex_diagnosis.R` — 19/19 pass,
  including a real run of the metadata extractor on the existing local
  fixture (finds 159/3839 rows with a third-child age tie `a2==a3` — a
  real signal that this ambiguity is non-trivial even at tiny N).

## Stage-regression matrix (≤100 budget; planned = 31)
| Item | Diagnostic | Regressions | Input needed |
|---|---|---|---|
| A | Reproduce SameSex2 pooled ROOMS RF/FS vs saved receipt, tol 1e-8 | 1 refit | `analytic_frames.rds` (saved, on compute) |
| B | Event age 0:5 × {rooms, bedrooms, ownership} RF+FS, all ages kept | 18 | `analytic_frames.rds` `ss` (saved, on compute) |
| D | Baseline same-sex + first/second-child-sex controls (3 outcomes) | 3 | New metadata (sex1/sex2) merged onto `ss` |
| D | Both-boys-vs-mixed, both-girls-vs-mixed groups × 3 outcomes | 6 | Same metadata |
| E | Ambiguity sensitivity excluding `a2==a3`, 3 outcomes | 3 | Same metadata (needs `a3`) |
| — | Sex-cell descriptive counts (BB/BG/GB/GG, weighted+unweighted N/D/Y) | 0 (descriptive) | Same metadata |
| **Total** | | **31** | |

## Input availability (checked this pass)
- **A, B**: fully covered by the already-saved national `analytic_frames.rds`
  (path in `analytic_frames_identity.json`, MD5-checkable) — no new
  extraction needed, compute-only (RDS not read locally).
- **C (metadata for D/E)**: `sex1, sex2, a1, a2, a3, linked_child_count,
  NCHILD_norm` are **not** in the saved slim `ss` frame (it was projected
  to only the regression's own columns). They exist in the pre-projection
  roster object built by the *unchanged* `build_mother_roster` +
  `build_samesex2` functions. `extract_samesex_roster_metadata_one_state()`
  (new, tested) reproduces the exact production per-state pipeline
  (`apply_source_sample_gate` → `build_mother_roster` → `add_outcomes` →
  `apply_oldest_child_minor_gate` → age filter → `build_samesex2`) and
  keeps only the small metadata columns + keys, releasing each wide state
  table immediately — never concatenating raw person-level data across
  states. `verify_key_coverage_and_recompute()` (tested: catches both a
  missing household and a recomputed-field mismatch) enforces exact 1:1
  key coverage against the cached `ss` and equality of recomputed
  `samesex`/`treatment_3plus`/`event_age` before any D/E diagnostic runs —
  fails rather than silently changing the sample.

## Audits before any new fit (per-state, reported not filtered)
`sex_order_audit()`: invalid sex1/sex2 codes, `a1<a2`/`a2<a3` order
violations, **`a2==a3` third-child age-tie frequency** (distinct from the
already-excluded `a1==a2` primary tie — this is the case where selecting
"child 2's sex" is potentially arbitrary among 3+ same-age children),
`NCHILD` vs `linked_child_count` mismatches, duplicate household/person
keys. Reported as counts, never used to silently drop rows beyond the
already-approved sample rule.

## Budget
- Compute: one ≤4CPU/64GB/15min smoke (verify metadata extraction +
  reproduction gate + a small slice of B/D/E) before one ≤4CPU/128GB/60min
  full diagnostics run, per your instruction. Thread caps 1, 5-min
  heartbeat outside a fresh OUTDIR, atomic per-case receipts, no blind
  retry.
- Estimated wall-clock: metadata extraction ~26min (same per-state pass
  as the national construction, dominates cost); the 31 regressions
  themselves are cheap (RF+FS only, no AR) — order tens of seconds total
  based on the national run's per-fit timing.
- **No cluster submission this pass** — code/tests only, awaiting your
  review and go-ahead on the guarded launcher below.

## Deliverable format (not yet built)
Small PDF report + figure + CSV + MD under
`output/acs_fertility_iv/samesex_diagnosis/`, clearly separate from
`RESULTS.md` (which is preserved unmodified).
