# Recovery smoke review — job 18271160 (2026-09-22)

**Result: COMPUTATIONAL GATE PASSED.** COMPLETED 0:0, elapsed 00:05:59, batch MaxRSS 5,084,664K (~4.85GiB).
- Stage 1 (real VT driver, STATEFIP_LIST=50): 19.068s, 18/18 `full_fit`.
- Stage 2 (scaled probe, VT frames replicated to national scale): 337.304s, gate PASS 6/6, Twin1 N=857,670 (target 857,607), SameSex2 N=603,985; actual RF/FS/IV nobs matched the independently-recomputed expected complete-case N for every case.

## VT fit_receipts (18/18) — checked
Nobs internally consistent (rf==fs==iv per case): **TRUE**. AR error-free (ar_n_errors==0, no ar_error field): **TRUE**. Named-coefficient count vs full V dimension match: **TRUE**. V symmetric: **TRUE**. V diagonal (instrument row) matches the table's reported SE: checked and consistent.

These are **internal** consistency checks against each receipt's own recorded fields — not an independent recomputation from raw VT source frames (that would require reading `analytic_frames.rds`/raw partitions, excluded from this small-artifact collection). Confirms self-consistent arithmetic, not external ground truth.

## Warnings
168 occurrences in the log of: `fixest: The VCOV matrix is not positive semi-definite and was 'fixed' (see ?vcov)`. No other errors/failures anywhere in the log. The numerical cause and magnitude of the non-PSD condition are unresolved -- not diagnosed here. fixest's automatic correction restores a valid, self-consistent V (symmetric, correctly dimensioned, diagonal matching the reported SE, as checked above), but that post-repair self-consistency does not quantify how much any specific V changed from an unfixed version, which was not computed. This is a real, unassessed caveat on **uncertainty quality** (SE/CI/AR) that remains provisional; it does not block the RF/FS point estimates from proceeding.

## Resource conclusion for national
- All-state construction (prior 18247804 evidence): ~26 min.
- VT's 6 full outcome×design fits: 5.5 min scaled-probe-equivalent time (337s for 6 realistic-N cases).
- Smaller event-age-3/5 subsets: proportionally faster.
- **Rough total estimate: ~40–60 min**, well inside the 3h/128GB allowance.
- VT's scaled probe does NOT measure the true national roster's memory footprint (real national has more distinct FE levels, more states, true clustering structure) — **128GB is retained**, not reduced, given this gap.
