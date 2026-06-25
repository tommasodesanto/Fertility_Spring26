# Markov room-moment collapse bug — investigation, fix, verification

Date: 2026-06-25. Scope: the renter room-distribution moments in the
one-market intergen model (`code/model/intergen_housing_fertility/`).
Method: systematic debugging (root cause → reproduce → minimal fix → verify).

## 1. Symptom

The implemented Markov statistic reported childless renter share rooms ≥6 as
**0.013** vs ACS target **0.138** — a large apparent miss — while a full
state-by-state enumeration of the same cached solution gave **0.124**
(≈ target). The renter mean rooms (4.455) was identical both ways.

## 2. Root cause (confirmed)

`compute_markov_statistics` (solver.py) collapses the income dimension *before*
the generic `compute_statistics` applies nonlinear operators to the renter
policy:

```
g_total  = sum_z g                                  # mass: linear, fine
hR_total = collapse_markov_policy(hR, g, z_weights) # = E_z[h_R] per (b,j,n,cs)
stats    = compute_statistics(g_total, ..., hR_total)
```

`compute_statistics` then computes renter threshold shares and the renter median
on `hR_total`. For a continuous renter policy that varies across income states
within a cell, this is a Jensen error:

\[
\mathbb{E}_z\!\big[\mathbf{1}\{h_R \ge t\}\big] \;\neq\; \mathbf{1}\{\mathbb{E}_z[h_R] \ge t\},
\qquad
\mathrm{median}_z(h_R) \;\neq\; \mathrm{median}(\mathbb{E}_z[h_R]).
\]

Means survive the collapse (\(\sum_b g_{tot}\,\mathbb{E}_z[h_R] = \sum_{b,z} g\,h_R\)),
which is why renter/owner mean rooms and the cross-tenure mean gap were correct.
Owner room moments are rung-based (\(H_{own}[ten-1]\), income-independent), so
they are also correct. Only **nonlinear-in-renter-policy** moments break.

## 3. Affected moments (cache, p=0.8144)

| moment | role | reported (buggy) | corrected | note |
|---|---|---:|---:|---|
| `prime30_55_childless_renter_share_rooms_ge6` | hard target in room sets (w≈25) | 0.013 | **0.124** | ≈ target 0.138; the "miss" was the bug |
| `prime_childless_renter_median_rooms` | hard target in room sets (w≈10) | 4.217 | **4.144** | mild |
| `renter25_45_all_cap_share` | diagnostic | 0.010 | **0.228** | |
| `renter25_45_current0_cap_share` | diagnostic | 0.008 | **0.089** | |
| `renter25_45_current1_cap_share` | diagnostic | 0.021 | **0.621** | one-child renters press hard on the cap |

Unaffected (verified identical): renter mean rooms 4.455, owner mean rooms
6.835, and the **active-target** `owner_minus_renter_mean_rooms` gap 2.380.
So the **current** diagnostic point's loss (target set
`candidate_replacement_young_old_roomgap_v1`, which weights only the mean gap)
was **not** corrupted. The bug is a landmine for any room-share/median-targeting
set and for the reported target-fit tables.

## 4. Fix

Minimal and surgical (solver.py): added `markov_renter_room_moments(g, hR, P)`
that evaluates the five renter threshold/median moments on the full
income-resolved distribution, and overwrote the five corrupted attributes in
`compute_markov_statistics` after the `compute_statistics` call. The shared
non-Markov `compute_statistics` is untouched (its renter data is genuinely
z-free, so its threshold moments are already correct — verified there is no
other collapse-then-threshold caller).

## 5. Verification

- `reproduce.py`: buggy-path reconstruction matches the stored cache values
  exactly for all five moments (root cause faithful), and the corrected path
  matches the full-enumeration values.
- Unit: `markov_renter_room_moments` on cached arrays returns the corrected set.
- Integration: `compute_markov_statistics` via the fixed path returns the
  corrected five moments and leaves means / mean gap identical to 1e-6.
- `compileall` OK; `cli smoke --quiet` re-solves end-to-end with no error.

## 6. Implications

- **Corrects the room-separation narrative.** The model does *not* badly miss
  renter ≥6 (true 0.124 ≈ target 0.138). The earlier claim that cap=6 makes
  renter ≥6 "mechanically unreachable / incoherent" was based on the buggy 0.013
  and is retracted: the cap delivers the at-cap mass that matches the target.
- **Re-read the target fit before any room calibration.** Any prior cluster run
  or panel that weighted `renter_share_rooms_ge6` or `renter_median_rooms` was
  optimizing against corrupted model moments; those losses/rankings are suspect.
- The young/old **ownership** diagnoses are unaffected (they come from the
  tenure policy/value ranking, not these forward stats).

## 7. Status / next

Working tree: `code/model/intergen_housing_fertility/solver.py` modified
(not committed). Suggested next: recompute the full target-fit table at the
current point with the fix, then proceed to the structural work (young-entry
owner advantage, old-age exit, bequest spec).
