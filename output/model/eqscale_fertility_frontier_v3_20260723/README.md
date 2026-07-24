# Fertility frontier v3: margin-split taste noise (entry vs continuation), 2026-07-23

Torch array `14683440` (smoke `14683429`), 147/147 completed, max GE residual
`7.0e-04`. Design: {shared kappa 15.91 (control); split kappa_entry 3 /
kappa_cont 0.3; split kappa_entry 10 / kappa_cont 0.3} x 7 psi_child x 7
delta_alpha_jump at the certified E3 winner, L4 literal parity, SSK power
scale, FL x HSV income. Author-approved probe of the gated
`kappa_fert_continuation` split (commit-gated, default bitwise-inert).

Pre-registered question: can margin-specific noise (noisy entry timing,
near-deterministic continuation) decouple the entry and continuation
margins that every v1/v2 configuration tied together?

**Answer: yes — the frontier rotates through the target neighborhood.**
Best v2 joint distance was 0.1233 (CF 1.26 at p0 0.17, 3+ share 0.08).
Best v3 cell (split 10/0.3, psi 0.5, djump 0.087): CF 1.835, childlessness
0.156 (chosen 0.072 / clock 0.084), joint distance **0.0309** — 4x closer —
with 3+ share 0.320, literal PP 1->2 flow 0.64, and mean age at first birth
27.3 arriving near their data values untargeted. The family ownership gap
rises monotonically in delta_alpha_jump (0.11 -> 0.40 across the row) while
fertility barely moves: the tilt and the fertility block are nearly
separable levers. All best cells sit at the psi grid MINIMUM (0.5) with the
larger entry noise: the optimum lies just outside the scanned box (lower
psi and/or kappa_entry ~ 12-16), where the remaining childlessness gap
(0.156 vs 0.188) plausibly closes.

Mechanism readout: with cheap continuation (kappa_cont 0.3) even psi = 0.5
makes parents progress (3+ share 0.32-0.40); entry noise spreads start
ages (mafb 25.8 -> 27.3 moving kappa_entry 3 -> 10), and late draws meet
the fecundity wall, so childlessness reappears half chosen / half clock.
Compare the shared-kappa control layer: CF capped at ~1.63-1.64 with p0
~0.07 everywhere — the v2 pathology, reproduced under the same externals.

Status: fixed-theta diagnostic cells, not calibrations. Net parameter cost
of the device: one (the split), to be identified by the first-birth timing
distribution (mean age, share 30+) with the parity-progression flows as
overidentification — the A6/A8 data pass supplies these. Next decision
(author's): refinement scan toward (psi < 0.5, kappa_entry 12-16) or fold
the split directly into the reconciled-system recalibration domain.
