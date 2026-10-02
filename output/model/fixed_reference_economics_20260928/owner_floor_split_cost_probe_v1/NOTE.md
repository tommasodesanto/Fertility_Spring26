# Owner-ladder floor probe — lead readout (2026-10-02)

Experimental. Fixed price and H0 at the quarter overnight point; no recalibration. The worker
stopped without a RECEIPT. The table below was extracted by the lead from `case_results.json`.

## Part A: remove the bottom of the owner ladder (current child need, need0)
| Owner sizes | Ownership 30–55 at 80% (whole-population rate) | Financing effect on completed fertility | Financing effect on entrant first births | Loss at 80% (fixed coordinates) |
|---|---|---|---|---|
| 2,4,6,8,10 (current) | 0.640 | −0.0103 | +0.23 pp | 51.6 |
| 4,6,8,10 | 0.615 | −0.0073 | +0.27 pp | 57.5 |
| 6,8,10 | 0.550 | +0.0009 | +0.43 pp | 284.3 |

The need1 arms give the same pattern. Cutting the bottom rungs removes the negative lock-in
effect but does not create a positive credit effect. It also costs ownership and fit.

## Part B: split transaction cost — FAILED, no results
The sandbox engine copy with a buyer cost failed its own identity gate: with psi_buy = 0 it should
reproduce the control, but 17.5% of policy index entries differ. The modification is wrong
somewhere, and no split-cost numbers exist. The DUE lookup was not reported.
