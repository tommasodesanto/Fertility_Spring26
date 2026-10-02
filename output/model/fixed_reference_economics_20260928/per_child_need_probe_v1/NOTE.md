# Per-child housing need probe — lead readout (2026-10-02)

Experimental mechanism test. Not a calibration, not adopted. Fixed price and H0, quarter-saving
overnight point (chain 54, loss 51.556). The reference arm reproduces the saved point bit for bit
(RECEIPT.md). Full rows in `results.csv`.

## Question
Credit does not move births because the parent housing need (2.49 rooms, flat in the number of
children at home) never exceeds the 6-room rental cap. If the need grows by one room per extra
child, does 100% financing start to raise births?

## Answer: no
| Row | Financing effect, flat need | Financing effect, +1 room/child | Difference |
|---|---|---|---|
| First-birth flow | -0.000242 (-0.49%) | -0.000232 | +0.00001 |
| Second-birth flow | -0.000220 | -0.000216 | +0.000004 |
| Third-birth flow | -0.000109 | -0.000104 | +0.000005 |
| Completed fertility | -0.0103 | -0.0100 | +0.0004 |
| Ownership 30-55 | +9.4 pp | +9.2 pp | -0.2 pp |
| Entrant (18-21) first-birth prob | +0.23 pp | +0.28 pp | +0.04 pp |

The per-child need works as intended on housing: renter parents at the cap rise to 81% (two
children) and 93% (three), and first_birth_rooms rises 1.17 -> 1.38. It lowers completed
fertility 2.10 -> 1.93 (no recalibration). Yet zero-down financing moves only 3-6 pp of capped
renter parents out of the cap, and births do not respond.

## What this rules out and what remains
- Ruled out (at this magnitude): "the need does not exceed the cap" as the reason credit is inert.
- New fact: with no down payment at all, most capped renter parents still choose to rent at the
  cap, although owning rung 6 gives 1.09 x (6 - 2.49) = 3.83 service rooms vs 3.51 renting at the
  same per-room cost. Something other than the down payment keeps them renting.
- Candidates, unverified: 6% selling cost plus income risk (100% LTV owners must sell after bad
  shocks), the discrete owner ladder (next rung is 8), location and moving choices, the tenure
  taste shock. Next diagnostic: for capped renter parents at phi = 1, decompose the value gap
  between renting at the cap and buying rung 6 or 8 into these terms.

## Caveats
Fixed price and H0, so variant arms carry nonzero renewal/housing residuals (reported, by
design). One point with a material misfit. +1 room per child is one illustrative magnitude.
