# Tenure-barrier probe — lead readout (2026-10-02)

Experimental mechanism test. Not a calibration and not adopted. Fixed price and H0 at the
quarter-saving overnight point. The B1 control reproduces the per-child-need probe bit for bit.
All 12 solves ran locally, because Torch SSH returned Permission denied at about 12:40. Full rows
are in `results.csv` and the receipt is `RECEIPT.md`.

## Question
Why don't capped renter parents buy when the down payment is removed?

## Part A (saved arrays, no solves)
Capped renter parents can afford to buy almost every owner size, even at 80% financing. Their
mean wealth is about 1.5 and their mean income state is about 4 of 0–8. Yet 36–60% of them choose
to rent. The down payment is not what stops them. Caveat: the feasibility replication uses the
income-augmented purchase screen. Positive choice probability on cells it marks infeasible is
exactly zero, which is a consistency check, not a proof.

## Part B: the financing effect (100% minus 80%), per-child need throughout except B6
| Arm | Completed fertility | First births | Entrant 18–21 first birth | Ownership 30–55 |
|---|---|---|---|---|
| B1 control | −0.010 | −0.0002 | +0.28 pp | +9.2 pp |
| B2 selling cost 0 | **+0.012** | **+0.0003** | **+1.22 pp** | +11.0 pp |
| B3 finer owner grid | −0.012 | −0.0003 | +0.27 pp | +9.4 pp |
| B4 finer grid + selling cost 0 | **+0.011** | **+0.0002** | **+1.55 pp** | +10.3 pp |
| B5 tenure shock at minimum | +0.004 | +0.0001 | +0.50 pp | +3.4 pp |
| B6 current need + finer grid | −0.013 | −0.0003 | +0.20 pp | +9.8 pp |

## Readout
- The 6% selling cost is the barrier. Without it, the financing effect on births turns positive,
  including for the matched 18–21 entrants. That is the first right-signed response we have
  obtained without changing preferences.
- A finer owner grid does nothing on its own.
- B5 is unreliable: its fixed-price housing residual is large (0.55).
- Mechanism (inferred, not decomposed): with a selling cost, credit pulls young childless
  households into small homes, and a later birth then forces a costly trade-up. Without the cost,
  owning carries no lock-in.
- Magnitudes are small (+0.011 completed fertility). Price is fixed: in general equilibrium the
  price response would dampen the effect.
- The B2 override sets `P.psi` to 0. That parameter also enters the old-age buyer estate floor
  $-(1-\psi)Q$, which is unlikely to matter for births.
- Six percent is a realistic transaction cost. The economic question is whether the model
  overstates lock-in. Candidate sources: four-year periods, the sale cost charged on the whole
  house at every size change, and no anticipatory starter-home choice.
