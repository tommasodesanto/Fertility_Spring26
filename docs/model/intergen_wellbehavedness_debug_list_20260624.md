# Intergen Well-Behavedness Debug List

Created after the June 24 external audit discussion. This is a working checklist,
not a calibration plan.

## Active Issues

1. First-look policy plots must show composed policies, not raw branch policies.
   For a renter who buys, the chosen owner branch must be evaluated at
   transaction-adjusted liquid wealth \(b-pH\), not at original \(b\).
2. Total-wealth policy axes must account for purchase cost. For a renter buying
   \(H\), post-choice liquidation wealth is
   \(b-pH+(1-\psi)pH=b-\psi pH\), not \(b+(1-\psi)pH\).
3. Non-fertile ages should not be plotted as zero fertility probabilities.
   Fertility probabilities are meaningful only in the active fertility window.
4. The PTI screen currently uses maximum financed debt, roughly
   \((q\phi+\tau_H)pH\), rather than actual mortgage debt after the transaction.
   This may mechanically exclude upper owner rungs even for cash-rich households.
5. Owner ladder support is mechanically thin. Given current family-space floors,
   \(\chi\), and PTI, 2-room and 10-room rungs are likely unavailable or
   unattractive by construction.
6. The bequest utility level effect needs audit. With \(\sigma=2\), multiplying
   negative CRRA bequest utility by child-dependent weight can penalize children
   at fixed estate and generate sharp wealth-fertility thresholds.
7. Current fertility architecture is one-shot completed-family-size choice, not
   sequential parity progression. Any 1-to-2 interpretation must be diagnostic
   only unless the state space is changed.
8. Child aging is geometric with one dependent stage. This can create early and
   late reversals in family-space needs and should be checked against intended
   lifecycle mechanics.
9. The liquid-wealth grid is much wider than the occupied support. This is not
   automatically an economic pathology, but it wastes resolution near relevant
   thresholds.
10. Global savings optimization needs an audit because the continuation value is
    non-concave across tenure, fertility, and ladder margins; golden-section
    search assumes unimodality within each conditional branch.

## Immediate Order

1. Fix and inspect composed first-look policy plots.
2. Add accounting/invariant diagnostics for budget residuals, probability sums,
   mass conservation, value monotonicity, and branch feasibility.
3. Run partial-equilibrium ablations for PTI, bequest normalization, owner
   service scaling, renter cap, and dense owner ladder.
4. Only then return to calibration target design.
