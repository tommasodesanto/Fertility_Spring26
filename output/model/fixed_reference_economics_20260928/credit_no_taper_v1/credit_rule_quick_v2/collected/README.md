# Two renter-debt rules: deadline result

Reference: **2007 stationary reference — block0506, September 28 verified export**. Same estimated parameters, including child-benefit psi; common asset price 0.7898695017462086. No recalibration, fertility normalization, new reference, or completed GE.

| Case | Completed fertility | Validation |
|---|---:|---|
| Frozen reference | 2.0999983368 | Previously authenticated reference/control |
| Rollover without age taper | 2.1026541568 | Household/cohort gates passed; exact repeat uncomputed |
| Strict zero renter saving debt and debt-clearing sales | 2.1008010303 | Raw lifecycle output only; rejected entrant feasibility |

Rollover raises completed fertility by 0.12647% and lowers ownership from 66.8165% to 66.5240% at fixed prices. Its counterfactual target loss is 20.7078939879, compared with the frozen fit loss 19.5813107601. These are diagnostic fit scores, not new estimates. All 14 fit rows are in `ours/target_fit.csv`, all 31 inherited parameter rows in `ours/parameters.csv`, and all 17 retained standard plots in `ours/standard_diagnostics/`. Household/fiscal/estate/operator checks are in `ours/gates.json`; no exact repeat was run.

The strict case rejects two age-18 renter entrant cells at wealth -0.2558139535. Their mass is 4.952233248535159e-6, or approximately 0.00802% of the retained entrant cohort (cohort mass 0.06173345618), against the unchanged 1e-12 feasibility tolerance. The preserved census is in `author/inherited_state_diagnostics/`. No mass, debt, value, or entry distribution was repaired. Its raw fertility number must not be interpreted as a validated counterfactual.

GE v3 job 18838220 stopped before any price solve because the remaining time could not fit its 180-second case allowance plus readout reserve. Both GE cases are uncomputed. The scientific deadline, September 30 01:00:25 UTC, was not extended; all numerical jobs are terminal. PE recovery job 18838216 completed the rollover task in 129 seconds and rejected strict entry in 106 seconds. Core results were communicated before the deadline; subsequent work only collected and verified evidence.

The strict model requires an explicit decision about negative-debt entrants before valid GE can be attempted. Changing entrant assets or granting repayment relief would be an additional economic change. The estate-counterparty and finite-grid caveats persist. The frozen reference remains intact.
