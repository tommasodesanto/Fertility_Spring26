# Selected 80% buyer diagnostics (saved policy, no model solves)

These are the best available **freshly postchecked experimental fits** selected on October 2, 2026: hard rule, Torch restart chain 11 (loss 97.01121981277964); quarter-saving rule, local original chain 54 (loss 51.55603604909936). Neither search is certified to have converged. The selection manifest SHA-256 is `19db8dfcc928c4cdea5810f70c4e62a9de65372bd038646d0bc2d49532913879`. Jobs 19022340 and 19022341 passed from the separate buyer v4 source; both `completed.json` files state `no_model_solve=true`.

The young buyer window is model ages 26, 30 and 34, labeled 25–34. Buyers here are renter-to-owner flows, not verified lifetime first-time buyers. The closing measure is an **implied net funding ratio** \(\max(0,Q-A)/Q\), where \(Q\) is the purchase price and \(A\) is net closing cash including proceeds from a sold home when relevant. The model does not separately identify gross mortgage debt and cash, so this is not observed mortgage loan-to-value.

| 80% fitted rule | Young renter-to-owner median | 90th percentile | Share above 80% | Share above 90% |
|---|---:|---:|---:|---:|
| Hard | 0.676803 | 0.783059 | 0 | 0 |
| Quarter-saving | 0.737094 | 0.872426 | 0.329628 | 0.046305 |

The financial-access comparison holds the first-birth origin state and all fitted inputs fixed and changes the financed share from \(\phi=0.8\) to \(1\). A household is counted as excluded if no owner product satisfies the rule's purchase screen, housing floor, transaction support, saving/estate floor and budget with positive consumption surplus. The denominator is the modeled **first-birth flow originating from renters**, not all renters and not households who express a preference to buy.

| 80% fitted rule | Share excluded at 80% but feasible at 100%, renter-origin first births | Same numerator / all first births |
|---|---:|---:|
| Hard | 0.618293 | 0.574459 |
| Quarter-saving | 0.074110 | 0.065114 |

At \(\phi=1\), the matched first-birth origin-renter states all have at least one financially feasible owner product under this map. The map does not evaluate continuation values or desired tenure, and it is **not** a causal first-birth response; the dated policy experiment supplies that response. Both observed birth-branch owner-choice audits find zero mass outside the pointwise financial map.

Each arm's `buyer_net_closing_ratio.json` contains weighted distributions, percentiles, CDF, buyer flow masses and negative-closing-cash checks. `matched_first_birth_financial_access.json` contains exact flow denominators, matched 80%/100% eligibility, and failure diagnostics. `completed.json` authenticates the selected source. Full calibration target-fit and parameter tables are under `../../collection/readout/{hard,quarter}/selected_root/`.
