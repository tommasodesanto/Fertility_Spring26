Independent review completed. The corrected main-checkout deliverables are:

- [README.md](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/README.md)
- [observer_contract.json](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json)
- [observer_contract.csv](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.csv)

Ten rows include eight housing/wealth/pooled-PSID restrictions and two approved validation observations. All early estimates and saved uncertainty match their builders. The target design is approved conditionally; this mapping remains **not certified for production SMM**.

Decisive findings for the lead:

1. Corrected the recovered worker's ELDCH description: it is the **oldest** resident own child, not youngest. The current model's recent-parent statistic is any positive dependent count minus never-parent ownership; it cannot reconstruct the exact ACS groups. The 3+ versus 1–2 resident-child room groups also need an explicit child-age/residence approximation.
2. Capped renter room means must use the full income-resolved policy before income collapse. A separate initial observer is needed; the existing historical observer is hard-pinned to old 2023 twelve-row measurement.
3. Current age masks are not literal empirical ranges: prime rows cover model cells 30–57, old wealth 78–85, young ownership 26–37. The contract gives implementable overlap weights, clearly labeled approximation rather than adoption.
4. Aggregate gross earnings formula and 18–65 cell coverage agree with the intended PSID concept. Old wealth uses model pension income versus PSID total family income; $1,000 cutoff and within-wave stock/flow dates remain unverified.
5. Initial PSID birth response should reuse the begin/finish dated branch with the stationary policy at both dates, retaining destination continuation births. The legacy sequential one-period helper suppresses them. Keep the pooled 0.7202462624 (SE 0.0852600513), its uncapped rooms and mapping caveats unchanged.
6. Annual bequests/wealth 0.0088 remains an external restriction with inherited synthetic scale 0.00044, never an empirical SE.

Verification passed: five ACS points and 1,000 paired bootstrap draws/SE/covariance reconstructed from saved metro components; aggregate PSID ratio from yearly totals; all early saved wealth bootstrap SEs and p90/median identity; SA coefficient and covariance arithmetic; actual pure age-index/child-bin functions; JSON/CSV equality; 34 recomputed source/evidence hashes. Large raw hashes were inherited and labeled, with file-stat checks; old quantiles were not reestimated. No model runs, raw-data passes, regression reruns, source edits, commits or pushes. The original recovered worker artifact remains untouched in the isolated checkout.
