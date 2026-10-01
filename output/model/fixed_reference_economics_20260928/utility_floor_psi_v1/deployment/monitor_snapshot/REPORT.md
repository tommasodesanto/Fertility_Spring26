# Provisional free-psi monitoring snapshot

The checkpoint snapshots are timestamped in `torch.json` and `local.json`; later cases are deliberately excluded. There are 37 completed Torch cases, of which 36 have passed computed/admissible GE reports and one is numerically rejected. Torch recorded 45 objective calls. Locally 80 completed cases all passed, with 88 objective calls. These are provisional fast-objective evaluations, not final selected postchecks or optimizer convergence.

Best original-weight loss across **all evaluated points**, including weight-profile chains, is 151.11345889770718 at local chain 1/0009_nm. The best Torch point is 166.81128932452887. The current local housing-weight profile prefers the lower-floor seed's room/ownership fit; its weighted best has original loss 191.31445435875415. Early-fertility weighted best has original loss 151.11345889770718; combined-profile best has original loss 191.31445435875415. The profiles have largely repeated the initial simplex: identical points across profiles do not demonstrate a successful weighting response. Psi remains 0.1355551166583114 at these best points. See `comparison.csv` for weighted/original loss, rooms, ownership, birth-room response and early fertility.

Passed-case median wall time is 121.5–219.2 seconds by Torch chain and 78.8–90.3 seconds by local chain. These are completed-case times from saved receipts; they include fast reporting but exclude deferred standard plots and selected-price repeat. Final native postchecks remain pending. All eight Torch chains were active with no fatal failure receipts.

The rejected case is Torch chain 2/0003_nm (`psi_child=.08`), after 10 lifecycle price evaluations. Every observed renewal residual had the same negative sign; search exhausted both external diagnostic price caps `[.0987336877,6.3189560140]`. The recorded reason is `uncomputed_price_unbracketed`, classified `inadmissible_numerical`, assigned the declared 1e12 objective penalty. It produced no valid equilibrium target fit and is not counted among passed GEs. This establishes a numerical bracket failure within configured caps, not economic infeasibility. The full compact source receipt is `root_rejection.json`; no widening or rerun was performed.

Each reported point has its full 14 original-weight target rows and 31 native parameter rows, plus the experimental weighted target table and exact parameter vector/source report path in `point.json`. No arrays or per-trial plots were collected.

Full tables: [comparison](comparison.csv); [best original fit](best_original_across_all/target_fit.csv) and [parameters](best_original_across_all/parameters.csv); [Torch best fit](best_torch_original/target_fit.csv) and [parameters](best_torch_original/parameters.csv).

| Local profile | Full original fit | Parameters | Experimental weighted fit |
|---|---|---|---|
| Original | [14 targets](base_control/target_fit.csv) | [31 parameters](base_control/parameters.csv) | [Weighted table](base_control/target_fit_weighted.csv) |
| Mean rooms and ownership ×4 | [14 targets](housing_levels4x/target_fit.csv) | [31 parameters](housing_levels4x/parameters.csv) | [Weighted table](housing_levels4x/target_fit_weighted.csv) |
| Early fertility ×4 | [14 targets](early4x/target_fit.csv) | [31 parameters](early4x/parameters.csv) | [Weighted table](early4x/target_fit_weighted.csv) |
| Both sets ×4 | [14 targets](both4x/target_fit.csv) | [31 parameters](both4x/parameters.csv) | [Weighted table](both4x/target_fit_weighted.csv) |
