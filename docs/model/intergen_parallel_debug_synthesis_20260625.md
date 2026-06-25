# Intergen Parallel Debug Synthesis

Date: 2026-06-25

This note summarizes the parallel diagnostic sprint on the June 2026
one-market intergenerational housing-fertility model. The underlying agent
reports and CSV outputs are under
`output/model/intergen_parallel_audits_20260625/`. These are diagnostics, not
calibration runs and not policy counterfactuals.

## Current Issue List

| Priority | Issue | Updated Read |
|---|---|---|
| `P0` | Young ownership near zero | Confirmed as primarily a value-ranking problem, not a static down-payment feasibility problem. For childless renter-origin states at ages 26--34, owner rungs are usually feasible, but renter branch values dominate. Owner-entry probability is essentially zero at ages 26--30 and only about `0.033` at age 34. |
| `P0` | Old ownership too high | Confirmed as near-absorbing same-rung ownership. Among old owners, about `0.976` stay on the same owner rung next period, `0.002` sell to rent, and `0.004` owner-downsize. About `0.471` of old owners are above the renter cap, but they also rarely owner-downsize. |
| `P0` | Same-renter \(c,h\) dips | Not an optimizer failure and not just coarse-grid noise. The dense \(b'\) audit found near-global policies; the fixed-price grid refinement to `Nb=240` shows dips improve in some slices but persist. \(b'(b)\) remains monotone, so the next object is the renter Bellman/intratemporal kink region. |
| `P0` | Renter large-room moment definition | New measurement warning. In the full cached 7D distribution, prime childless renter weak \(h_R\ge6\) mass is about `0.124` versus target `0.138`, but the implemented Markov moment reports `0.013` because the nonlinear threshold is applied after income-state policy collapse. This must be audited before treating the renter large-room miss as economic. |
| `P1` | Renter cap and room separation | Raising `hR_max` from `6` to `8` or `10` hits the implemented renter \(\ge6\) target but compresses the owner-renter room gap from about `2.38` to about `2.10` and lowers ownership at fixed price. The current room gap is materially support-driven. |
| `P1` | Owner ladder has no clean starter rung | The 2-room owner rung is utility-dead; the 4-room rung is weak/floor-like for family states; the 6-room rung is the first robust owner product. Alternative ladders create some starter-rung mass but barely move young ownership at the saved price and tend to raise old ownership. |
| `P1` | Bequest child-level penalty | Confirmed. With \(\sigma=2\), current \(B(a,n)=-\theta_0(1+\theta_n n)/(\theta_1+a)\), so an extra child is a large negative level shift at low estates. At live theta, \(B(0,1)-B(0,0)=-43.25\). |
| `P1` | Gross terminal estate convention | Terminal bequests use \(b+pH\), not net liquidation value \(b+(1-\psi)pH\). At terminal states, mean gross estate is `4.448` versus net estate `4.130`. |
| `P1` | Wealth grid and feasibility cliffs | Liquid grid is overwide: at mass cutoff `1e-6`, only `34/60` grid points are occupied, with support `[-3,15.333]` versus grid `[-35,100]`. Purchase/down-payment cliffs are reachable: `2,955` occupied threshold cases straddle adjacent feasible/infeasible grid cells. |
| `P2` | Tenure smoothing interpretation | `tenure_choice_kappa>0` is a logit mixture in Bellman and forward distribution, not just an argmax tie-breaker. However, hard argmax nearly clears locally at the saved price, so clearing does not mechanically rely on smoothing. Large smoothing, especially `0.05`, materially changes ownership and rung use. |
| `P2` | Fertility terminology | Code/output are internally consistent: one-shot completed-family-size choice, `tfr = 2 * mean_completed_fertility`, `childless_rate = parity_share_0`. Risk is terminology such as `parity_progression_1to2`, `second birth`, and `housing_increment_1to2`, which can invite a sequential hazard interpretation. |
| `P2` | Documentation drift | `CALIBRATION_STATUS.md` and `code/model/intergen_housing_fertility/README.md` are the right orientation. Parent `code/model/README.md` still reads as if the center-periphery path were the current codebase. |

## To-Do List: Diagnostics And Potential Solutions

| Priority | Diagnostic | Potential Solution |
|---|---|---|
| `1` | Audit the renter Bellman/intratemporal housing block around the persistent kink regions: age 26--34, \(z=1.0\), childless renter, especially \(b\in[0.3,0.8]\) and \(b\in[3.0,3.3]\). Compare raw \(Q_R(b')\), intratemporal \(c,h_R\), cap/floor terms, and continuation value pieces before interpolation/export. | If raw objectives switch across local optima, treat the dips as threshold-saving economics and document them. If raw objectives are smooth but exported policies kink, fix interpolation/export. |
| `2` | Audit the renter large-room moment definition in the Markov-income path. Recompute threshold shares before income-state policy collapse and compare with target construction. | If confirmed, fix the moment extractor so nonlinear room-bin shares are computed on full state-by-income distributions. Re-score current points before changing model primitives. |
| `3` | Young-owner value-gap decomposition by age, wealth, income, and rung. Decompose renter advantage into user-cost/rent pricing, \(\chi\), transaction cost, continuation value, and family-size option value. | Add a separate young-entry force if needed: parental down-payment transfers, rent markup or landlord wedge, tax wedge, capital-gains expectation, or a better starter-owner product. |
| `4` | Old-owner retention decomposition. Turn off or vary sale cost, \(\chi\), terminal gross-estate convention, maintenance/user-cost objects, and old-age shocks one at a time at fixed theta/price. | Add an old-age exit or downsizing mechanism: maintenance cost, health/mobility shock, liquidity/expense risk, accessibility shock, or net-estate bequest convention. |
| `5` | Re-score current and frontier points after the renter threshold-moment audit. | Avoid calibrating to a moment computed after inappropriate income-state collapse. |
| `6` | Bequest formula ablations with a rescaled normalized bequest. Current zero-estate normalization at live theta is a large shock: TFR rises from `1.71` to `2.84`. | Use zero-estate normalization only with reparameterized \(\theta_0,\theta_n,\theta_1\), or use a separable warm-glow plus child motive. Use net liquidation estate unless gross house bequests are deliberate. |
| `7` | Owner-ladder redesign probes, but only after young value decomposition. | Candidate menus: `[4,5,6,8,10]` or `[4.5,5,6,8,10]`; pair with owner-size costs or tenure-specific room utility, because starter rungs alone barely move young ownership. |
| `8` | Renter support redesign after moment audit. | Replace the hard cap with a soft cap, size-specific rental supply/price wedge, or scarce large-rental product so large rentals can exist without collapsing the owner-renter room gap. |
| `9` | Numerical hygiene patch and smoke test. | Replace `1e-10` owner residual-service floors and `-1e10` cliffs with explicit feasibility masks; refit the liquid wealth grid to a conservative support such as a buffer below `-8`, dense over `[-4,18]`, and modest upper tail sentinels. |
| `10` | Tenure-smoothing robustness around the final preferred primitives. | Treat `kappa_t=0.01` as a structural/taste smoothing device if retained; do not call deterministic `tenure_choice` the realized policy without also reporting mixture probabilities. |
| `11` | Fertility documentation cleanup. | Rename or qualify `parity_progression_1to2` as a completed-parity conditional share; call `housing_increment_1to2` a two-plus-versus-one completed-family-size housing response proxy, not a second-birth hazard. |
| `12` | Active-strand documentation cleanup. | Add a top note to `code/model/README.md` pointing current June 2026 intergen work to `code/model/intergen_housing_fertility/` and `code/model/run_intergen_model.py`. |

## Recommended Work Order

1. Fix or verify the renter large-room moment definition before any new
   calibration search. This is the most immediate target-object risk.
2. Audit the renter Bellman/intratemporal kink regions. Grid refinement says
   the same-renter dips persist, so the next question is whether the raw
   intratemporal/renter objective itself has stable local switches.
3. Decompose young owner value gaps and old owner retention in parallel. The
   evidence now says these are distinct mechanism failures: young entry is
   value-ranking, while old ownership is near-absorbing same-rung retention.
4. Only after those objects are clean, consider model changes: young-entry
   transfer or rent wedge, old-age exit/downsize shock, rental large-unit price
   wedge, bequest normalization, and a starter-owner menu.
