# Revision 1 note (October 3, 2026, 07:06–07:25 New York)

Same Fable 5.1 session, one correction pass on `ECONOMIC_MEMO.md`. Zero solves;
no numbers, figures or code changed. `ERRATA.md` applied. Pre-revision files are
in this folder with `pre_revision_manifest.json`.

## Corrected claims

| # | Pre-revision claim | Correction | Verification |
|---|---|---|---|
| 1 | Early fertility "not a credit phenomenon at all"; "the null is structural"; moment "not identified by the financial block" | Replaced with local insensitivity at the examined states and executed comparisons; indirect channels (budget, housing menu, transaction cost, wealth, price) stated as open; causal test defined as a current-soft credit-policy counterfactual with solved price and standard packet | Inference correction; no source needed |
| 2 | 0.676 "ceiling" as impossibility; mean-age row "forbids" earlier births | Now a conditional bound given the current 18–21 first-birth flow and 22–25 hazard (`analysis/native_grid_analysis.py` lines 123–131); earlier timing could be offset by later births; no impossibility proven; age definition \([25,26)\) with post-birth weight 0.875 per `ERRATA.md` | Script lines re-read; ERRATA.md |
| 3 | 95–97% screen pass read as "binds nobody" | Populations and denominators stated for every screen: closing screen, ending-floor feasibility (14–24% fail), cash test (85–91% fail), owners at floor (4–9% of young realized owners); pass rates not read as effects; near-bound \(\kappa\) flag described as interval-width screen, not proof of deterministic choice | `analysis/out/native_grid_analysis.json` |
| 4 | Hard/quarter rules "tighten the ending floor, not the closing test"; no binding closing test ever run | Hard rule \(A\ge(1-\phi)Q\), quarter \(A+0.25S\ge(1-\phi)Q\); historical kernel rejects renter buyers with `bg_b < dpn` | `purchase_rules_overnight_v1/local_runtime/local_plan.json` lines 148–150, 3075–3077; `purchase_rules_overnight_v1/engines/quarter/refactor_lab/engine/kernels.py` lines 314–328 (rejection) and 912–918 (quarter floor add-on); read this session |
| 5a | "temporary tax-shock paths" | Temporary 100% financing (\(\phi\)) paths, one-date change, 48/64 dates | `purchase_rules_overnight_v1/mechanism_deployment/README.md` line 4 and `accepted_*_temporary_*.json` |
| 5b | "halving of the wealth target" | 35.6% decrease (6.926584 to 4.458387); entrant distribution unchanged | `ERRATA.md`; `input_snapshot/new_wealth/RESULTS.md` |
| 5c | Grid diagnosis cited as covering the new winners; zero endpoint mass as convergence | Diagnosis was at the older chain 16 / case 0046 point; new winners' grid sensitivity is an unresolved check; endpoint mass is not a certificate | `asset_grid_diagnosis_v1/README.md` |
| 6 | Recommended rebuilding the empirical target on model support; literature cited as evidence | Reframed as a diagnostic beside the unrestricted statistic; target changes need measurement design and identification; literature pointers marked unverified | Inference correction |
| 7 | Hacamo (2021) cited as if read | Marked as from project reference notes, not re-opened | Memory note only |

## Remaining uncertainty

- Whether the early-fertility target is reachable under the full scored system is open; only the conditional bound was computed.
- Numerical sensitivity of chain 15 and chain 13 to asset-grid resolution was not tested.
- The empirical early-fertility income gradient and age-specific first-birth hazard shape were not verified in this session.
- No current-soft credit-policy counterfactual exists; the one-line clean test described in the memo (Section 9, item 3) needs lead authorization.
