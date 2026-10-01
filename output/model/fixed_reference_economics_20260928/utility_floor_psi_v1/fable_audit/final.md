I have what I need. Writing the audit now.

# Verdict

The free-psi restart is correctly wired (candidate psi reaches the solver, the derived coefficient, the 31-row gate and the final check), and no economic gate was relaxed. It will not deliver credible calibrated numbers by morning for two reasons that are verifiable from receipts, not conjecture. First, the renewal-price closure turns the parenthood floor into a housing-price lever: a one-coordinate +0.1 step in h_P (all nine other coordinates identical) raised mean rooms five times more than the first-birth room differential, because the price fell 4.9% to restore replacement fertility. Second, the search budget is far too small for ten-dimensional Nelder-Mead, and a stale hard-coded price start roughly doubles the cost of every evaluation. Expect every chain to stop on the two-hour reserve near its seed or a single simplex vertex.

# Ranked findings

**1. Verified mistake (operational): stale price anchor on every evaluation.** `run_psi.py:156` hard-codes `selected_price=.40056…` and passes it as the start for every candidate and again for the postcheck at `run_psi.py:161`. Both seeds have roots near 0.67 to 0.68 (`utility_floor_nm_v1/chain_1/results/0020_nm/phase_b_ge/price_search.json` shows six fresh solves: start 0.40, expand to 0.64 and 1.03, then three secant steps). Live cases confirm 6 to 8 lifecycle solves and 79 to 90 s per GE locally; a start within a few percent of the root would need 2 to 3. The fast-objective docstring (`fast_objective.py:76-78`) states the anchor is fixed by design, so this is a conscious choice with a large cost. A seed-specific constant start, or the last accepted root, is a numerical start only; bracket, caps and gates are unchanged.

**2. Verified mistake (planning): the budget cannot run a ten-dimensional search.** Each chain allows 200 calls and 7200 s minus a 900 s reserve (`run_psi.py:21-22,83`). At observed ~85 s per GE, local chains get about 70 evaluations; the initial simplex alone is 11. Torch chains run one thread; their first four lifecycle solves took under 130 s including initialization (`deployment/attempt2/native_start_receipts.json`), so roughly 30 to 45 GEs per chain. Snapshot of local chains after 11 GEs each: best points are the seed itself (chains 8, 12) or one simplex vertex (chain 0: h_P vertex; chain 1: tenure-kappa vertex). No psi move has yet produced a best. Also the normalized simplex steps for the two kappas are 0.0003 and 0.0007 of their [0.02, 50] spans (`run_psi.py:45-58`), so those coordinates are effectively unexplored.

**3. Verified fact: the renewal closure makes price a strong lever on housing moments.** Equations are below. From `chain_0/results/best_so_far.json` versus the identical-except-h_P seed (`utility_floor_nm_v1/chain_2/results/0012_nm/.../target_fit.csv`):

| Δ for h_P +0.1 | value |
|---|---|
| price q | −0.033 (−4.9%) |
| mean_rooms | +0.216 |
| first_birth_rooms | +0.041 |
| ownership_30_55 | +0.012 |
| recent_parent_ownership | +0.009 |
| all fertility moments | ≤ 0.001 |

The static Stone-Geary prediction for the differential is α·Δh_P = 0.073 and for mean rooms about 0.02 (parent share ~0.3); the realized room response is ten times that, so the price channel dominates. Along the 0020_nm search path, renewal moves from +0.16 to −0.20 while q rises 2.5×, an implied fertility-price elasticity near −0.4: any fertility-shifting parameter requires a large price move. Torch receipts show the same coupling for psi: chains at psi 0.3 (and chain 4 at 0.2) were still expanding beyond 0.40×1.6³ ≈ 1.64 at their fourth solve, while psi 0.08 chains bracketed below 0.64. No scored target disciplines q, rent or price-to-income (`plan.json:54-139`). Whether this "distorts" is the author's call; that it couples the blocks is established.

**4. Verified fact: first_birth_rooms versus mean_rooms is a structural tension under this floor.** Utility services are h−h_P for renters (`kernels.py:764-766`) and χ(h−h_P) for owners (`kernels.py:852,861`); the floor applies whenever children are at home (`shared.py:255-290`). Hitting 1.465 statically needs h_P ≈ 2.0 before dynamic leakage (realized pass-through was 56% to 89% at the two seeds), i.e. h_P at or above the 2.3 bound (`plan.json` round2 lanes, `h_P_bounds_provenance`), with mean rooms rising further. The one-period-later measurement is confirmed in `run_e5f_transition_calibration.py:1330-1388` (control held childless by construction; treated may have a continuation birth, 36% did at the seed).

**5. Verified fact: early_fertility is nearly invariant and reweighting it adds a constant.** Model values are 0.535 to 0.553 at four points spanning h_P 1.03 to 1.89, chi, beta, kappas and tenure dispersion. The ×4 profile adds exactly 22.5 to the loss (chains 10, 12). Shares at age 25 are 53/38/9/0 percent for 0/1/2/3+ children; two birth windows (18-22, 22-26) cap the moment, and the CPS target counts births before model entry at 18. Treat as a definitional check, not an identification target.

**6. Verified fact: loss is concentrated in two rows, which the profiles do not touch.** Standardized residuals (√contribution) at the best points:

| row | chain_1 best (151.1) | chain_0 best (175.6) |
|---|---|---|
| first_birth_rooms | −6.0σ | −8.9σ |
| recent_parent_ownership | −2.7σ | −7.7σ |
| mean_rooms | +6.0σ | +4.4σ |
| ownership_30_55 | −4.4σ | +2.3σ |
| cps_childlessness | −4.7σ | +0.4σ |
| early_fertility | −2.6σ | −2.7σ |

The housing_levels4x and early4x profiles reweight rows that are not the dominant misfit, so with ~70 evaluations they diagnose a tradeoff only if the search converges, which it will not. Cross-profile comparison under the original objective is free: every case row stores `base_loss` and `base_residual` (`run_psi.py:100`).

**7. Verified facts on the free-psi wiring.** Binding via `inputs.bind` setattr (`inputs.py:107`); coefficient `(1−κ)·psi` from the candidate (`runner.py:39`); exact 31-row equality gate (`phase_b_pilot.py:33-42`); the "Child benefit drift" gate compares the candidate to itself (`phase_b_pilot.py:78,276`) so it never restricts psi; receipts show psi 0.08/0.2/0.3 and matching coefficients. The only defect is cosmetic: the coefficient's status string still says "fixed benefit" (`phase_b_pilot.py:165-166`).

**8. Fiscal and hidden fixed objects.** PAYGO balances identically (residual 1e-13 at both seeds) because entry, mortality, age masses and labor income are exogenous; no hidden free fiscal object. Consequential fixed objects: rental cap 6 rooms (drives recent_parent_ownership through forced ownership), owner max 10 rooms, 2% rate for saving and mortgages with no spread, LTV 0.8, selling cost 6%, equivalence scale ((2+0.7m)/2)^0.7, theta1, alpha 0.733, sigma 2. H0 only sets N and is inert for every per-household moment.

**9. Entry wealth: no evidence it causes room overshoot.** `chain_0/results/input_contract.json` entry report: raw quintile ratios [−2.22, −0.05, 0.10, 0.35, 3.10] of annual gross income, censored and scaled by λ=0.363 to [0, 0, 0.038, 0.128, 1.127]; 66.7% of entrants at exactly zero; mean 0.26 years of income; clipped mass 0. Aggregate wealth is under target at every point and the top quintile is shrunk 2.75×, which depresses young ownership if anything. Units are wealth over annual gross labor income (`inputs.py:71-73`); whether the empirical bins used that denominator is unverified.

**10. Failure taxonomy.** Attempt1 TypeError (string+list), the missing parent directory, the unused copied variable and the NFS visibility race were all orchestration failures before native solves. The only numerical failures were `uncomputed_price_unbracketed` at perturbed starts, which are cap hits, not infeasibility. Four one-expression patches in one evening is itself a schedule risk.

# Equations, identification, plan

**Closure (per `phase_b_pilot.py:93-131`).** Unknowns q (price) and N (household count). Births B(q), entry E (exogenous, 0.0617 per normalized household), normalized demand D(q), supply S(q)=H0·(u·q/r̄)^ξ with ξ=0.63, H0=6.2935. Equations: B(q)=2.1·E (renewal root, tolerance 1e-6) and N·D(q)=S(q) (N solved trivially). So TFR=2.1 is imposed by q for every candidate; the "initial_normalization" row is not a target. The older normalization solved psi for replacement with q clearing housing; I do not recommend restoring it, only note that under the current closure psi and q are nearly collinear in fertility, so psi is identified only through composition and through housing moments' response to q.

**Utility (per `shared.py`, `kernels.py`, `child_preferences.py:11-33`, `fertility_nested.py:89`, `household.py:1488-1515`).** Per period, with m children at home, s = χ·(h−h_P·1{m≥1}) for owners and s = h−h_P·1{m≥1} for renters, Q = c^α s^(1−α), e(m)=((2+0.7m)/2)^0.7:
u = (Q/e(m))^(1−σ)/(1−σ) + ψ·m^(1−κ), with σ=2, α=0.733. First birth subtracts the one-time utility cost first_birth_fixed_cost. Bequest utility is θ0·(θ1+b)^(1−σ)/(1−σ) (which gating spec is active I did not trace). Birth and tenure choices are logit with scales kappa_fert, kappa_fert_continuation, tenure_choice_kappa.

**Identification map (structural, with receipt evidence where available).**

| parameter | primary moments | note |
|---|---|---|
| h_P | first_birth_rooms, recent_parent_ownership (via rental cap), mean_rooms (via q) | FD above: rooms 5× differential |
| chi | ownership_30_55, mean_rooms | owner services premium |
| tenure_choice_kappa | ownership, recent_parent_ownership, mean_rooms | +0.002 step cut loss 15.8 in chain 1 |
| beta_annual, theta0 | wealth_earnings, bequest_wealth, old_dispersion | theta0 is bequest weight, not a child cost |
| first_birth_fixed_cost | childlessness, mean_age | extensive margin |
| psi | q (strongly), childlessness vs fixed cost | level absorbed by q |
| curvature | exactly_one, 3+ share | affects m≥2 only |
| kappa_fert, kappa_fert_continuation | timing dispersion, exactly_one | near-unexplored by NM |

Count is 10 free versus 10 scored, but early_fertility is near-invariant, psi and q are collinear, and curvature, continuation kappa and theta0 all act on parity progression; local rank deficiency is likely. Feasible zero-call rank check: (a) rebuild the 8-column Jacobian from the probe residuals in `utility_calibration_round1_v1/deployment/attempt2/floor_selected_verified/run/cases.json` (the summary at `identification.json` never stored J because the ninth probe failed); (b) pool all `cases.json` rows across fixed-psi and free-psi chains (same target fingerprint, gated at `run_psi.py:32`) and fit a local linear map from the 10 parameters plus q to the 10 base residuals; report singular values.

**Weights.** Each weight is 1/σ² with implied σ: childlessness 0.0053, exactly_one 0.0061, mean_age 0.085 yr, wealth_earnings 0.36, bequest 0.00044, mean_rooms 0.088, ownership 0.021, first_birth_rooms 0.085, recent_parent_ownership 0.0061, early_fertility 0.10. The ×4 profiles halve σ on their rows. Four validation rows carry zero weight and are retained in every table.

**Minimal corrections (propose only).** Seed-specific or last-root price start at `run_psi.py:156,161`; tolerate 1e-10 instead of exact equality in the postcheck (`run_psi.py:164`) and record the difference; fix the coefficient status label; tighten kappa bounds to interpretable ranges so near_bound flags and simplex steps are meaningful. Exploratory choices, not corrections: the ×4 profiles, psi bounds, h_P upper bound.

**Fastest next plan.**
1. Tonight, zero model calls: the two rank checks above; per-profile best under the original objective from stored `base_loss`; read the 17 verified PNGs in the floor packet for rooms by age and family status; check the early_fertility target construction for pre-18 births and the first_birth_rooms empirical window.
2. Let the running chains finish (no restart). Report each selection with its postcheck status and its original-objective loss.
3. One-dimensional sweeps at the chain_1 seed with a warm price start: h_P in {1.5, 1.9, 2.3} and psi in {0.08, 0.2, 0.3}, plus one fixed-q counterfactual at the h_P vertex. About 7 GEs, roughly 7 to 10 minutes locally with a warm start, stop on any unbracketed root. This isolates the price channel and the floor pass-through directly.
4. Author decision on what disciplines q (score a rent or price-to-income moment, change the absorber, or accept the coupling). Only then one Gauss-Newton round with the measured Jacobian: 10 probes plus 3 proposals per round, about 40 GEs for three rounds.

**Falsification tests.** Price channel: falsified if the fixed-q GE at the h_P vertex reproduces the +0.216 rooms, or if rooms do not co-move with q across pooled cases after controlling for h_P and chi. Floor mechanics: falsified if the parent-minus-childless room gap does not rise about 0.73 per unit h_P. Young over-housing or anticipation: falsified if model rooms of childless households aged 22-34 match ACS/AHS in the PNG profiles (control mean was 6.16 rooms at the seed, above the population mean). Wealth mapping: falsified further if the empirical bin denominator matches annual gross income. Early fertility definitional: falsified if the CPS target excludes births before 18.

**Source identity.** The executed solver is the indexed `small_credit_lab` copy under `small_credit_replication_v1/arms/indexed/source`, hash-identical to `code/model/refactor_lab/engine` at launch (`source_pins.json` pairs, runtime cross-check at `runner.py:200-210`). It is not `intergen_eqscale_seq_optimized`, which the model README calls the production reference. Observers come from `code/model/tools` (hashes recorded in each `observers.json`), bind-mounted from a frozen September 28 snapshot on Torch but read from the live tree locally; `authenticate_frozen` is called, but I did not trace that it compares observer hashes to a pinned manifest. Confirm by diffing `observers.json` source hashes between one local and one Torch case once Torch results sync; none are available locally yet.

**Unresolved.** Memory files are symlinked outside the working directory and unreadable here. The with-A comparison fingerprint (loss 18.13) was not independently reverified. Torch progress beyond the start receipts is unknown. The empirical window of the first_birth_rooms target and the bequest gating spec were not traced.
