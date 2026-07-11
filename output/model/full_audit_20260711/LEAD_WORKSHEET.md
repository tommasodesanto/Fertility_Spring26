# Lead worksheet — full audit 2026-07-11 (interim, superseded by FULL_AUDIT_20260711.md)

Facts verified by the lead directly (file:line), before merging subagent output.

## Reproduction (lead-run battery)

| run | config | loss | stored | gap | residual | strict |
|---|---|---|---|---|---|---|
| current_repro_exact | Nb=120, eq10, tol 1e-4 | 14.780020700 | 14.780020700 | 0.000e+00 | 2.906e-05 | yes |
| current_repro_repeat | fresh process, same | 14.780020700 | 14.780020700 | 0.000e+00 | 2.906e-05 | yes |
| housing_repro_exact | Nb=120, eq10 | 13.937130540 | 13.903465319 | +3.367e-02 | 6.597e-05 | yes |
| current_eq30 | max_iter_eq=30 | 15.051377067 | — | +0.2714 | 2.657e-05 | yes |
| current_tol25e6 | eq40, tol 2.5e-5 | 15.051267606 | — | +0.2712 | 3.022e-06 | yes |
| current_pinit_lo | p_init=0.655 | 14.946606375 | — | +0.1666 | 5.415e-05 | yes |
| (battery continuing: pinit_hi, Nb=240 x2, housing repeat/eq30, current Nb240 eq30) | | | | | | |

Key readings:
- Current-bound candidate is bit-reproducible locally (two fresh processes identical to stored to all printed digits). The local dirty tree is evidently code-identical to the torch snapshot for this path.
- Equilibrium-acceptance slack: at the SAME theta, tight equilibrium gives 15.051 (two independent tight configs agree). The production evaluation (max_iter_eq=10, tol_eq=1e-4) reports 14.780. Slack ~0.27 loss units (~1.8%). A price move of 4.6e-05 (0.0066%) moves the loss 0.27, spread over tenure/room moments (owner-renter room gap +0.080 contrib, own_rate_2534 +0.059, renter rooms +0.051, young wealth +0.049, parent-room gap +0.043, own_family_gap +0.029, own_rate -0.021...). Loss differences below ~0.3 between candidates are inside solver-acceptance noise; the third decimal of 14.780 has no meaning; winner-selection across ~24k records partially selects on favorable slack (winner's-curse).
- Housing-relaxed candidate does NOT reproduce bit-exactly locally (gap +0.034, residual 6.6e-5 vs 3.1e-5) — consistent with a platform-sensitive near-tie in deterministic argmax choices; battery2 checks local determinism.

## Spec items verified by lead (line reads)

1. Override merge order: local_panel.py:848-853 {base, **extra(profile+fixed_spec), **income, **theta} — Rouwenhorst income wins over the profile's stale 5-point z-grid. VERIFIED. Wave-2 global stage builds extra the same way (run_wave2_lateral.py:64-76).
2. Rouwenhorst construction local_panel.py:1038-1076: rho_period=rho_annual^4=0.85000; period innovation var = sigma_a^2*sum rho^{2k} (correct aggregation); sigma_log stationary = 0.2310; Kopecky-Suen recursion; binomial weights; E[z]=1 level normalization. VERIFIED construction.
3. MATCHED values are reverse-engineered from the OLD ad hoc grid: rho_annual = 0.85^(1/4) EXACTLY = 0.9601845894; sigma chosen so stationary log sd = old grid's 0.23102 EXACTLY (computed: weights [.1,.2,.4,.2,.1] on log[0.6..1.4] give sd 0.23102). NOT an external calibration. Sommer-Sullivan anchor in compare tool: rho=0.90, sigma=0.20 -> stationary sd 0.459 = 2x the model. Income risk understated by half vs the tool's own literature anchor. FINDING (interpretation/paper-claims tier).
4. Bequests solver.py:5604-5620: normalized form exact; B(0,n)=0; monotone; sigma=1 branch correct; estate tax (rate 0 default) applied to gross estate BEFORE utility; scale theta0*max(1+theta_n*nk,0); nk = encoded child count (get_completed_fertility). Call sites 2070/2438: estate = b_grid + p*H (GROSS of transaction cost psi; heq=(1-psi)pH used only for sale/continuation paths). Evaluated exactly off-grid (no interp). VERIFIED.
   - normalize_bequest_utility DEFAULT False (5611 getattr) and NOT in production_profile_overrides; set only by run_local_polish arg (local_panel.py:578) and the drivers' fixed_spec. Staleness trap for any tool re-solving a record without it.
5. Owner services kernels.py:915-918: ht = hsv - h_bar; clip 1e-10; ht *= chi => chi*(H - hbar) residual form. Matches adopted spec (July 9 F2 decision: paper adopts code). VERIFIED code-side.
6. Finance: parameters.py:80-82 defaults q=(1.04)^4-1=0.16986, delta=1-0.98^4=0.07763, tau_H=0.04 (i.e., OLD spec r=4%, delta=2%). Combined overrides q=0.08243, delta=0.04328 land via apply_overrides; user_cost_rate and R_gross RECOMPUTED after overrides (parameters.py:332-333). VERIFIED. Stale-default trap: tools omitting fixed_spec run r=4%/delta=2% — economically large.
7. Supply: solver.py:1471 H_s = H0*(owner_user_cost/r_bar)^xi_supply; eta_supply override aliases to xi_supply (parameters.py:282-284). VERIFIED eta=1.75 lands. r_bar default 0.16.
   - H0 BOUND FINDING: at p_eq=0.6898, user_cost_rate=0.16571, supply = 9.9997*(0.6898*0.16571/0.16)^1.75 = 5.55 = the delivered mean-rooms moment vs target 5.78. H0 pinned at upper bound 10; the [1,10] box mechanically caps the rooms target the parameter exists to hit (needed ~10.4 at this price). Bound mis-set.
8. Timing: solver.py:3602-3679: g_current = realize_current_cross_section(...) when flag True (getattr default True — safe direction); moments AND fast-path clearing stats computed on g_current; mass conservation asserted atol 1e-10 (3631-3633); wealth moments = beginning-of-period b conditioned on current tenure (3679, deliberate convention). Legacy path only when flag explicitly False.
9. Income profile parameters.py:452-471: entry ages get age_values[0]=0.650 (M1 fixed); step profile; normalized to working-age mean 1; pension balanced PAYGO tau*avg_income*(J_R/(J-J_R)) flat in z; income = (1-tau_pay)*w*profile (after-tax).
10. Age windows solver.py:4614-4616 round((age-18)/4): own_rate_2534 -> model ages 26,30,34 (years 26-37 vs data 25-34); 65-75 -> ages 66,70,74 (years 66-77: still drops 65, includes 75-77). Unavoidable on a 4y grid; needs paper disclosure; NOT a code bug.
11. Loss formula: contributions = weight*(model-target)^2 raw units; verified against stored contributions (160*0.18356^2=5.3914 etc.); sum = stored loss for both candidates.
12. Selection integrity: record_selection_loss = inf unless strict (local_panel.py:1151-1155); reducer filters strict, min loss per arm and cross-arm. No non-strict record can win. Wave-2 raises RuntimeError if the global-stage min-loss record were non-strict (liveness, not correctness, risk).
13. Bounds patching: run_overnight_combined.py:116 and run_wave2_lateral.py:52 monkey-patch lp.GLOBAL_DE_BOUNDS with arm bounds. bounds_for_arm('baseline_pattern') returns PRODUCTION_SEARCH_BOUNDS unchanged (replacements={} for that arm). For relaxed arms, meta['source_controlled_bounds'] is mislabeled (reports patched bounds) while meta['production_profile_spec']['bounds'] reports the source box — contradictory metadata inside relaxed-arm tasks. The extracted report tables list the PATCHED (actual) bounds — correct.
14. Both winning candidates share beta=0.8411423842686492 (beta_annual=0.9577, INTERIOR — first time the 0.94 floor is slack; bequests active theta0=0.132 interior), c_bar_n, theta0 to 16 digits; H0 differs by exactly 1.44: both are pattern-polish descendants of the SAME wave-1 seed; wave-2 LHS never escaped that basin in either arm. Effective search breadth around the reported optimum is ONE basin.
15. At-bound parameters (current arm): c_bar_0=1.28 (upper, exact), chi=1.15 (upper, exact), h_bar_0=1.0 (lower, exact), tenure_choice_kappa=0 (lower, user-preferred), H0=9.9997 (~upper). Preference-relaxed arm results (c_bar_0 cap 1.80) NOT extracted locally — blocked question for cluster.

## Relaxed-housing decomposition (from report summary.json)

Improves: owner-renter room gap (2.663->0.640), renter rooms (1.358->0.331), young wealth (1.590->0.757), own_rate_2534 (0.657->0.003, hits 0.347 vs 0.341), parent 3+ room gap, owner ge6 share, family gap, aggregate rooms. Worsens: old liquid wealth median 1.029->0.379 vs target 2.230 (contrib 1.156->2.743), old own 0.948->0.963 (5.39->6.32), own_rate 0.653->0.697 (0.605->1.467), tfr 1.983->2.097 (0.085->0.640), housing_increment overshoots, childless worsens. Net -0.877.
Reading: relaxing chi (1.15->1.20) + h_bar_0 (1.0->0.92) buys the room/young-tenure block at the price of near-universal ownership and liquid-broke elderly. The scalar win is a weight artifact (old-wealth median weight 0.8 vs room-gap 12). Informative about the chi-bound shadow price; misleading as a candidate. DO NOT ADOPT.

## Blocked pending Torch SSH (kinit/login refresh)

- Local tree vs snapshot diff (mitigated: bit-exact reproduction implies code identity for the evaluated path).
- Full cases.jsonl forensics (23,830 strict records; duplicate-evaluation and seed-distinctness checks at record level).
- lateral_summary.json per task (global_completed counts — LHS truncation check).
- cross_arm_summary.json (arm bests incl. the UNEXTRACTED preference-relaxed arm best).
- Wave-1 seed root (stage3_onehour_from_17_5912) records.
