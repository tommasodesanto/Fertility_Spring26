# Paper Specification-Consistency Audit (2026-07-11)

Audited object: `latex/intergenerational_housing_fertility_paper_draft.tex` (dirty working tree,
529 insertions / 694 deletions vs HEAD `fe092ca`) against the corrected combined-spec code
(`code/model/intergen_housing_fertility/` + `tmp/overnight_combined_20260710/run_overnight_combined.py`).

## Headline

The rewritten model and calibration sections DO describe the corrected combined-spec code: the
Rouwenhorst income process, the normalized bequest function, gross estates, the finance block,
the eta=1.75 supply closure with estimated H0, the residual owner-service form chi(H-hbar), the
{2,4,6,8,10}/cap-6 menu, and the 15-moment/14-parameter system are all correct in the draft.

The paper's numeric tables are NOT the overnight `current_bound_best` (loss 14.7800). They are an
exact transcription of the deliberately frozen "circulation candidate"
(`output/model/intergen_circulation_candidate_20260711/`, record `pattern_i003_d08_m`, case 96 of
`source/task13_cases.jsonl`, rank loss **16.13191908615311**, residual **2.3904e-05**, Nb=120).
Every paper number I could test — all 14 parameter estimates, all 15 model moments, the residual
2.39e-5, the mechanism numbers 0.178 / 0.130 / 0.307 / 2.37x, mean age at first birth 30.9, the
family-stock table (2.4/16.9/43.5/37.2/80.7/74.8), and both quantitative figures (byte-identical
SHA1 to the frozen packet) — matches that frozen candidate exactly. So the draft is internally
consistent and reproducible against a documented artifact, but it reports a candidate whose loss
is 9.1% worse than the strict overnight best now in hand.

## Paper-vs-code discrepancy table

| # | Severity | Paper location (tex lines) | Paper says | Code / audited object | Verdict |
|---|---|---|---|---|---|
| D1 | MAJOR | 986–1010 (parameter table), 1062–1088 (fit table), 1090–1091, 1159–1214, 1216–1253 | Frozen circulation candidate: loss 16.1319, e.g. `H̄=7.51`, `β_ann=0.949917`, `ψ_child=0.35`, `κ_F=1.0`, `θ0=0.2215`, `θ_n=1.0711`, fit residual 2.39e-5 | Overnight `current_bound_best` (reproduced by lead, loss **14.7800**, residual 2.906e-5, strict): `H0=9.9997` (at upper bound 10), `β_ann=0.9577`, `ψ_child=0.269`, `κ_fert=1.891`, `θ0=0.1318`, `θ_n=0.7679`, `h̄_0=1.0` (lower bound), fit table materially different (old-age own 0.948 vs 0.889; childless 0.171 vs 0.220; rooms 5.551 vs 5.188) | Deliberate freeze (documented in `output/model/intergen_circulation_candidate_20260711/README.md`), but every quantitative claim in the paper is conditional on a dominated candidate. Lead must decide: refresh tables or state the freeze in the draft. |
| D2 | MAJOR | 994–1007, footnote 1049–1051 | Parameter table reports point estimates; footnote only: "Several estimates lie at or near the limits of the search region." Identification prose (1035–1044) discusses `ψ_child`, `κ_F`, `χ_O` as identified estimates | `evidence/free_parameters.csv`: `c_bar_0=1.28` AT upper bound, `psi_child=0.35` AT upper bound, `kappa_fert=1.0` AT lower bound, `chi=1.15` AT upper bound, `tenure_choice_kappa=0` at lower bound (disclosed), `h_bar_0=1.10` near lower bound (2% of range). 4 of 14 estimates sit on undisclosed search bounds | Paper presents bound-pinned values as interior estimates; only κ_T's bound is named. Violates the project's own reporting standard and materially changes the identification interpretation. |
| D3 | SUSPECTED | 1062–1105, 886–917 (computation) | All reported quantitative results; no grid-size (Nb) disclosure anywhere | Frozen candidate solved at Nb=120 only (README: "J=17, Nb=120"); production profile designates Nb=240 verification (`production_profile.py:15–16`, `PRODUCTION_SEARCH_NB=120, PRODUCTION_VERIFY_NB=240`); known project finding: Nb=60→240 loss drift, dp overstated | Paper numbers are search-grid numbers not re-verified at Nb=240. Decisive test: one Nb=240 solve of the frozen theta; compare the 15 moments. |
| D4 | MINOR | 431–449, table row 955–956 | "earnings depend on age"; "age profile" — no values | `production_profile.py:23–24`: breaks {22,26,34,46,58}, values {0.650,0.850,1.000,0.985,0.935}; `parameters.py:466–470` normalizes the working-life profile to mean one | Age-profile values and mean-one normalization undisclosed; income process not reproducible from the paper alone. |
| D5 | MINOR | 527–533, 975–977, 1007 | `H^S(r)=H̄(r/r̄)^η`, "r̄ is the reference rental user cost"; `H̄=7.51` | `solver.py:1647`: `Hs[i] = P.H0[i] * (r[i] / P.r_bar[i]) ** P.xi_supply[i]`; `calibration.py:1006`: `"r_bar": np.array([0.16])` | r̄=0.16 never stated; H̄=7.51 is uninterpretable/unreproducible without it. |
| D6 | MINOR | 886–917 | Computation section: no Nb, no grid point count, no search/verify distinction, no search bounds for the 14 estimated parameters | `production_profile.py:15–16, 39–53` (bounds), overnight driver adds `("H0", 1.0, 10.0)` (`run_overnight_combined.py:87,156`) | Reproducibility omission; also blocks readers from seeing D2. |
| D7 | MINOR | 743–823 (entry section) | Entry block fully specified structurally (κ_E, W̄^E, M, G_E) but no calibrated values reported; policy results (±3.4%/±4.9%, −18.8%) depend on them | `parameters.py:142`: `kappa_entry = P.kappa_loc` (=2.0 default); benchmark runs use `population_closure="normalized"` with fixed inflow `E_total=1/J` (`parameters.py:129–131`); entry margin activated only in the outer policy protocol | Entry-block calibration values undisclosed; policy numbers not reproducible from the paper. Cluster-side policy configs unverifiable locally (blocked). |
| D8 | MINOR | 174–180, 894 vs code | Paper: `q` = bond price, return `q^{-1}-1`, conversion `q=(1+i_ann)^{-4}`; "the bond pays b/q" (671) | Code: `q` = 4-year NET RATE, `P.q=(1.02)^4-1` (spec override, `run_overnight_combined.py:160`), `R_gross=1+q`, `user_cost_rate=q+delta+tau_H` (`parameters.py:154`) | Economically identical (paper q^{-1}-1 = code q = 0.08243); pure symbol collision. Internally consistent; flag only as a replication trap. |
| D9 | MINOR | 1062–1088 | Fit table has no weight column | `calibration.py` target set `candidate_replacement_post_audit_v1` weights range 0.8–160 (e.g. old_age_own_rate 160, own_rate 100) | The fit is weight-driven (old-age ownership alone contributes 2.48 of 16.13); weights materially aid interpretation and are the project reporting standard. |

## Verified items (paper = code)

- Income process: 5-state Rouwenhorst, annual rho 0.960 / sd 0.0645, rho_4y = 0.8500 exactly
  (`local_panel.py:1042–1068`); paper grid `{0.613,0.773,0.974,1.227,1.546}` and weights
  `(0.0625,0.25,0.375,0.25,0.0625)` (tex 908–910) reproduce the code construction to 3 decimals;
  mean-one grid normalization stated (tex 449) and coded (`z_grid /= z_weights @ z_grid`).
- Bequests: tex 503–513 displays exactly `theta0(1+theta_n n)[((theta1+max{W,0})^{1-σ}-theta1^{1-σ})/(1-σ)]`
  = `solver.py:5604–5620` with `normalize_bequest_utility=True` (spec override). Gross estate
  `W^B=b+Ph` (tex 501–502) = `solver.py:2066,2070` (`bequest_utility_vec(b_grid + hv,...)`,
  `hv = p_hat*H_own`, no ψ^s, no tax). "no tax at death" correct: `estate_tax_rate=0.0` default, not overridden.
- Finance: i_ann=2% (tex 942 / spec `q=(1.02)^4-1`), delta_ann=1.1% (tex 944 / spec
  `delta=1-(1-0.011)^4`), tau_p=4×1% (tex 900, 946 / `tau_H=0.04`), phi=0.80 financed share with
  down payment `(1-phi)PH_k` (tex 555–560) = `solver.py:2060` `dp_arr=(1-phi_ncs)*hcost`; debt floor
  `b' ≥ -phi P H_k` (tex 564–566) = `solver.py:2061` `bmo=-phi_ncs*hcost`; sale proceeds
  `(1-ψ^s)Ph`, ψ^s=6% (tex 574–577, 949) = `solver.py:2053`, `P.psi=0.06`.
- Supply: functional form matches `solver.py:1647/5484`; eta=1.75 (tex 961 / spec override);
  H0 estimated in [1,10] against ACS 5.780 target weight 6 (tex 975–977 /
  `run_overnight_combined.py:156–158`, `ROOMS_TARGET=5.779970481941968`).
- Owner services: paper displays the residual form `chi_O[h - h̄(n,s)]` (tex 485–494) = code
  `kernels.py:915–918` `ht_c = hsv - owner_h_bar_scale*hbc; ht_c = owner_service_premium*ht_c`
  (`owner_service_premium = P.chi`, `owner_h_bar_scale = 1.0`). July 9 decision implemented.
- Menu/lifecycle: H_own {2,4,6,8,10} (tex 957 / `production_profile.py:17`); renter cap 6
  (tex 959 / :18); J=17 four-year periods, entry 18, last period 82–85 (tex 889–891); retirement 66,
  fertile window 18–42 (tex 891, `calibration.py:17–20`); child-stage maturity 2/9 (tex 911–912 /
  `stage_durations=18/4=4.5`); `E_0=1/J` (tex 912 / `parameters.py:129`); balanced pension rule
  (tex 912–913 / `parameters.py:424–431`); wealth grid [-12,30], core [-5,7], tail knot 15
  (tex 906–907 / `production_profile.py:28–32`).
- Units crosswalk: footnote tex 600–606 discloses index n∈{0,1,2}, n=1 ≙ 1–2 children, n=2 ≙ 3+,
  two children per index unit, CPS measure = 2E[n]; code `diagnostics.py:368`:
  `"tfr": 2.0 * mean_completed_fertility`. Present and correct.
- 15 moments / 14 parameters: paper fit table (tex 1071–1085) lists exactly the 15 moments of
  `candidate_replacement_post_audit_v1` + `aggregate_mean_occupied_rooms_18_85`; all target values
  match the code target set; parameter count 14 matches the searched vector (13 production bounds +
  H0). All model moments match the frozen candidate `packet/target_fit.csv` to displayed precision.
- Entrant wealth: "wealth-to-income distribution, childless renters 25–35, PSID" (tex 963–964) =
  `calibration.py:33–56` PSID quintile nodes via `external_entry_wealth_overrides()` inside
  `base_overrides` (`calibration.py:1022`).
- Timing: fit and figures produced under the repaired post-decision measurement
  (`use_postdecision_current_distribution=True`, `production_profile.py:166`); paper demand equation
  (tex 851–867) integrates current-choice probabilities over the beginning-of-period composition,
  and figure captions state post-choice integration — consistent.
- Figures: `latex/figures/quant_{lifecycle_equilibrium,decision_rules}_repaired_nb120.png` are
  byte-identical (SHA1) to the frozen-candidate packet copies.
- Mechanism/stock numbers: zeta=0.178, r=0.130, total 0.307, 2.37x (tex 1181–1185) match
  `mechanism/rental_cap_mechanism_summary.json`; family-stock table (tex 1195–1206) matches
  `family_stock/family_stock_holders_model.json`; mean age at first birth 30.9 (tex 1135) matches
  `packet/moments.json` `mean_age_first_birth=30.9065`.

## Required-changes list (for the draft)

1. Decide candidate: either refresh all calibration tables/figures/mechanism/policy numbers to the
   overnight strict best (14.7800) or add one sentence stating the reported candidate is a frozen
   2026-07-11 circulation candidate and a better strict optimum exists (D1).
2. Add bound status to the parameter table (or its notes): c̄_0, ψ_child, χ_O at upper bounds,
   κ_F at the lower bound, h̄_0 within 2% of its lower bound; state the search bounds (D2, D6).
3. State the solution grid (Nb=120 search / Nb=240 verification) and whether the reported numbers
   are verified at Nb=240; run that verification (D3).
4. Disclose the income age-profile breaks/values and the working-life mean-one normalization (D4).
5. Report r̄=0.16 next to H̄, otherwise the supply scale is uninterpretable (D5).
6. Report the entry-block calibration (κ_E, benchmark π̄^E, W̄^E, M) used by the policy
   counterfactuals, or state they are set in the policy protocol (D7).
7. Optional: add a symbol note that the quantitative code's `q` is the four-year net rate
   (q_paper^{-1}-1) (D8); add the weight column to the fit table (D9).

## Blocked / out of scope

- Cluster-side policy-run configs behind the +3.4%/+4.9%/−18.8% numbers (SSH unavailable).
- Whether the benchmark markov path actually enforces the paper's entry/scale equation
  (known "entry-margin plumbing" gotcha) — lead's economics task.
- Measurement-window fidelity of individual moments (65–75 etc.) — target-audit specialist.
