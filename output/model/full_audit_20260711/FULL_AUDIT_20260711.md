# FULL AUDIT — Intergenerational housing–fertility model, combined specification
**Date:** 2026-07-11 · **Auditor:** lead (Fable) + 11 specialist agents + 46 adversarial verifiers · **Scope:** code, targets, calibration, search, results, paper
**Audited object:** dirty working tree at HEAD `fe092ca` (the tree that produced the July 10 overnight; bit-identity with the Torch snapshot established by exact reproduction) plus the isolated drivers in `tmp/overnight_combined_20260710/`.
**Rule followed:** no production file, target, or result was modified. All audit artifacts live in `output/model/full_audit_20260711/`.

---

## 0. Verdict

**1. Is the corrected code solving the intended model?**
Yes on all five audited specification dimensions — the Rouwenhorst income process (construction and 4-year aggregation exact), the normalized bequest function (term-by-term exact, B(0,n)=0, gross estate b+pH), the finance block (q=(1.02)^4−1, δ₄=1−0.989⁴=0.043279309, recomputed user cost, (1−φ) down payment, −φPH debt floor), the supply closure (H_s=H0·(uc/r̄)^1.75, η alias verified), the housing menu ([2,4,6,8,10], cap 6, J=17, Nb=120), and the July 9 timing repair (post-decision distribution used for BOTH moments and clearing; verified end-to-end by new identity tests at 1e-15 precision, with the legacy flag reproducing exactly the old one-period lag).
**BUT** the intended model, as solved, contains one FATAL numerical artifact: households injected onto Bellman-infeasible entry nodes (20.4% of each cohort at this candidate) receive exactly uniform (1/3,1/3,1/3) fertility probabilities from the softmax over −1e10 sentinels (solver.py:2305-2312; mass by parity at j=0 splits into exact thirds). Up to ~21% of all model births are this mechanical default; TFR net of the j=0 artifact is ~1.57 vs the reported 1.983 (upper-bound correction acknowledged). The fertility-block fit (tfr, childless_rate) and the fertility parameters (kappa_fert, psi_child, c_bar_n, h_bar_n) are estimated against contaminated moments. This is a pre-existing latent defect (the July 2 entry-value gotcha) that the *intended corrections amplified*: the age-18 income fix lowered entry income to 0.65 while c_bar_0 sits at its 1.28 cap, making far more entry nodes infeasible.

**2. Is the 14.780 calibration reproducible?**
Bit-exactly, from the local dirty tree, through the same evaluation path the overnight used (loss gap 0.0, residual matching to printed digits, twice, fresh processes). **Not** from committed source: five load-bearing model files are modified vs HEAD (HEAD would crash on missing kwargs, and HEAD's override merge order would silently run the OLD income process), and every driver is untracked. And the number itself is soft in three measured ways: (i) equilibrium-acceptance slack — tight solving at the same theta gives **15.051** (production config: max_iter_eq=10/tol 1e-4 leaves ±0.27); (ii) grid drift — at the production *verification* grid Nb=240 the candidate scores **16.587** (16.340 tight); the record was reported without the "verify at 240" step the July 6 protocol requires; (iii) the housing candidate shows ±0.034 cluster-vs-Mac platform sensitivity. Loss differences below ≈0.3 at Nb=120 are solver noise; the third decimal of 14.780 has no meaning.

**3. Is it economically usable for circulation?**
Not as a final calibration. Grounds: the FATAL fertility artifact (both candidates, and the paper's frozen 16.13 candidate, share it); the H0=10 bound mechanically caps the new rooms target (needs ≈10.4); the "current-bound" winner is a projection-plus-polish descendant of a wave-1 *diagnostic-arm* seed clipped onto the production box (its c_bar_0/χ/h̄₀ corners are where the clip landed, then polish stayed); weights are ad hoc while an unused July 6 bootstrap-SE file implies inverse-variance weights different by up to ~500×, under which the two candidates' ranking flips; and the Nb=240 verification was skipped. It **is** usable as a disclosed working point for internal iteration.

**4. Is the 13.903 relaxed result informative or misleading?**
Informative as a shadow-price diagnostic, misleading as a candidate. It buys the room-gap/young-ownership block (own_rate_2534 hits: 0.347 vs 0.341) by pushing ownership toward universal (aggregate 0.697, old-age 0.963) and collapsing old liquid wealth to 0.379 vs target 2.230 — a catastrophic miss on a moment weighted 0.8 vs the room gap's 12. The scalar win is a weighting artifact. Same fertility artifact, same drift (15.494 at Nb=240). Keep diagnostic-only, as already labeled.

**5. What must be withheld?**
(i) The paper's policy magnitudes (+3.4%/+4.9% births, +0.66%/−18.80% price) — fixed-theta results at the stale 16.13 candidate, from a battery whose README flags "Baseline validation: FAIL" (immaterial numerically) and whose closure-variant sensitivity (up to 2.4× on the grant price effect) is undisclosed; (ii) any TFR/childlessness *fit* claim until the fertility artifact is quantified/removed — the near-perfect fertility fit is partly artifact composition; (iii) any "global search over the stated bounds" claim; (iv) all Nb=120-only numbers presented as final (report the 240 verification alongside); (v) the 14.78-vs-13.90 ranking (flips under SE weights and is within combined solver noise + grid drift).

---

## 1. Honest loss table (lead reproduction battery, this audit)

| Run | Grid | Equilibrium config | Loss | vs stored | Residual | Strict |
|---|---|---|---:|---:|---:|---|
| current stored (cluster) | 120 | eq10, tol 1e-4 | 14.780021 | — | 2.906e-05 | yes |
| current local ×2 (fresh processes) | 120 | eq10 | 14.780021 | +0.000000 | 2.906e-05 | yes |
| current tight | 120 | eq30 / eq40+tol 2.5e-5 | 15.051 | +0.271 | 2.7e-05 / 3.0e-06 | yes |
| current p_init 0.655 | 120 | eq10 | 14.947 | +0.167 | 5.4e-05 | yes |
| current p_init 0.725 | 120 | eq10 | 14.854 | +0.074 | 1.07e-04 | **no** (683s) |
| **current verification grid** | **240** | eq10 | **16.587** | **+1.807** | 7.4e-05 | yes |
| current verification tight | 240 | eq30 | 16.340 | +1.560 | 3.0e-05 | yes |
| housing stored (cluster) | 120 | eq10 | 13.903465 | — | 3.140e-05 | yes |
| housing local ×2 (deterministic) | 120 | eq10 | 13.937131 | +0.034 | 6.597e-05 | yes |
| housing tight | 120 | eq30 | 14.025 | +0.121 | 2.2e-05 | yes |
| **housing verification grid** | **240** | eq10 | **15.494** | **+1.591** | 6.4e-05 | yes |

Grid-drift decomposition (current, 120→240): own_rate_2534 0.2506→**0.1956** (+1.04 of the +1.81), renter rooms +0.29, young wealth 0.543→0.561 (+0.16), owner≥6 share +0.16, family gap +0.16; own_rate improves 0.653→0.632 (−0.29). At the honest grid the young-ownership miss is 14.5 points, not 9. Equilibrium-slack decomposition: a price move of 4.6e-05 (0.0066%) moves the loss 0.27 across the tenure/room moments — households snap across discrete thresholds; the loss is a function of (theta, solver path), not theta alone.

---

## 2. Findings register (post-adversarial-verification, lead-adjudicated)

All 23 FATAL/MAJOR findings survived two-lens adversarial verification (46/46 CONFIRMED); severities below are the lead's final calls after the verifiers' adjustments. Full details, evidence quotes, and decisive tests: `workflow_findings.json`, `workflow_return.json`, and the per-dimension worksheets.

### FATAL
- **DIAG-1 — Mechanical fertility on infeasible nodes.** solver.py:2305-2312 (softmax), 3710-3716 (entry injection). Intended: Gumbel choice over feasible options. Implemented: uniform (1/3,1/3,1/3) on all-sentinel nodes, populated by 20.4% of each entering cohort; ≈21% of births (gross upper bound) are artifacts; TFR net of j=0 artifact ≈1.57 vs 1.983. Contaminated: tfr, childless_rate, own_family_gap composition, fertility-income gradient, and the four fertility parameters. Affects BOTH overnight candidates and the paper's 16.13 candidate (verified same pathology at the relaxed candidate: invalid mass 1.74%, p_birth-zero share 43.2%). Lead-verified at the cited lines; mass-by-parity exact thirds is the smoking gun. Decisive test (specified, not run — requires authorization since it changes evaluation code): zero out fert rows where V≤−1e9 in a KFE copy and re-extract moments; expect TFR −0.2 to −0.4 and childlessness to rise past target.

### MAJOR — results interpretation
- **DIAG-2 — 20.4% of entrants injected onto infeasible nodes; 1.72% of stationary mass lives on V=−1e10 states** following default kernel policies (young-cohort moments partly reflect non-optimizing behavior). Root cause of DIAG-1; the binding channel is fertility.
- **DIAG-3 — 40.3% of childless fertile-window mass has p_birth = 0 exactly** (birth option budget-infeasible; correct Gumbel treatment of −∞, not a bug): the extensive fertility margin is a feasibility cliff, not a smoothed choice; kappa_fert is identified off the other 60%; housing-cost counterfactuals will jump as mass crosses the cliff.
- **GRID — Nb=120→240 drift +1.8 / +1.6** (this audit, §1): the reported candidates were never verified at 240, violating the July 6 protocol; own_rate_2534 is the dominant drifter.
- **SLACK — equilibrium-acceptance wobble ±0.27** at fixed theta (this audit, §1): with ~23.8k strict records ranked by the noisy evaluator, selection partially rewards favorable slack (winner's curse).
- **S1 — Both winners descend from ONE wave-1 seed** (beta, c_bar_n, theta0 identical to 16 digits; remaining differences decompose exactly into pattern steps, e.g. ΔH0=2·0.08·9=1.44). The wave-2 LHS "global" stage (~1,416 draws) never beat that lineage. Moreover the shared seed almost surely came from wave-1's **housing_relaxed (diagnostic-only) arm** and was silently clipped onto the production box by np.clip in `global_unit_from_theta` — h̄₀=1.0 exactly at the production lower bound while the housing winner sits at h̄₀=0.921<1, representable only in the relaxed box. The current-bound corners (c_bar_0, χ, h̄₀) are clip-then-polish artifacts confirmed as *box-constrained* optima by 362+ polish evals, but "global search over the stated bounds" is unsupported. (Cluster-side seed-path confirmation blocked.)
- **T-DATA-2 — Bootstrap SEs exist (code/data/moment_standard_errors/, July 5-6, reproduces 12/14 targets <5e-7) and were not used.** Implied inverse-variance weights differ from the hand weights by up to ~500× (relative); the old-age wealth-gap moment (SE=68% of target) is near-optimally weighted while own_family_gap is relatively down-weighted ~500×. The current-vs-housing ranking flips under SE weights. No inference can hang on the current loss.
- **T1 — Hard-targeted old-age wealth MEDIAN is grid-discrete** (uninterpolated b-grid value; piecewise-constant in theta). One node move shifts the objective by 0.54–1.24 — more than the 0.88 gap between the candidates. Violates the standing "never hard-target medians" rule (memory: grid-discrete medians break surrogate/BO too).
- **T2/F2 — H0=9.9997 at the [1,10] upper bound; rooms target unreachable inside the box** (needs ≈10.4 at p_eq; supply factor (uc/r̄)^1.75=0.555 with the stale r̄=0.16 pivot inherited from the old finance regime). The one newly *searched* parameter cannot do its stated job. (Materiality verifier argued MINOR on pure loss share — 2.1% — but the lead keeps MAJOR: it defeats the stated identification design of the new target.)
- **T3 — theta_n still effectively unidentified (D11 unresolved).** Disciplined by a single ~0.8%-of-loss moment; the two near-equal candidates differ 0.768 vs 1.109 in theta_n while that moment moves 0.33%; theta1=0.01 makes the bequest motive saturate by W≈0.5 so theta_n acts as a level jump, not a saving slope. theta0 interior does not rescue it.
- **INC-1 — Income-process dispersion is reverse-engineered, not calibrated.** rho_annual=0.85^(1/4) exactly; sigma chosen to hit the old ad hoc grid's stationary log-sd 0.23102 exactly (lead-verified arithmetic; agent-verified to 3e-17). Half the dispersion of the repo's own Sommer–Sullivan anchor (0.46) and 2–3.5× below standard estimates. Disclosed as "matched" in the draft, but the empirical justification is continuity with a grid that never had one. First-order for the five wealth/income-conditioned targets and for risk-sensitive policy results.
- **P2 — 4 of 14 reported estimates sit on undisclosed bounds in the paper's tables** (16.13 candidate: c_bar_0, psi_child, χ upper; kappa_fert lower; h̄₀ within 2%); identification prose treats them as interior.
- **PS-1/P1 — The paper's entire quantitative section is the frozen 16.13 candidate** (`output/model/intergen_circulation_candidate_20260711`, documented and reproduced — a deliberate freeze that predates the overnight best by ~10.5h, not an error). Same targets/weights, so 16.13 vs 14.78 are comparable; several fit narratives flip sign at 14.78 (own_rate under→over; childlessness over→under); the at-bound set changes. Fit prose needs rewriting, not renumbering.
- **PS-2 — Headline policy magnitudes** (abstract/intro/§7/conclusion) come from the phase9b fixed-theta battery at 16.13; the same table's closure variants give price effects differing up to 2.4× (fixed-scale +0.66% vs scaled +1.56% vs published +1.34%), undisclosed; the battery README itself labels them "fixed-theta diagnostics, not recalibrated SMM losses."

### MAJOR — reproducibility/process (no effect on the computed numbers)
- **PROV-1/2/3 — Not reproducible from committed source.** Five modified tracked files are load-bearing (incl. the merge-order flip that lets income beat the profile, and the one-line `.reshape(-1)` H0-override fix at parameters.py:299 without which the scalar H0 search crashes); all drivers untracked; results git-ignored. HEAD would compute a *different economy* (old income process) where it ran at all. Wave-1 metadata marks relaxed arms `diagnostic_only`; wave-2 and the reducer drop that marker (reducer takes an unconditional strict min per arm).
- **BQ-1 / INC-2 / PS-3 — The stale-defaults trap family.** `normalize_bequest_utility` defaults False; parameters.py defaults keep the OLD finance regime (q=0.1699, δ=0.0776 — user cost 1.73× the combined one) and the OLD income process; `production_profile_overrides()` carries NONE of q/δ/η/H0/bequest-flag/Rouwenhorst. Eight checked-in re-solve tools and every CLI subcommand except local-polish silently evaluate the old economy on a combined-spec record; `build_intergen_mechanics_packet.py` applies the combined spec only via an opt-in flag and ignores the flag stored in the record; per-record best.json does not store the overrides. Any post-hoc figure, diagnostic, or policy number is one omitted flag away from the wrong model. (PS-4 — the dt_cp-strand plotters — downgraded MINOR: wrong strand, produced no current paper artifact.)
- **BOUNDS-METADATA — under runtime patching, `meta['source_controlled_bounds']` reports the *patched* box (mislabeled) while `production_profile_spec.bounds` reports the source box** — contradictory metadata inside relaxed-arm tasks; the actual box is correct in `meta['bounds']` and in the report CSVs.
- **T-DATA-1 — The new rooms target is built by an uncommitted builder edit,** regenerated ~3h before the overnight, byte-matches the driver constant, and its own CSV metadata says "candidate; re-audit before SMM use." It is also **not national ACS**: it is the 42-metro MMS matched sample of household heads (18-85, literal ROOMS, HHWT) — consistent with the other room targets but not with a "national" description.

### Selected MINOR / notes (full list in workflow_findings.json — 29 MINOR, 13 SUSPECTED)
- Age windows on the 4-year grid: own_rate_2534 ≡ model ages 26/30/34; "65-75" ≡ 66/70/74 (drops 65, includes 75-77). Unavoidable grid mapping; disclose in the paper (weights 80 and 160 sit on these).
- Two different "childless" definitions coexist across moments (childless-in-household incl. empty-nesters for rooms/young-wealth vs never-parent for the old-age gap) — matches the data cuts per the data audit, but document it.
- Room block vs ownership block use different head/structure universes (PERNUM==1 all structures vs RELATE==1 & UNITSSTR 3:10).
- tfr/childless targets are citation-level (CPS June 2024 H2 cell), no raw artifact in repo.
- Value-function "monotonicity violations" (168) are sentinel contamination on 0.57% of mass, not preference violations.
- Convergence basin asymmetry: from p_init above equilibrium the solver ran 683s and missed strict tolerance.
- 16.13-era ledger entries BQ-2/HF-4 in the untracked issues file are stale (describe pre-combined code).
- panel-path `write_best_case_diagnostics` solves without profile/fixed-spec overrides (another staleness surface).
- F1 (owner services paper-vs-code) — **withdrawn by the lead**: the dirty-tree draft displays χ_O[h−h̄(n,s)] (tex:489), matching kernels.py:915-918. The finance auditor read a stale copy. The verifiers' A/B stands as documentation: the alternative form would shift ownership moments 9-12pp, so the display choice is substantive and now consistent.

### VERIFIED (independently checked and correct — highlights)
Rouwenhorst construction and 4-year aggregation; E[z]=1 normalization; age-18-21 income = 0.650 (M1 fixed) with one income array feeding Bellman/KFE/moments/entry/pension; balanced PAYGO pension; normalized bequest exact to 1e-13 with correct σ=1 guard, gross estate, dormant estate tax ordering, encoded-n scaling; (1+q) symmetric on assets and debt; (1−φ) down payment with −φPH floor on all owner rungs incl. stayers; renter b′≥0; menu/cap in force; rent=(q+δ+τ_H)p; supply identity holds on all repro records; post-decision timing end-to-end (three identity tests at ~1e-15; legacy flag = exact one-period lag = the pre-repair defect, now gated); objective = Σw(m−t)² exactly, both stored losses recomputed to the exact float; strict-selection airtight at all three levels; 13/15 target values byte-reproduced from in-repo builders; both stored candidate records exactly re-solved (§1).

---

## 3. Provenance manifest (summary)

HEAD `fe092ca` (2026-07-10, main; pushed). Working tree dirty: ~45 modified + ~80 untracked paths; load-bearing for the candidate: `intergen_housing_fertility/{local_panel,solver,calibration,parameters,cli}.py` (modified), `production_profile.py` + tests (untracked), `tmp/overnight_combined_20260710/*` (7 files, untracked), `code/model/tools/run_intergen_combined_recalibration.py` + 2 submit scripts (untracked), report extracts (git-ignored). SHA256 of all 37 relevant files: `provenance/hashes.json`. Diffs-vs-HEAD summarized per file in `provenance/PROVENANCE_MANIFEST.md`. Torch snapshot diff **blocked** (SSH); mitigated by bit-exact reproduction of the current-bound record through the identical code path. Repo driver ≠ overnight drivers (caps 160 evals/60 min vs 850/330; no jitter; no arm bounds): the committed driver is a *test* variant, not what ran.

## 4. Target-fit tables

### 4.1 current_bound_best — stored, Nb=120 (loss 14.780020699972585, residual 2.906e-05, p=0.68976)
| Moment | Target | Model | Gap | Weight | Contribution |
|---|---:|---:|---:|---:|---:|
| `old_age_own_rate` | 0.764261 | 0.947826 | +0.183565 | 160.0 | 5.391359 |
| `prime30_55_childless_owner_minus_renter_mean_rooms` | 2.418762 | 1.947705 | -0.471057 | 12.0 | 2.662734 |
| `young_childless_renter_liquid_wealth_to_annual_gross_income_2535` | 0.179226 | 0.543225 | +0.364000 | 12.0 | 1.589948 |
| `prime30_55_childless_renter_mean_rooms` | 3.805288 | 4.281094 | +0.475806 | 6.0 | 1.358349 |
| `old_nonhousing_wealth_to_income_median_6575` | 2.230461 | 1.028539 | -1.201922 | 0.8 | 1.155692 |
| `own_rate_2534` | 0.341166 | 0.250554 | -0.090612 | 80.0 | 0.656843 |
| `own_rate` | 0.575472 | 0.653252 | +0.077780 | 100.0 | 0.604970 |
| `prime30_55_childless_owner_share_rooms_ge6` | 0.596131 | 0.741074 | +0.144942 | 25.0 | 0.525208 |
| `aggregate_mean_occupied_rooms_18_85` | 5.779970 | 5.551108 | -0.228863 | 6.0 | 0.314269 |
| `own_family_gap` | 0.167662 | 0.230799 | +0.063137 | 45.0 | 0.179382 |
| `prime30_55_parent_3plus_minus_1to2_mean_rooms` | 0.367700 | 0.242956 | -0.124744 | 8.0 | 0.124488 |
| `old_parent_childless_nonhousing_wealth_to_income_gap_6575` | 1.007450 | 1.251428 | +0.243978 | 2.0 | 0.119051 |
| `tfr` | 1.918000 | 1.983159 | +0.065159 | 20.0 | 0.084914 |
| `housing_increment_0to1` | 0.664435 | 0.687187 | +0.022752 | 14.0 | 0.007247 |
| `childless_rate` | 0.188000 | 0.171317 | -0.016683 | 20.0 | 0.005566 |

### 4.2 current_bound_best — audit re-solve at Nb=240 (loss 16.587070, residual 7.4e-05)
| Moment | Target | Model | Gap | Weight | Contribution |
|---|---:|---:|---:|---:|---:|
| `old_age_own_rate` | 0.764261 | 0.946603 | +0.182342 | 160.0 | 5.319792 |
| `prime30_55_childless_owner_minus_renter_mean_rooms` | 2.418762 | 1.938778 | -0.479983 | 12.0 | 2.764609 |
| `young_childless_renter_liquid_wealth_to_annual_gross_income_2535` | 0.179226 | 0.561413 | +0.382188 | 12.0 | 1.752810 |
| `own_rate_2534` | 0.341166 | 0.195616 | -0.145550 | 80.0 | 1.694794 |
| `prime30_55_childless_renter_mean_rooms` | 3.805288 | 4.329793 | +0.524505 | 6.0 | 1.650634 |
| `old_nonhousing_wealth_to_income_median_6575` | 2.230461 | 0.955318 | -1.275143 | 0.8 | 1.300791 |
| `prime30_55_childless_owner_share_rooms_ge6` | 0.596131 | 0.762005 | +0.165874 | 25.0 | 0.687851 |
| `aggregate_mean_occupied_rooms_18_85` | 5.779970 | 5.539805 | -0.240166 | 6.0 | 0.346077 |
| `own_family_gap` | 0.167662 | 0.254068 | +0.086406 | 45.0 | 0.335970 |
| `own_rate` | 0.575472 | 0.631728 | +0.056256 | 100.0 | 0.316468 |
| `prime30_55_parent_3plus_minus_1to2_mean_rooms` | 0.367700 | 0.231082 | -0.136618 | 8.0 | 0.149316 |
| `old_parent_childless_nonhousing_wealth_to_income_gap_6575` | 1.007450 | 1.257837 | +0.250387 | 2.0 | 0.125388 |
| `tfr` | 1.918000 | 1.993138 | +0.075138 | 20.0 | 0.112916 |
| `housing_increment_0to1` | 0.664435 | 0.698650 | +0.034215 | 14.0 | 0.016390 |
| `childless_rate` | 0.188000 | 0.162247 | -0.025753 | 20.0 | 0.013265 |

### 4.3 housing_relaxed_best — stored, Nb=120 (loss 13.903465318628617, residual 3.140e-05) — DIAGNOSTIC ONLY
| Moment | Target | Model | Gap | Weight | Contribution |
|---|---:|---:|---:|---:|---:|
| `old_age_own_rate` | 0.764261 | 0.963057 | +0.198796 | 160.0 | 6.323147 |
| `old_nonhousing_wealth_to_income_median_6575` | 2.230461 | 0.378936 | -1.851525 | 0.8 | 2.742517 |
| `own_rate` | 0.575472 | 0.696577 | +0.121104 | 100.0 | 1.466621 |
| `young_childless_renter_liquid_wealth_to_annual_gross_income_2535` | 0.179226 | 0.430328 | +0.251102 | 12.0 | 0.756627 |
| `prime30_55_childless_owner_minus_renter_mean_rooms` | 2.418762 | 2.187736 | -0.231026 | 12.0 | 0.640477 |
| `tfr` | 1.918000 | 2.096905 | +0.178905 | 20.0 | 0.640141 |
| `prime30_55_childless_owner_share_rooms_ge6` | 0.596131 | 0.724145 | +0.128014 | 25.0 | 0.409691 |
| `prime30_55_childless_renter_mean_rooms` | 3.805288 | 4.040040 | +0.234752 | 6.0 | 0.330652 |
| `aggregate_mean_occupied_rooms_18_85` | 5.779970 | 5.585081 | -0.194889 | 6.0 | 0.227891 |
| `housing_increment_0to1` | 0.664435 | 0.759333 | +0.094898 | 14.0 | 0.126080 |
| `old_parent_childless_nonhousing_wealth_to_income_gap_6575` | 1.007450 | 1.247288 | +0.239839 | 2.0 | 0.115045 |
| `own_family_gap` | 0.167662 | 0.212293 | +0.044632 | 45.0 | 0.089640 |
| `childless_rate` | 0.188000 | 0.147860 | -0.040140 | 20.0 | 0.032225 |
| `own_rate_2534` | 0.341166 | 0.346889 | +0.005723 | 80.0 | 0.002620 |
| `prime30_55_parent_3plus_minus_1to2_mean_rooms` | 0.367700 | 0.371080 | +0.003380 | 8.0 | 0.000091 |

### 4.4 housing_relaxed_best — audit re-solve at Nb=240 (loss 15.494040, residual 6.4e-05)
| Moment | Target | Model | Gap | Weight | Contribution |
|---|---:|---:|---:|---:|---:|
| `old_age_own_rate` | 0.764261 | 0.964347 | +0.200086 | 160.0 | 6.405478 |
| `old_nonhousing_wealth_to_income_median_6575` | 2.230461 | 0.148007 | -2.082454 | 0.8 | 3.469291 |
| `own_rate` | 0.575472 | 0.699565 | +0.124093 | 100.0 | 1.539904 |
| `prime30_55_childless_owner_minus_renter_mean_rooms` | 2.418762 | 2.110925 | -0.307837 | 12.0 | 1.137163 |
| `tfr` | 1.918000 | 2.096983 | +0.178983 | 20.0 | 0.640697 |
| `young_childless_renter_liquid_wealth_to_annual_gross_income_2535` | 0.179226 | 0.405975 | +0.226749 | 12.0 | 0.616983 |
| `prime30_55_childless_renter_mean_rooms` | 3.805288 | 4.083898 | +0.278610 | 6.0 | 0.465742 |
| `prime30_55_childless_owner_share_rooms_ge6` | 0.596131 | 0.712459 | +0.116328 | 25.0 | 0.338307 |
| `aggregate_mean_occupied_rooms_18_85` | 5.779970 | 5.574401 | -0.205570 | 6.0 | 0.253554 |
| `housing_increment_0to1` | 0.664435 | 0.796184 | +0.131749 | 14.0 | 0.243010 |
| `old_parent_childless_nonhousing_wealth_to_income_gap_6575` | 1.007450 | 1.320624 | +0.313174 | 2.0 | 0.196156 |
| `own_family_gap` | 0.167662 | 0.222233 | +0.054571 | 45.0 | 0.134011 |
| `childless_rate` | 0.188000 | 0.142185 | -0.045815 | 20.0 | 0.041981 |
| `prime30_55_parent_3plus_minus_1to2_mean_rooms` | 0.367700 | 0.330283 | -0.037416 | 8.0 | 0.011200 |
| `own_rate_2534` | 0.341166 | 0.343817 | +0.002651 | 80.0 | 0.000562 |

## 5. Parameter / bound tables

beta is stored as the 4-year factor; beta_annual = beta^(1/4) = **0.957673** for both candidates (interior; the ≥0.94 external floor is slack for the first time — bequests active).

### 5.1 current_bound_best (bounds = production box + H0∈[1,10])
| Parameter | Lower | Estimate | Upper | Flag |
|---|---:|---:|---:|---|
| `H0` | 1.0000 | 9.999672 | 10.0000 | near bound (2%) |
| `alpha_cons` | 0.4000 | 0.673106 | 0.9500 |  |
| `beta` | 0.9400 | 0.841142 | 0.9950 |  |
| `c_bar_0` | 0.0800 | 1.280000 | 1.2800 | AT upper |
| `c_bar_n` | 0.0500 | 0.455979 | 1.5000 |  |
| `chi` | 0.4000 | 1.150000 | 1.1500 | AT upper |
| `h_bar_0` | 1.0000 | 1.000000 | 6.0000 | AT lower |
| `h_bar_jump` | 0.0500 | 1.475894 | 2.5000 |  |
| `h_bar_n` | 0.0200 | 0.984107 | 2.0000 |  |
| `kappa_fert` | 1.0000 | 1.890630 | 12.0000 |  |
| `psi_child` | 0.0000 | 0.269080 | 0.3500 |  |
| `tenure_choice_kappa` | 0.0000 | 0.000000 | 0.1200 | AT lower |
| `theta0` | 0.0000 | 0.131830 | 2.0000 |  |
| `theta_n` | 0.0000 | 0.767920 | 1.5000 |  |

### 5.2 housing_relaxed_best (bounds = housing-relaxed box: χ≤1.60, h̄₀≥0.25, h̄_jump≤3.50, h̄_n≤3.00)
| Parameter | Lower | Estimate | Upper | Flag |
|---|---:|---:|---:|---|
| `H0` | 1.0000 | 8.559672 | 10.0000 |  |
| `alpha_cons` | 0.4000 | 0.629106 | 0.9500 |  |
| `beta` | 0.9400 | 0.841142 | 0.9950 |  |
| `c_bar_0` | 0.0800 | 1.280000 | 1.2800 | AT upper |
| `c_bar_n` | 0.0500 | 0.455979 | 1.5000 |  |
| `chi` | 0.4000 | 1.200617 | 1.6000 |  |
| `h_bar_0` | 0.2500 | 0.921047 | 6.0000 |  |
| `h_bar_jump` | 0.0500 | 1.662749 | 3.5000 |  |
| `h_bar_n` | 0.0200 | 0.825707 | 3.0000 |  |
| `kappa_fert` | 1.0000 | 2.004678 | 12.0000 |  |
| `psi_child` | 0.0000 | 0.287000 | 0.3500 |  |
| `tenure_choice_kappa` | 0.0000 | 0.000000 | 0.1200 | AT lower |
| `theta0` | 0.0000 | 0.131830 | 2.0000 |  |
| `theta_n` | 0.0000 | 1.109009 | 1.5000 |  |

At-bound reading (lead, §2 and ECONOMIC_ASSESSMENT.md): H0 upper = mis-set box (mechanical); c_bar_0 upper = young-wealth damper at its cap (weak identification + missing mechanism; preference-relaxed arm result not extracted — blocked); χ upper & h̄₀ lower = the room-gap/tenure conflict (one parameter, two jobs); tenure_choice_kappa lower = user-imposed external restriction that should leave the search vector.

## 6. Code-versus-mathematics specification table

| Object | Intended | Implemented | Where | Status |
|---|---|---|---|---|
| Income process | 5-state Rouwenhorst, ρ_a=0.9601845894, σ_a=0.0645373326, 4-year periods, ρ₄≈0.85 | exact; σ_ε,₄²=σ_a²Σρ^2k; ψ=σ_log√(N−1); binomial weights; E[z]=1 | local_panel.py:1038-1076 | VERIFIED (provenance of values: INC-1) |
| Income profile | ages 18-21 not at peak | 0.650 at entry; step profile; working-mean normalized; single array everywhere | parameters.py:452-471, 408-421 | VERIFIED |
| Pension | balanced, consistent | τ·avg_income·(J_R/(J−J_R)), flat in z, Bellman=KFE | parameters.py:424-437 | VERIFIED |
| Bequest | θ₀(1+θ_n n)[((θ₁+W)^{1−σ}−θ₁^{1−σ})/(1−σ)], B(0,n)=0, W=b+pH, tax-then-utility | exact (1e-13); max{W,0} clip; encoded n; dormant tax correctly ordered | solver.py:5604-5620, 2063-2070 | VERIFIED (flag default: BQ-1) |
| σ=1 limit | log form | correct unreachable guard (σ=2 fixed) | solver.py:5612-5616 | VERIFIED |
| q | (1.02)⁴−1=0.08243216 | override lands; R_gross, user_cost recomputed; symmetric on b≷0 | parameters.py:80,332-333 | VERIFIED (defaults trap: BQ-1 family) |
| δ | 1−0.989⁴=0.043279309 | override lands; used in owner flow cost and rent | parameters.py:81,332 | VERIFIED |
| Down payment | (1−φ)pH, φ=0.80 financed share | exact; debt floor −φpH all rungs incl. stayers; no path uses φ as DP share | kernels/solver (audit C) | VERIFIED |
| Supply | η=1.75, H0 searched to hit E[H]=5.77997 | H_s=H0(uc/r̄)^η; η alias OK; **H0 at bound; r̄=0.16 stale pivot** | solver.py:1471; parameters.py:146,282-284 | VERIFIED code / MAJOR design (T2/F2) |
| Owner services | χ(H−h̄) (adopted July 9) | χ·max(H−h̄,ε) | kernels.py:915-918; paper tex:489 agrees | VERIFIED (F1 withdrawn) |
| Menu/grid | [2,4,6,8,10], cap 6, Nb=120, J=17 | in force in the overnight chain and both records | production_profile.py; audit C | VERIFIED |
| Timing | post-decision moments+clearing, horizon-0 incl. same-period birth-purchases | g_current for both; identity tests 1e-15; horizon-0 = identity | solver.py:3602-3679; TIMING_AUDIT.md | VERIFIED |
| Objective | Σ w(m−t)², 15 moments | exact; stored losses recomputed to the float | calibration.py; recompute_loss.py | VERIFIED |
| Fertility choice | Gumbel over feasible options | uniform on all-infeasible nodes, populated by entrants | solver.py:2305-2312, 3710-3716 | **FATAL (DIAG-1)** |

## 7. Target provenance and measurement

Full table: `target_audit/TARGET_DATA_AUDIT.md` + `target_audit/TARGET_MODEL_AUDIT.md`. Summary: 13/15 values byte-reproducible in-repo (bootstrap harness reproduces 12/14 to <5e-7); tfr/childless are CPS-2024 citation-level; the new rooms target is MMS-42-metro (not national), from an uncommitted builder edit, self-flagged "candidate". Model counterparts verified line-level for all 15 (age windows, conditioning, annual-income denominators fixed as of the overnight code — the target-object audit doc is stale on this). Known measurement conventions to disclose: 4-year age-window mapping; two childless definitions; head/structure universe differences; wealth-moment timing = beginning-of-period b conditioned on current tenure.

## 8. Numerical diagnostics packet

`diagnostics/current_bound/` (standard package plot set via diagnostics.write_diagnostics; no new graph types) + `diagnostics/DIAGNOSTICS_REPORT.md`. Clean: mass accounting exact (total 1.0, per-age 1/17 to 5e-16), no b-grid edge bunching, no negative consumption, strict clearing, choice probabilities in [0,1]. Pathologies: DIAG-1/2/3 (§2); 0.57% of mass on V<−1e6 nodes explains all 168 apparent monotonicity violations; renter-cap bunching persists. Relaxed candidate: same pathologies (1.74% invalid mass, 43.2% p_birth=0).

## 9. Calibration-search integrity

Design (verified from the local drivers; scratch record-level forensics blocked): wave 1 = 4 arms × 6 jitters (0..0.15), pattern/NM, 850 evals/330 min, seeds 2026071100+ID, seeded from the stage3 chain (17.59-era); wave 2 = 3 arms × 8 tasks, 60-draw LHS over each arm's cube (units[0] = wave-1 cross-arm strict best) then pattern polish (720 evals), seeds 2026071200+ID; smoke script present; reducer strict-min per arm and cross-arm. Strict gating airtight (in-run, select_seed, reducer); non-strict failure mode is a crash, not masking. Bounds patching real but correctly reflected in actual-bounds metadata and report CSVs; `source_controlled_bounds` mislabeled under patching; relaxed arms' diagnostic_only marker dropped after wave 1. Budget arithmetic says the LHS stage fit its deadline (unverifiable exactly without lateral_summary.json). **Effective breadth: one basin** — both winners are polish descendants of the single (probably diagnostic-arm, box-clipped) wave-1 seed; ~23,830 strict records overstate exploration around the optimum. 15/15 moments present in every scored record; no evidence any non-strict record won.

## 10. Identification

Nominal: 15 moments, 14 parameters. Effective: tenure_choice_kappa at its externally preferred corner (should be fixed, not searched); H0 gradient removed by the mis-set bound; c_bar_0/χ/h̄₀ at corners → interior parameter count ≈9-10; theta_n effectively unidentified (T3); weights ad hoc with influence concentrated (old_age_own_rate = 36-46% of both losses). Fresh one-sided FD Jacobian at the candidate (tight equilibrium, 2% box steps): `diagnostics/jacobian_current_bound.json` — see JACOBIAN ADDENDUM below.

## 11. Paper / figure staleness inventory

Full inventory: `paper_audit/RESULTS_STALENESS.md`; required-changes list: `paper_audit/PAPER_SPEC_AUDIT.md`. Bottom line: the model section correctly describes the corrected spec (incl. Rouwenhorst values, normalized bequest display with max{W,0}, finance block, supply closure, residual owner services, menu, timing, units crosswalk, 15/14 system); ALL quantitative content is the single frozen 16.13 vintage (no mixing, no 6.31-era numbers, figures byte-identical to the frozen copies); required changes: refresh to whichever candidate survives the fixes (fit prose flips sign in places), disclose the 4 at-bound estimates + bounds + Nb + verification-grid results, disclose income-process provenance (INC-1), disclose closure-variant sensitivity of policy numbers or withhold them, disclose the fertility-artifact caveat on tfr/childlessness, fix the rooms-target description (MMS 42-metro matched sample, not national ACS).

## 12. Blocked items (need Torch login refresh — `kinit`, then `ssh torch`)

1. Byte diff of `/scratch/td2248/projects/Fertility_Spring26_20260710_combined_spec` vs the local tree (mitigated: bit-exact repro).
2. `cases.jsonl` record-level forensics (duplicate evals, per-task seed traces, 23,830-record census).
3. `lateral_summary.json` per wave-2 task (`global_completed` — LHS truncation check).
4. `cross_arm_summary.json` — including the **unextracted preference-relaxed arm best** (does c_bar_0 leave 1.28 when allowed to 1.80? psi_child? kappa_fert?).
5. Wave-1 per-task `overnight_launch.json` (source_seed_path — confirms the diagnostic-arm-seed inference in S1).
6. Slurm accounting for 13313018/13314033/13314034 (zero-failure claim).

## 13. Recommended actions (NOT implemented — awaiting authorization)

Ordered by value per unit risk: (1) fix the infeasible-node fertility default (zero the artifact births or reallocate entrant mass to nearest feasible nodes) — this is model-critical numerics; spec first, then re-estimate; (2) widen H0 bound (~[1,20]) and re-polish; (3) drop tenure_choice_kappa from the search vector (fix 0); (4) replace the grid-discrete median target with an interpolated quantile or a smooth moment; (5) adopt SE-based weights from `code/data/moment_standard_errors/`; (6) re-report at Nb=240 with tight equilibrium as the headline number (search at 120 remains fine); (7) commit the production reality: the working tree, the drivers, the profile with income+finance+bequest spec inside `production_profile_overrides()` and a fingerprint check in `validate_production_profile`; (8) make every re-solve tool read the spec from the record (store overrides in best.json). Items (1)-(6) require re-estimation; the current candidates then become seeds, not results.

## 14. Artifact index

- `LEAD_WORKSHEET.md` — lead line-verification log (this audit's primary evidence trail)
- `ECONOMIC_ASSESSMENT.md` — the economics of the misses, bounds, and the relaxed arm
- `repro/*.json` + `scripts/reproduce_candidate.py`, `scripts/run_battery*.sh` — the reproduction battery
- `scripts/audit_jacobian.py`, `diagnostics/jacobian_current_bound.json` — identification
- `spec_audit/{INCOME,BEQUEST,FINANCE_SUPPLY,TIMING}_AUDIT.md`, `scripts/timing_identity_test.py`, `scripts/check_income_process.py`
- `target_audit/{TARGET_MODEL_AUDIT,TARGET_DATA_AUDIT}.md`, `scripts/recompute_loss.py`
- `search_audit/SEARCH_INTEGRITY.md`
- `paper_audit/{PAPER_SPEC_AUDIT,RESULTS_STALENESS}.md`
- `provenance/{PROVENANCE_MANIFEST.md,hashes.json}`
- `diagnostics/current_bound/` + `diagnostics/DIAGNOSTICS_REPORT.md`
- `workflow_findings.json`, `workflow_return.json` — full findings + 46 adversarial verdicts
- `FABLE_CROSSCHECK_PROMPT.md` — independent-adjudication package
