# Calibration Status

Updated: `2026-07-06` (grid-convergence finding; mechanics-first plan; Phase 1 mechanics memo DONE)

## July 6 afternoon: Phase 1 mechanics mission COMPLETE (memo delivered)

All five probes (M1-M5) done at the resolved grid, full runs on torch
(snapshot Fertility_Spring26_20260706_mechanics; bit-identical smoke vs
local). Deliverable: output/model/fable_size_mapping_audit_20260701/
MECHANICS_MEMO_20260706.md (one page per probe + Phase-2 decision sheet).
Headlines that change interpretation: (1) the ownership jump 0.46->0.75 is
buying AT the fertility deadline (age 42 = last fertile period; buy region
expands there; buyers arrive ~10 wealth cells above their dp) — not slow
saving crossing the constraint, and not bequests (calibrated theta0=4.7e-4);
(2) the b~0 young-renter mass is impatience, NOT a trap (dp -5% => zero
saving response below the threshold; relief harvests the near-threshold band
and the b~0 share GROWS 0.467->0.502); (3) at Nb=120/240 the tax2 dTFR is
+0.004..+0.009 at ALL THREE calibrations — the canonical +0.021 was mostly
grid noise, so the July-4 "weight fragility" claim dissolves; capitalization
-12.2% invariant to grid, weights, and ladders; (4) cap-pinned deriv share
0.734 at 120/240 (grid STRENGTHENED it) and ladder-INVARIANT (0.72-0.73
across five owner ladders) => the mechanism is renter cap-to-floor geometry;
D3/D5 (estimate or externally pin hR_max) is the highest-value Phase-2
decision; (5) the "dead rungs 8/10" premise dissolves — they hold 29% of ALL
owners (parent stock) and pushing them away costs -0.126 TFR; (6) M5
childlessness gradient anti-data at 120, structural (D12 unchanged).
Canonical Nb=120 full target-fit table in the memo (loss 7.4207; worst
misses: old-age own 0.890 vs 0.764, own_2534 0.196 vs 0.341, family gap
0.349 vs 0.168; housing_increment_0to1 HITS). ROUTING NOTE (Tommaso,
standing): real runs go to torch even when minutes-scale locally; local is
for smokes/analysis only.

## July 6 CRITICAL: the canonical calibration is not grid-converged

The canonical Nb=60 record (loss 6.0014895774) scores 7.43 at Nb=120 and 7.17
at Nb=240 at the SAME theta — first-order grid convergence; the 6.00 fit was
partly tuned to grid noise. Distribution/tenure-flow moments drift 3-18% from
60->120 (ownership-flow moments dominate the loss drift; the median moment is
grid-discrete but low-weight). Root causes isolated: under-resolved renter mass
at the borrowing constraint and buy thresholds; a stable +5-12% overstatement
of the effective down payment (threshold node-snap + the audit's sentinel
smear); grid redesign at Nb=60 CANNOT fix it (owners live at b<0, old wealth at
2-6, renters at 0-2 — no 60-node placement covers all; verified by a failed
dense-band A/B). At Nb=120 grid error ~= data bootstrap SEs; Nb=240 converged.
Protocol going forward: search at Nb=120, verify at Nb=240. ALL July-5/6
overnight results (transfer recalib 12.5, plateau retarget 10.7, tax-dTFR
collapse) are Nb=60-based and therefore diagnostic-direction-only pending
re-measurement at 120. At the resolved grid the REAL economic misses worsen
(own_family_gap 0.34 vs 0.17; old-age own 0.90 vs 0.76) while the h0 housing
increment basically HITS (0.653 vs 0.664 — the Nb=60 "miss" was a grid
artifact). Live plan: output/model/fable_size_mapping_audit_20260701/
PLAN_mechanics_first_20260706.md (mechanics at 120 -> model decisions -> ONE
recalibration at 120 -> verify at 240 -> policy tables). Evidence: COPY
output/model/{nb_refine_check,grid_redesign_ab,nb240_ref}.json,
grid_pathology_diag.log, canonical_diag_packet_20260706/.

## July 4 experimentation day (loop session; decision page delivered)

Full record: output/model/fable_size_mapping_audit_20260701/
DECISION_PAGE_20260704.md (FINAL). Headlines: the two re-weighted
calibrations CONVERGED (W-mech 23.6121, W-equal 20.03 under their own
unit-free objectives; DE sterile); under honest weights the paper-core
block is fittable — own_family_gap HIT (0.176/0.149 vs 0.168, from
0.296), own_rate hit, own_rate_2534 half-closed, young wealth 0.28-0.29
— at the price of old-age ownership (0.92-0.93 vs 0.764; missing exit
margin, ledger D11) and size segmentation. Mechanism stability: cap-
pinned concentration survives (0.60-0.63 vs 0.68), capitalization
identical; the tax->TFR level is the one weight-fragile headline
(+0.021 -> +0.005) and its disciplining moment (housing_increment_0to1)
is measurement-broken in both directions (h0 excludes purchases, h1
uncontrolled — ledger D13). Childlessness is structural with an
anti-data income gradient (D12). CORRECTION: the July-2 claim that the
fertility margin is the lowest income type was a column misread — the
margin is MIDDLE-income renters (z=1.0 carries 52% of derivative mass).
Ledger D1-D13 prepared; D7/D10-D13 await decisions.

## July 2 late evening: code audit, econometric review, frontier runs (ultracode session)

Full session record: `output/model/fable_size_mapping_audit_20260701/`
(CODE_AUDIT_REPORT.md, CALIBRATION_ECONOMETRIC_REVIEW.md, updated ledger
D7/D8, HANDOFF). Headlines only: (1) a 35-agent adversarially-verified
code audit found 18 correctness findings — most importantly a MEASUREMENT
tier on active targets (tenure/room stats measure the pre-choice state,
~4y-young effective windows; the "65-75" windows drop 65 and add 75-77 on
the w=160 moment; [the parity-top-coding finding was WITHDRAWN 2026-07-03:
the household-doubling convention (tfr = 2x parity; parity 2 = fertility 4)
is in the code and covers it — see the lead correction in
CODE_AUDIT_REPORT.md]; own_family_gap measures ALL parents, not new parents;
housing_increment_0to1 excludes the purchase margin; three target values
lack in-repo provenance) — none applied, all ledger-gated; (2) the
econometric review shows the ad hoc weights give the young-wealth moment
the LOWEST effective influence of all 14 (0.1%; four moments carry 83%),
5-6 of 13 parameters sit at silent bounds, and the loss has a ~5e-4
equilibrium-acceptance noise floor; (3) four cluster runs (48 tasks,
jobs 12307288/89 on the June snapshot, 12307458/60 on the NEW
`Fertility_Spring26_20260703_wealthfork_frontier` snapshot) are in
flight: R1 seeded+unseeded DE insurance under the current objective, R2
young-wealth mean→weighted-median frontier (0.09996729, influence-
preserving weight 38.6), R3 relaxed-bounds frontier (beta≥0.90, chi≤1.40,
c_bar_0≤2.0, psi_child≤0.70, h_bar_0≥0.50) — ALWAYS set
INTERGEN_N_HOUSE=5 (submit default 6 is a different menu; caught by
smoke); (4) a safe efficiency pack is committed in the COPY (bitwise-
identical regression at Nb=60; big speedups specced, not applied);
(5) D7 fork evidence complete: the July tax-alone zero-crossing was
partly a wealth-grid artifact (Nb=240: ~3% mass/one node just under the
family dp); transfers act as near-substitutes for the tax at the
affordability margin; candidate re-target moments built (liquid
0.203/0.043; SCF-2022 era fact 1.14/0.385); lambda (birth-cohort
retention) evidence: 0.60-0.75, which shrinks the phase9b scale
multiplier 1.88→~1.49.

## July 2 evening: entry/scale margin + young-wealth diagnosis (no model change; D7/D8 parked)

The wealth-moment mission (MISSION_wealth_moment_20260703.md; deliverables
in `output/model/fable_size_mapping_audit_20260701/`, mirrored in the copy)
completed all tracks. Track 0 established that the population-closure
machinery is UNREACHABLE in the canonical 5-state markov-income solve path
(every published intergen run was fixed-population; phase9 smoke:
bit-identical solves under both closure strings), and that the entry value
is ill-defined under the PSID entrant wealth distribution (~26% of entrant
weight sits at V=-1e10 infeasibility penalties; the stationary 25-35
renter distribution holds ZERO negative-wealth mass — the model erases the
data's negative-net-worth tail). The corrected counterfactual protocol runs
driver-level with no model change (`COPY/audit_drivers/
phase9b_entry_margin_scale_loop.py`, node-level entry logits, W^E
calibrated to qbar*=0.5): tax2 conclusions survive nearly unchanged
(dp -12.21 -> -11.96%, total births +1.25%); grant capitalization doubles
(+1.34 -> +3.14%, births +9.1%); the OLD-OWNER LOCK-IN WEDGE rewrites: TFR
damage -0.058 -> -0.004 per capita but city scale -5.3% (births -5.5%).
Scale responses are a LOCAL-BIRTH FEEDBACK (B0/E0), robust to the
undisciplined (qbar*, kappa_E) but hinging on local_birth_entry_weight=1.0
— see ledger D8. Track 1-2 (WEALTH_MOMENT_DIAGNOSIS.md): the 0.179 target
is PSID-SHELF NETWORTH2R/INCFAMR (non-housing net worth, weighted mean of
ratios, N=7163, weighted median 0.0999) — builder verified line-by-line;
the model's 0.350 overshoot is generated by saving behavior between entry
and 25-35 (injection is data-consistent and unclipped; beta is bound-pinned
and near-inert for young wealth — 5-point grid probe; the ywlock repair
breaks the room-gap/fertility/old blocks, +40 loss). 47% of 25-35 renters
hold less than the starter down payment; 26% could afford the family-rung
dp but rent (user-cost margin). Track 3: PARENTAL_TRANSFER_DESIGN_MEMO.md
prepares a type-flagged dp-share transfer (pi_T~0.22, tau_p~0.5, both
external, verified vs Engelhardt-Mayer/NAR/SHED/Charles-Hurst) — ledger D7,
NOT implemented, NOT decided.

## July 2 marginal-mass + rental-cap program (no model change; decisions parked)

The July-2 mission (isolated copy; deliverables in
`output/model/fable_size_mapping_audit_20260701/`, mirrored in the copy's
`mission_docs_mirror/`) completed both deliverables. MARGINAL_MASS_MAP.md:
the "68% cap-pinned" fertility-derivative concentration reproduces exactly
(0.6806), is robust to `kappa_fert` (not a logit artifact), and is governed
by the cap-to-floor distance (0.836 at `hR_max=5.5` to 0.219 at 8.0), the
jump-floor geometry (`linear_only` collapses it to 0.27), and the
young-wealth overshoot (0.39 at the data-consistent young-wealth theta);
the susceptible mass sits at the top of the renter wealth distribution,
ages 34-42, carried by the MIDDLE income types (z=1.0: 52% of derivative
mass; corrected 2026-07-04 — the July-2 'lowest income type' line misread
the location column); exact cohort-weight reconciliation puts
56-62% of the PE cap-relief TFR response on childless renters and 37-40%
on childless OWNERS (switch-to-rental option — exposure is tenure-wide).
SOFT_CAP_DESIGN_MEMO.md: within-cell ACS regression (filters replicated
verbatim from the family-size-supply builder; 1.28M renter heads) finds a
ZERO marginal rent premium above family size (m_hat 0.09, insignificant;
no bunching) — family-sized rental scarcity is availability, not price;
verified literature (DUE cap = 1.32 x owner-grid min external; KMV 2020
estimates the cap inside the calibration; Halket-Pignatti screening
microfoundation). Menu recalibration chains CLOSED: menu_dense converged
at 6.763 (above the canonical 6.0015; same floor geometry), menu_linonly
10.143 with `h_bar_n`/`chi` still bound-pinned and the parity-gap ordering
still inverted (data 1.80 vs 0.83) — the linear spec is structurally
rejected by the existing target pair. All decisions (D1-D6, esp. D3 cap
formalization and D5 estimating `hR_max` with the 0.1377 renter rooms>=6
moment) are PREPARED, NOT TAKEN, in the mission folder's
DECISIONS_LEDGER.md.

## July 1 size-mapping audit (fixed theta, no change to active model code)

A full support/measurement/cost-unit audit of the one-market model ran in the
isolated copy `~/Desktop/Projects/Fertility/Fertility_Spring26_fable_size_mapping_audit_20260701`
(baseline reproduced exactly: loss `6.0014895774`, `p=0.6759498581`). Final
report: `docs/model/fable_size_mapping_diagnosis_20260701.md`; full evidence in
`output/model/fable_size_mapping_audit_20260701/`. Verdict (decision rule B):
property tax always delivers `dp<0` (−12.2%/−21.9%, mechanical via the
rental-rate supply curve), `dYoungOwn>0`, and `dTFR>0` (~+0.02, via cheaper
birth-state ownership options at all fertile ages), but NEVER
`dPr(H>=6 | fertility-marginal young)>0` under any owner support, measurement
map, or cost-unit deformation — young buyers always take the cheapest rung,
which is below the parity-1 floor (4.52 at the current best theta). The
July-1 lambda-deformation "revival" was (i) the anchor-2 normalization, i.e.
an unfunded large-home subsidy, and (ii) biased by a Bellman-only hook
implementation: `owner_price_multiplier_by_rung` and the old-owner sale wedge
never entered the forward distribution or wealth statistics (TFR bias up to
0.19); the prior sale-wedge failure is therefore unreliable evidence, and a
consistency patch lives in the audit copy (commit `6c89ea2`). The renter cap
is the strongest structural fertility lever (hR_max 6→8: baseline TFR
+0.071) but is load-bearing for the family-ownership fit (own_family_gap
0.296→0.035). Recommended next step: a funded young(-renter)-targeted rebate
object (threshold accounting: it bridges the H=6 down-payment gap for ~30% of
the fertility-marginal band), implemented symmetrically in Bellman+KFE; not
segmentation, not sequential fertility, no blind recalibration.

## June 2026 One-Market Intergenerational Strand

A separate one-market intergenerational housing-fertility strand is active under
`code/model/intergen_housing_fertility/`. It is the current no-location
quantitative model for the June 2026 intergenerational strand: simplified
relative to the older spatial center-periphery model by dropping location. It
is **not** yet a final production calibration.

June 30 morning checkpoint. The beta-and-chi-capped strict fixed-income chain
has improved again. As of the morning scan, wave 3 (`12043716`) is still
running, so this is the best-so-far across the completed strict run, completed
overnight waves 1--2, and the current wave-3 checkpoints. Current best:

- Loss: `6.0014895774`.
- Market residual: `8.60e-05`.
- Equilibrium price: `p = 0.6759498581`.
- Result JSON:
  `/scratch/td2248/projects/Fertility_Spring26_20260629_fixed_income_beta_chicap_3h/code/cluster/results_intergen_housing_fertility_intergen_fixed_income_beta_chicap_localpolish_overnight_w3_20260629/task_10/best.json`.
- Algorithm label/case: `local_polish_pattern`, `pattern_i000_d07_p`, case
  `23`.

Full current target table:

| Moment | Target | Model | Gap | Weight | Contribution |
|---|---:|---:|---:|---:|---:|
| `tfr` | 1.918000 | 1.969671 | 0.051671 | 20.000000 | 0.053398 |
| `childless_rate` | 0.188000 | 0.257704 | 0.069704 | 20.000000 | 0.097173 |
| `own_rate` | 0.575472 | 0.524712 | -0.050761 | 100.000000 | 0.257666 |
| `own_family_gap` | 0.167662 | 0.296228 | 0.128566 | 45.000000 | 0.743817 |
| `housing_increment_0to1` | 0.664435 | 0.572957 | -0.091477 | 14.000000 | 0.117153 |
| `old_parent_childless_nonhousing_wealth_to_income_gap_6575` | 1.007450 | 0.746537 | -0.260913 | 2.000000 | 0.136151 |
| `prime30_55_childless_renter_mean_rooms` | 3.805288 | 3.668796 | -0.136492 | 6.000000 | 0.111781 |
| `prime30_55_childless_owner_share_rooms_ge6` | 0.596131 | 0.607019 | 0.010888 | 25.000000 | 0.002964 |
| `old_nonhousing_wealth_to_income_median_6575` | 2.230461 | 2.002945 | -0.227516 | 0.800000 | 0.041411 |
| `young_childless_renter_liquid_wealth_to_annual_gross_income_2535` | 0.179226 | 0.350023 | 0.170798 | 12.000000 | 0.350061 |
| `prime30_55_childless_owner_minus_renter_mean_rooms` | 2.418762 | 2.165314 | -0.253447 | 12.000000 | 0.770826 |
| `old_age_own_rate` | 0.764261 | 0.874370 | 0.110109 | 160.000000 | 1.939852 |
| `own_rate_2534` | 0.341166 | 0.214583 | -0.126583 | 80.000000 | 1.281860 |
| `prime30_55_parent_3plus_minus_1to2_mean_rooms` | 0.367700 | 0.257373 | -0.110327 | 8.000000 | 0.097376 |

Full current free-parameter table. The JSON stores `theta.beta` as the
four-year-period discount factor; the table reports the annualized value used
for the search bound, with the period value in the note.

| Parameter | Lower | Estimate | Upper | Flag |
|---|---:|---:|---:|---|
| `beta_annual` | 0.940000 | 0.940020 | 0.995000 | near lower; `theta.beta=0.780816` |
| `alpha_cons` | 0.400000 | 0.619019 | 0.950000 |  |
| `c_bar_0` | 0.080000 | 1.279999 | 1.280000 | near upper |
| `c_bar_n` | 0.050000 | 0.093973 | 1.500000 |  |
| `h_bar_0` | 1.000000 | 1.000998 | 6.000000 | near lower |
| `h_bar_jump` | 0.050000 | 2.151760 | 2.500000 |  |
| `h_bar_n` | 0.020000 | 1.367641 | 2.000000 |  |
| `psi_child` | 0.000000 | 0.350000 | 0.350000 | near upper |
| `kappa_fert` | 1.000000 | 1.026725 | 12.000000 |  |
| `tenure_choice_kappa` | 0.000000 | 0.000001 | 0.120000 | near lower |
| `chi` | 0.400000 | 1.148915 | 1.150000 | near upper |
| `theta0` | 0.000000 | 0.000467 | 2.000000 | near lower |
| `theta_n` | 0.000000 | 0.897232 | 1.500000 |  |

After this readout, wave 3 (`12043716`) was cancelled at user request to save
compute. Three short diagnostic experiments were launched from the same
best-so-far seed:

- Coordinate-probe / local pattern diagnostic: Slurm `12061934`, run tag
  `intergen_beta_chicap_coordprobe_20260630`, one task, `30` eval cap,
  initial/min normalized step `0.01`. Purpose: check local directions around
  the best point rather than keep polishing blindly.
- Seeded global-DE population probe: Slurm `12061935`, run tag
  `intergen_beta_chicap_globalde_probe_20260630`, `12` tasks, `140` eval cap
  per task, population `18`, normal
  `candidate_replacement_post_audit_v1` objective. Purpose: test whether a
  different population algorithm can escape the stale local basin.
- Hard young-wealth-lock global-DE probe: Slurm `12061936`, run tag
  `intergen_beta_chicap_youngwealthlock_globalde_20260630`, `12` tasks,
  `140` eval cap per task, population `18`. This uses the same target values
  and strict beta/chi bounds but overrides
  `young_childless_renter_liquid_wealth_to_annual_gross_income_2535` weight
  from `12` to `12000` at runtime. Purpose: ask whether the model can match
  young childless renter liquid wealth at almost any cost to the rest of the
  objective.

These three diagnostics completed cleanly (`0:0`). The coordinate probe
(`12061934`) and normal seeded global-DE population probe (`12061935`) did not
beat the current seed; their best is still the same loss `6.0014895774` with
young childless renter wealth `0.350023` versus target `0.179226`. The hard
young-wealth-lock probe (`12061936`) found a best stress-objective candidate at
`/scratch/td2248/projects/Fertility_Spring26_20260629_fixed_income_beta_chicap_3h/code/cluster/results_intergen_housing_fertility_intergen_beta_chicap_youngwealthlock_globalde_20260630/task_5/best.json`.
It reduces young childless renter wealth to `0.219795`, much closer to the
target, but at normal-objective loss `46.1513` and stress-objective loss
`65.8817`. The tradeoff is severe: `tfr=2.5167`, childless rate `0.0741`,
childless renter rooms `4.7890`, childless owner-renter room gap `1.0401`,
old ownership `0.9441`, and young ownership `0.1925`. Interpretation: under
the current strict beta/chi bounds and product/entry setup, matching young
wealth is feasible only by breaking the rest of the calibration badly.

Policy diagnostic on the hard young-wealth-lock candidate. To test whether a
candidate with more mass near the young-wealth target generates stronger
fertility responses, fixed-theta policy packets were run from the stress best.
Scratch outputs:

- Standard LTV-95 and other policy packet:
  `/scratch/td2248/projects/Fertility_Spring26_20260629_fixed_income_beta_chicap_3h/output/model/policy_youngwealthlock_20260630/`.
- Extreme credit packet:
  `/scratch/td2248/projects/Fertility_Spring26_20260629_fixed_income_beta_chicap_3h/output/model/policy_youngwealthlock_extreme_credit_20260630/`.
- Local summary CSV copies:
  `output/model/policy_youngwealthlock_20260630/policy_summary_ltv95_and_other.csv`
  and
  `output/model/policy_youngwealthlock_20260630/policy_summary_extreme_credit.csv`.

Baseline stress candidate has `tfr=2.516689`, childlessness `0.074067`,
ownership `0.481139`, young ownership `0.192549`, and young childless renter
wealth `0.219795`. LTV-95 credit relief raises ownership strongly but does not
raise fertility: new-parent/all-parent LTV-95 gives `tfr=2.512969`
(`-0.003720`) and ownership `0.551347` (`+0.070209`); universal LTV-95 gives
`tfr=2.506412` (`-0.010277`) and ownership `0.667498` (`+0.186359`). Full
parent down-payment removal (`parent_ltv100`) gives only `tfr=2.524488`
(`+0.007799`) while ownership rises to `0.595606` (`+0.114468`) and young
ownership to `0.348645` (`+0.156095`). Universal full down-payment removal
lowers fertility to `2.505615` (`-0.011074`) while ownership rises to
`0.772863` (`+0.291724`). Interpretation: even in the stress candidate,
credit policies mainly move tenure/ownership; the fertility response is tiny
and not robustly positive.

June 29 beta-and-chi-capped overnight chain launch. The active strict
fixed-income calibration array `12041214` is still running, but has already
improved to a current best loss `6.1271523` at
`/scratch/td2248/projects/Fertility_Spring26_20260629_fixed_income_beta_chicap_3h/code/cluster/results_intergen_housing_fertility_intergen_fixed_income_beta_chicap_localpolish_3h_20260629/task_10/best.json`.
That current best has residual `4.40e-05`, `beta_annual=0.94`,
`chi=1.1490`, and `tenure_choice_kappa=0.0`; both substantive constraints are
binding. The full current target table before the overnight chain begins is:

| Moment | Target | Model | Gap | Weight | Contribution |
|---|---:|---:|---:|---:|---:|
| `tfr` | 1.918000 | 1.981515 | 0.063515 | 20.000000 | 0.080684 |
| `childless_rate` | 0.188000 | 0.254918 | 0.066918 | 20.000000 | 0.089561 |
| `own_rate` | 0.575472 | 0.528112 | -0.047360 | 100.000000 | 0.224297 |
| `own_family_gap` | 0.167662 | 0.305278 | 0.137616 | 45.000000 | 0.852222 |
| `housing_increment_0to1` | 0.664435 | 0.580835 | -0.083600 | 14.000000 | 0.097845 |
| `old_parent_childless_nonhousing_wealth_to_income_gap_6575` | 1.007450 | 0.661651 | -0.345799 | 2.000000 | 0.239154 |
| `prime30_55_childless_renter_mean_rooms` | 3.805288 | 3.659281 | -0.146007 | 6.000000 | 0.127908 |
| `prime30_55_childless_owner_share_rooms_ge6` | 0.596131 | 0.604409 | 0.008278 | 25.000000 | 0.001713 |
| `old_nonhousing_wealth_to_income_median_6575` | 2.230461 | 2.002945 | -0.227516 | 0.800000 | 0.041411 |
| `young_childless_renter_liquid_wealth_to_annual_gross_income_2535` | 0.179226 | 0.352275 | 0.173050 | 12.000000 | 0.359354 |
| `prime30_55_childless_owner_minus_renter_mean_rooms` | 2.418762 | 2.175818 | -0.242944 | 12.000000 | 0.708262 |
| `old_age_own_rate` | 0.764261 | 0.874196 | 0.109935 | 160.000000 | 1.933703 |
| `own_rate_2534` | 0.341166 | 0.214499 | -0.126667 | 80.000000 | 1.283572 |
| `prime30_55_parent_3plus_minus_1to2_mean_rooms` | 0.367700 | 0.263137 | -0.104562 | 8.000000 | 0.087466 |

The current parameter table is:

| Parameter | Lower | Estimate | Upper | Flag |
|---|---:|---:|---:|---|
| `beta_annual` | 0.940000 | 0.940000 | 0.995000 | near lower |
| `alpha_cons` | 0.400000 | 0.619316 | 0.950000 |  |
| `c_bar_0` | 0.080000 | 1.280000 | 1.280000 | near upper |
| `c_bar_n` | 0.050000 | 0.090733 | 1.500000 |  |
| `h_bar_0` | 1.000000 | 1.000000 | 6.000000 | near lower |
| `h_bar_jump` | 0.050000 | 2.151170 | 2.500000 |  |
| `h_bar_n` | 0.020000 | 1.364754 | 2.000000 |  |
| `psi_child` | 0.000000 | 0.350000 | 0.350000 | near upper |
| `kappa_fert` | 1.000000 | 1.029219 | 12.000000 | near lower |
| `tenure_choice_kappa` | 0.000000 | 0.000000 | 0.120000 | near lower |
| `chi` | 0.400000 | 1.149004 | 1.150000 | near upper |
| `theta0` | 0.000000 | 0.000000 | 2.000000 | near lower |
| `theta_n` | 0.000000 | 0.894375 | 1.500000 |  |

An overnight continuation chain was queued to start only after `12041214`
finishes. A small setup job scans all completed `task_*/best.json` records and
writes the true best seed before each wave. Jobs:
setup `12043711` -> wave 1 `12043712` -> setup `12043713` -> wave 2
`12043714` -> setup `12043715` -> wave 3 `12043716`. Each wave is a
24-task local-polish array with `03:55:00` Slurm wall time, internal
`INTERGEN_MINUTES=225`, `INTERGEN_LOCAL_MAX_EVALS=2000`, same strict bounds
`beta_annual in [0.94,0.995]` and `chi in [0.4,1.15]`, target set
`candidate_replacement_post_audit_v1`, `J=17`, `Nb=60`, five Markov income
states, `H_own=[2,4,6,8,10]`, `hR_max=6.0`, and `max_iter_eq=10`. Result roots
are under
`/scratch/td2248/projects/Fertility_Spring26_20260629_fixed_income_beta_chicap_3h/code/cluster/`
with run tags
`intergen_fixed_income_beta_chicap_localpolish_overnight_w{1,2,3}_20260629`.
Monitor with
`squeue -j 12041214,12043711,12043712,12043713,12043714,12043715,12043716`.

June 29 fixed-income / beta-and-chi-capped local-polish launch. At user
request, the fixed-income `chi<=1.15` run `12037478` was cancelled after about
`2h22m`; partial outputs were preserved. Its best partial candidate before
cancellation was task `10`, loss `8.5950410`, residual `5.95e-05`,
`beta_annual=0.9236`, `chi=1.15`, `own_rate_2534=0.1820`, aggregate
ownership `0.4978`, old ownership `0.8669`, and young childless renter wealth
`0.4124`. This confirmed that the optimizer was working but still wanted both
a high owner-service premium and an annual discount factor below an acceptable
lifecycle range.

A stricter diagnostic run was then staged from the same fixed-income code,
with search bounds changed to `beta_annual in [0.94,0.995]` and
`chi in [0.4,1.15]`. The seed is the partial best from the cancelled run,
clipped to `beta_annual=0.94` and `chi=1.15`:
`/scratch/td2248/projects/Fertility_Spring26_20260629_fixed_income_beta_chicap_3h/output/model/best_seed_beta_chicap.json`.
Slurm smoke job `12041178` completed cleanly (`0:0`) with empty stderr and
strict residual `1.90e-06`; the clipped seed scores loss `14.8929053`.

The full stricter optimizer run is Slurm job `12041214`, array `1-24%24`,
wall time `03:10:00`, internal budget `INTERGEN_MINUTES=175`, target set
`candidate_replacement_post_audit_v1`, `J=17`, `Nb=60`, five Markov income
states, `H_own=[2,4,6,8,10]`, `hR_max=6.0`, `max_iter_eq=10`. Launch command
from the scratch `code/cluster` directory:
`env INTERGEN_PYTHON=/share/apps/anaconda3/2025.06/bin/python INTERGEN_RUN_TAG=intergen_fixed_income_beta_chicap_localpolish_3h_20260629 INTERGEN_TARGET_SET=candidate_replacement_post_audit_v1 INTERGEN_J=17 INTERGEN_NB=60 INTERGEN_INCOME_STATES=5 INTERGEN_N_HOUSE=5 INTERGEN_MAX_ITER_EQ=10 INTERGEN_MINUTES=175 INTERGEN_LOCAL_MAX_EVALS=1000 INTERGEN_SEED_BASE=2026062970 INTERGEN_SEED_THETA_JSON=/scratch/td2248/projects/Fertility_Spring26_20260629_fixed_income_beta_chicap_3h/output/model/best_seed_beta_chicap.json INTERGEN_LOCAL_MIN_STEP=0.001 sbatch --parsable --array=1-24%24 --time=03:10:00 --mem=4G submit_intergen_housing_fertility_local_polish.sh`.
Early health check: all 24 tasks running; all tasks wrote `best.json` and
`cases.jsonl`; after 63 total cases the best was task `11`, loss `12.7981`,
residual `6.15e-06`, with `beta_annual=0.94` and `chi=1.15` still binding.
Monitor with `squeue -j 12041214`; result root:
`/scratch/td2248/projects/Fertility_Spring26_20260629_fixed_income_beta_chicap_3h/code/cluster/results_intergen_housing_fertility_intergen_fixed_income_beta_chicap_localpolish_3h_20260629/`.

June 29 fixed-income / chi-capped local-polish cluster test. Claude's separate
workspace model copy at
`/Users/tommasodesanto/Desktop/intergen_housing_fix_claude_work/code/model/`
identified and fixed the age-18 income-profile bug, enabled working-age mean
income-profile normalization, and tightened the default wealth grid. For the
cluster test, a staged copy was created under
`tmp/intergen_fixed_income_chicap_model/` and synced to Torch scratch
`/scratch/td2248/projects/Fertility_Spring26_20260629_fixed_income_chicap_3h/`.
The only additional staged search restriction is `chi <= 1.15` in
`intergen_housing_fertility/local_panel.py`; this is a calibration/search
restriction, not a deep model change. Slurm smoke job `12037475` completed
cleanly (`0:0`) with empty stderr and reproduced the constrained seed at loss
`20.0712674`, strict market residual `2.28e-05`, `chi=1.15`,
`own_rate_2534=0.0947`, and young childless renter wealth `0.655`.

The full optimizer run is local-polish, not global-DE-only random search:
Slurm job `12037478`, array `1-24%24`, wall time `03:10:00`, internal budget
`INTERGEN_MINUTES=175`, `INTERGEN_LOCAL_MAX_EVALS=1000`, target set
`candidate_replacement_post_audit_v1`, `J=17`, `Nb=60`, five Markov income
states, `H_own=[2,4,6,8,10]`, `hR_max=6.0`, `max_iter_eq=10`, and seed
`/scratch/td2248/projects/Fertility_Spring26_20260629_fixed_income_chicap_3h/output/model/best_seed_chicap.json`.
The exact launch was from the scratch `code/cluster` directory with:
`env INTERGEN_PYTHON=/share/apps/anaconda3/2025.06/bin/python INTERGEN_RUN_TAG=intergen_fixed_income_chicap_localpolish_3h_20260629 INTERGEN_TARGET_SET=candidate_replacement_post_audit_v1 INTERGEN_J=17 INTERGEN_NB=60 INTERGEN_INCOME_STATES=5 INTERGEN_N_HOUSE=5 INTERGEN_MAX_ITER_EQ=10 INTERGEN_MINUTES=175 INTERGEN_LOCAL_MAX_EVALS=1000 INTERGEN_SEED_BASE=2026062950 INTERGEN_SEED_THETA_JSON=/scratch/td2248/projects/Fertility_Spring26_20260629_fixed_income_chicap_3h/output/model/best_seed_chicap.json INTERGEN_LOCAL_MIN_STEP=0.001 sbatch --parsable --array=1-24%24 --time=03:10:00 --mem=4G submit_intergen_housing_fertility_local_polish.sh`.
Early health check: all 24 tasks were running and writing `best.json` and
`cases.jsonl`; task 1 improved from `20.07` to `18.47` within the first minute,
and an early cross-task scan found best loss `17.4002` in task 4. Monitor with
`squeue -j 12037478`; collect from
`/scratch/td2248/projects/Fertility_Spring26_20260629_fixed_income_chicap_3h/code/cluster/results_intergen_housing_fertility_intergen_fixed_income_chicap_localpolish_3h_20260629/`.

June 29 wide-bound 12-hour search completion. The 24-worker, 3-wave
wide-bound run on Torch completed cleanly: global-DE jobs
`11987205 -> 11987207 -> 11987209` and local-polish jobs
`11987206 -> 11987208 -> 11987210` all exited `0:0`. The final best remained
the wave-2 local-polish record, not wave 3:
`/scratch/td2248/projects/Fertility_Spring26_20260628_widebounds_12h/code/cluster/results_intergen_housing_fertility_intergen_age18_widebounds_local_polish_12h_w2_20260628/task_5/best.json`.
It has algorithm label `local_polish_nelder-mead / nm_shrink_524_08`, target
set `candidate_replacement_post_audit_v1`, `J=17`, `Nb=60`, five Markov income
states, `H_own=[2,4,6,8,10]`, `hR_max=6.0`, linear interpolation, strict
market convergence, rank loss `7.070863941410668`, market residual
`1.0967853758464806e-05`, and price `p=0.6459035095984136`. The local final
packet is `output/model/intergen_widebounds_final_20260629/packet/`; the pulled
best is `output/model/intergen_widebounds_final_20260629/best.json`; the exact
complete readout is
`output/model/intergen_widebounds_final_20260629/packet/complete_readout.md`.
The final target fit is:

| Moment | Target | Model | Gap | Weight | Contribution |
|---|---:|---:|---:|---:|---:|
| `prime30_55_childless_owner_minus_renter_mean_rooms` | 2.418762 | 1.975764 | -0.442998 | 12.000000 | 2.354961 |
| `young_childless_renter_liquid_wealth_to_annual_gross_income_2535` | 0.179226 | 0.617265 | 0.438040 | 12.000000 | 2.302544 |
| `old_age_own_rate` | 0.764261 | 0.860090 | 0.095829 | 160.000000 | 1.469321 |
| `own_family_gap` | 0.167662 | 0.243257 | 0.075595 | 45.000000 | 0.257159 |
| `old_nonhousing_wealth_to_income_median_6575` | 2.230461 | 1.666396 | -0.564065 | 0.800000 | 0.254535 |
| `old_parent_childless_nonhousing_wealth_to_income_gap_6575` | 1.007450 | 0.721238 | -0.286211 | 2.000000 | 0.163834 |
| `prime30_55_parent_3plus_minus_1to2_mean_rooms` | 0.367700 | 0.226530 | -0.141169 | 8.000000 | 0.159431 |
| `tfr` | 1.918000 | 1.881069 | -0.036931 | 20.000000 | 0.027278 |
| `own_rate_2534` | 0.341166 | 0.359450 | 0.018284 | 80.000000 | 0.026745 |
| `childless_rate` | 0.188000 | 0.224567 | 0.036567 | 20.000000 | 0.026743 |
| `housing_increment_0to1` | 0.664435 | 0.628067 | -0.036368 | 14.000000 | 0.018517 |
| `prime30_55_childless_renter_mean_rooms` | 3.805288 | 3.831950 | 0.026662 | 6.000000 | 0.004265 |
| `prime30_55_childless_owner_share_rooms_ge6` | 0.596131 | 0.608677 | 0.012546 | 25.000000 | 0.003935 |
| `own_rate` | 0.575472 | 0.571478 | -0.003995 | 100.000000 | 0.001596 |

The estimated free parameters are:

| Parameter | Lower | Estimate | Upper |
|---|---:|---:|---:|
| `beta_annual` | 0.880000 | 0.907737 | 0.995000 |
| `alpha_cons` | 0.400000 | 0.701129 | 0.950000 |
| `c_bar_0` | 0.080000 | 0.991115 | 1.280000 |
| `c_bar_n` | 0.050000 | 0.275608 | 1.500000 |
| `h_bar_0` | 1.000000 | 1.556840 | 6.000000 |
| `h_bar_jump` | 0.050000 | 2.157678 | 2.500000 |
| `h_bar_n` | 0.020000 | 1.133341 | 2.000000 |
| `psi_child` | 0.000000 | 0.141933 | 0.350000 |
| `kappa_fert` | 1.000000 | 3.176379 | 12.000000 |
| `tenure_choice_kappa` | 0.000000 | 0.002677 | 0.120000 |
| `chi` | 0.400000 | 1.280970 | 2.400000 |
| `theta0` | 0.000000 | 1.168628 | 2.000000 |
| `theta_n` | 0.000000 | 0.710935 | 1.500000 |

Read: this is a real improvement over the previous deterministic-kappa best
(`11.7463`), but it does not solve the housing/wealth mechanism. The fit now
matches aggregate and young ownership well, but still leaves young
childless-renter liquid wealth far too high, the childless owner-renter room
gap too small, and old-age ownership too high. `tenure_choice_kappa` is no
longer exactly zero, but remains economically close to deterministic tenure
choice.

June 27 wealth-unit update. The live solver now reports explicit annual-gross
liquid-wealth moments for young households:
`young_all_liquid_wealth_to_annual_gross_income_2530`,
`young_childless_liquid_wealth_to_annual_gross_income_2535`, and
`young_childless_renter_liquid_wealth_to_annual_gross_income_2535`, each with a
weighted median and sample mass. These are diagnostics for the wealth-target
audit and are exposed through `calibration.extract_moments`. The active
`candidate_replacement_roomgap_14moment_tfr192_v1` target set now uses
`young_childless_renter_liquid_wealth_to_annual_gross_income_2535 = 0.17922556`
instead of the old period-income statistic. The current review packet shows the old
`young_liquid_wealth_to_income` statistic is a 4-year period after-tax-income
ratio (`0.353` at `output/model/intergen_current_review/quick/`), while the
annual-gross young-childless-renter object is about `1.10`, compared with the
PSID candidate mean `0.179`. Conclusion: do not interpret historical losses
that targeted the old wealth statistic as comparable to the active target set.
The companion audit is
`output/model/intergen_current_review/wealth_units_audit/README.md`, and the
issue ledger is
`docs/model/intergen_housing_block_audit_plan_20260626.md`.
The broader active-target object audit is
`docs/model/intergen_target_object_audit.md`. It flags the same denominator
problem for the old-age nonhousing wealth-to-income moments and records that
room/ownership targets are mostly clean measurement objects, while `tfr` should
be read as completed-fertility-equivalent and `housing_increment_1to2` as an
additional-child housing-demand proxy rather than a sequential second-birth
hazard.

June 27 target-system follow-up. The active TFR-1.92 roomgap target set keeps
14 hard moments but replaces the young wealth key:
`young_liquid_wealth_to_income` ->
`young_childless_renter_liquid_wealth_to_annual_gross_income_2535`, target
`0.17922556`, weight `12.0`. With external entry wealth below, the active
search now has 13 free parameters and 14 hard moments. Historical target sets
that still contain `young_liquid_wealth_to_income` should be read as old
diagnostic systems unless revised explicitly.

June 27 entry-wealth closure update. The intergen calibration base overrides now
use an externally calibrated entrant liquid-wealth distribution rather than
searching `b_entry_fixed`. The distribution is in annual gross family-income
units, estimated from the same PSID young childless renter sample used for the
wealth target: weighted quintile-bin means of `NETWORTH2R / INCFAMR` equal
`[-2.51940697, -0.07907025, 0.10228762, 0.35287169, 3.03955200]` with
approximately 20% weight each; weighted mean `0.17922556`, weighted median
`0.09996729`. Solver mode `entry_wealth_mode="income_ratio_distribution"` maps
these ratios into model liquid-wealth units using entrant annual gross income by
location/income state, then linearly scatters the mass over adjacent wealth-grid
nodes. This keeps counterfactual house-price feedback through prices and
down-payment thresholds, but does not mechanically index initial assets to house
prices. The legacy scalar `b_entry_fixed` path remains available, but the active
local/global intergen search vector drops `b_entry_fixed` and now has 13 free
parameters.

June 27 20:10 EDT post-audit 3-hour calibration sample launch. After the
Slack target-object audit closed all 14 moments, a new target set was added as
`candidate_replacement_post_audit_v1`. It keeps 14 hard moments against the
13-parameter search vector, changes the fertility block to completed-fertility
targets (`tfr=1.918`, `childless_rate=0.188`), replaces the old PSID
`housing_increment_1to2` hard target with the ACS/MMS parent room gap
`prime30_55_parent_3plus_minus_1to2_mean_rooms=0.36769955881`, keeps the young
annual-gross renter wealth target, and fixes the old-age nonhousing wealth
model statistics to annual-income mean/median ratios. Local checks passed:
target count `14`, weight count `14`, free parameters `13`, compile clean, and
a one-evaluation local `global-de-panel` completed with strict market
convergence. A fresh Torch scratch source snapshot was synced to
`/scratch/td2248/projects/Fertility_Spring26_20260627_postaudit_run` (ignore
the interrupted partial folder ending in `_postaudit`). Remote target/compile
checks passed. Slurm preflights also passed with empty stderr: global-DE job
`11927758`, panel job `11927759`, each one evaluation, both exit `0:0`, loss
`36.9806095`, market residual `1.41e-05`. The real 12-worker sample is running:
global-DE job `11927776` and local/random-panel job `11927777`, arrays
`1-6%6` each, wall time `03:10:00`, internal `INTERGEN_MINUTES=175`,
`J=16`, `Nb=60`, `income_states=5`, `INTERGEN_N_HOUSE=5`
(`H_own=[2,4,6,8,10]`), `max_iter_eq=10`, `interp_method=linear`,
`use_pti_constraint=False`, seeded from the June 26 global-DE best record
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib/code/cluster/results_intergen_housing_fertility_intergen_fixedstats_seeded_globalde_12h_mem4_w3_20260625/task_10/best.json`.
Result roots are
`/scratch/td2248/projects/Fertility_Spring26_20260627_postaudit_run/code/cluster/results_intergen_housing_fertility_intergen_postaudit_globalde_sample3h_20260627/`
and
`/scratch/td2248/projects/Fertility_Spring26_20260627_postaudit_run/code/cluster/results_intergen_housing_fertility_intergen_postaudit_panel_sample3h_20260627/`.
Monitor with `squeue -j 11927776,11927777`. Collect after completion from
the scratch copy's `code/model` with
`python tools/collect_intergen_panel_results.py --results-dir ../cluster/results_intergen_housing_fertility_intergen_postaudit_globalde_sample3h_20260627`
and the analogous panel directory.

June 27 late completion/readout. The post-audit 3-hour sample completed cleanly:
all global-DE job `11927776` and local/random-panel job `11927777` array tasks
exited `0:0`, with empty Slurm stderr in the checked tasks. The global-DE
branch evaluated `5894` successful cases and found the best record
`de_g045_i011` in task `1`, case `1030`, saved at
`/scratch/td2248/projects/Fertility_Spring26_20260627_postaudit_run/code/cluster/results_intergen_housing_fertility_intergen_postaudit_globalde_sample3h_20260627/task_1/best.json`.
The local/random-panel branch evaluated `6553` successful cases but never beat
the warm start (`36.9806`). The local review folder is
`output/model/intergen_postaudit_sample3h_review_20260627/`. Re-solving the
best record with the live `candidate_replacement_post_audit_v1` target set,
`J=16`, `Nb=60`, `income_states=5`, `H_own=[2,4,6,8,10]`, `hR_max=6.0`,
`max_iter_eq=10`, and `interp_method=linear` wrote
`output/model/intergen_postaudit_sample3h_review_20260627/best_globalde_task1_packet/`
with rank loss `16.2934`, market residual `4.15e-05`, and equilibrium price
`p=0.694225`. Complete readout files are
`complete_readout.md`, `best_globalde_task1_packet/target_fit.md`,
`parameter_bounds.md`, `all_moments.md`, and `shortlist_top10.md` in the local
review folder. Main unresolved fit issues remain economic rather than a blank
or crashing solve: the owner-renter room gap is too small, old-age ownership is
too high, the new high-child parent room gap is too high, aggregate ownership is
too high, and young childless renter wealth-to-income is too high.

June 27 timing-convention update. The active intergen model timing has been
changed from adult entry at `age_start=22` with first fertile age `26` to adult
entry at `age_start=18` with fertile ages `18,22,26,30,34,38,42`. The standard
live horizon is now `J=17`, preserving terminal age `82`; active CLI,
local/global panel, cluster wrapper, and solver-audit defaults were updated from
`J=16` to `J=17`. A small `J=17`, `Nb=30`, one-case informed smoke solved
cleanly with strict market convergence and residual `9.20e-07`; this is only a
breakage check, not a calibration result. Target implications: completed
fertility and completed childlessness can remain the same statistical objects,
but any first-birth timing target would need to be redefined on the new age
grid. The externally calibrated entrant-wealth distribution is still sourced
from young childless renters ages 25--35; with entry now at age 18, that object
should be re-extracted or explicitly reinterpreted before a serious production
calibration. A new literal deterministic/modal policy overlay tool was added at
`code/model/tools/plot_intergen_policy_overlay.py`; for the current June 27
best cache it wrote the first essential overlay at
`output/model/intergen_postaudit_sample3h_review_20260627/best_globalde_task1_packet/policy_overlays/deterministic_policy_overlay_age30_lowhigh_renter_ownerH4.png`.
The same old best theta was also re-solved once under the new `age_start=18`,
`J=17` timing without recalibration; packet
`output/model/intergen_postaudit_sample3h_review_20260627/best_globalde_task1_packet_age18_j17_resolve/`
has loss `30.5741`, residual `9.04e-06`, price `0.705161`, `tfr=2.257`,
childlessness `0.113`, young ownership `0.581`, aggregate ownership `0.762`,
and old ownership `0.947`. This is not a valid new benchmark; it only proves
that the prior `J=16` best is stale once early adult fertility is allowed. Its
matching overlay is
`output/model/intergen_postaudit_sample3h_review_20260627/best_globalde_task1_packet_age18_j17_resolve/policy_overlays/deterministic_policy_overlay_age30_lowhigh_renter_ownerH4_age18_j17.png`.

June 28 00:40 EDT age-18/J17 overnight calibration launch. A fresh dirty-source
scratch snapshot was synced to
`/scratch/td2248/projects/Fertility_Spring26_20260628_age18_12h`. Remote
compile and target-count checks passed: `candidate_replacement_post_audit_v1`
has 14 targets, 14 weights, and the active search has 13 free parameters. The
remote lifecycle check confirmed `age_start=18`, `J=17`, terminal age `82`, and
fertile ages `18,22,26,30,34,38,42`. The start point is the stale-but-useful
June 27 post-audit global-DE record copied to
`/scratch/td2248/projects/Fertility_Spring26_20260628_age18_12h/output/model/intergen_postaudit_sample3h_review_20260627/best_globalde_task1_best.json`.
Exact Slurm preflights passed with empty stderr: global-DE job `11940278` and
local/random-panel job `11940282`, each one seeded full-grid/full-target
evaluation, both exit `0:0`, loss about `30.574`, market residual
`2.52e-05`. The full overnight search uses the current hand/documented weights
rather than optimal SMM weights; optimal weighting is deferred until the target
covariance/standard-error file is built.

The active 12-hour chain is three dependent `03:55:00` waves on `cpu_short`,
4G per task, with 24 concurrent tasks in each active wave: 12 global-DE array
tasks plus 12 seeded local/random-panel array tasks. Common settings:
`INTERGEN_TARGET_SET=candidate_replacement_post_audit_v1`, `INTERGEN_J=17`,
`INTERGEN_NB=60`, `INTERGEN_INCOME_STATES=5`, `INTERGEN_N_HOUSE=5`
(`H_own=[2,4,6,8,10]`), `INTERGEN_MAX_ITER_EQ=10`,
`INTERGEN_MINUTES=225`, `interp_method=linear`, and
`use_pti_constraint=False`. Submission command pattern from the scratch
`code/cluster` directory was:
`env <common settings> INTERGEN_RUN_TAG=<tag> INTERGEN_GLOBAL_EVALS_PER_TASK=3000 INTERGEN_SEED_BASE=<base> sbatch --parsable --array=1-12%12 --time=03:55:00 --mem=4G submit_intergen_housing_fertility_global_de.sh`
and analogously
`env <common settings> INTERGEN_RUN_TAG=<tag> INTERGEN_CASES_PER_TASK=3000 INTERGEN_SEED_BASE=<base> sbatch --parsable --array=1-12%12 --time=03:55:00 --mem=4G submit_intergen_housing_fertility_twohour_panel.sh`,
with `--dependency=afterok:<previous-wave>` on waves 2 and 3. Global-DE jobs:
`11940313 -> 11940314 -> 11940315`, seed bases `2026062800`,
`2026064800`, `2026066800`; local/random-panel jobs:
`11940316 -> 11940317 -> 11940318`, seed bases `2026072800`,
`2026074800`, `2026076800`. Result roots under the scratch copy are
`code/cluster/results_intergen_housing_fertility_intergen_age18_postaudit_globalde_12h_w{1,2,3}_20260628/`
and
`code/cluster/results_intergen_housing_fertility_intergen_age18_postaudit_panel_12h_w{1,2,3}_20260628/`.
Monitor with
`squeue -j 11940313,11940314,11940315,11940316,11940317,11940318`.

June 28 09:30 EDT search-pipeline correction. The remaining third wave of the
overnight chain was canceled at user request after about 1h16m of wave-3 wall
time (`11940315`, `11940318`; waves 1 and 2 had already completed cleanly).
The cancellation preserved all scratch outputs. The best available candidate at
the time of cancellation was still the global-DE wave-2 task-5 record
`de_g050_i021`, not partial wave 3: rank loss `24.7570822`, market residual
`3.80e-06`, price `0.664626`, `tfr=1.92894`, childlessness `0.20927`,
young ownership `0.47671`, aggregate ownership `0.67312`, old ownership
`0.92981`. The local/random-panel branch never beat the warm start, confirming
that it is not a real local optimizer and should not receive substantial compute
unless redesigned.

A true local-polish command was added to the active code:
`python -m intergen_housing_fertility.cli local-polish`. It optimizes in the
same normalized bounded parameter cube as global-DE, checkpoints every case to
`cases.jsonl`/`best.json`, and currently supports bounded Nelder-Mead and
coordinate pattern search. The Slurm wrapper is
`code/cluster/submit_intergen_housing_fertility_local_polish.sh`. Local compile
and a tiny low-grid CLI smoke passed. The updated code and best-DE seed were
synced to the same scratch snapshot; the seed file is
`/scratch/td2248/projects/Fertility_Spring26_20260628_age18_12h/output/model/intergen_age18_current_best_globalde_20260628.json`.
Exact Slurm preflight passed with empty stderr: local-polish job `11959255`,
one task, two full-grid/full-target evaluations, exit `0:0`, reproducing seed
loss `24.7570822`.

A 2-hour local-polish array is now running as job `11959283`, array `1-12%12`,
wall time `02:05:00`, internal `INTERGEN_MINUTES=115`, target/grid
`candidate_replacement_post_audit_v1`, `J=17`, `Nb=60`, `income_states=5`,
`INTERGEN_N_HOUSE=5`, `max_iter_eq=10`, 4G per task. Tasks 1--6 use
Nelder-Mead with normalized initial steps
`0.025,0.040,0.060,0.085,0.115,0.150`; tasks 7--12 use coordinate pattern
search with the same step grid. Monitor with `squeue -j 11959283`. Result root:
`/scratch/td2248/projects/Fertility_Spring26_20260628_age18_12h/code/cluster/results_intergen_housing_fertility_intergen_age18_postaudit_local_polish_2h_20260628/`.
Early health check after the first few evaluations shows the polish already
improving the seed from `24.7570822` to `23.2019100` in task 8
(`pattern_i000_d10_m`), with `tfr=1.92515`, childlessness `0.20987`, young
ownership `0.43177`, aggregate ownership `0.64774`, and old ownership
`0.92031`.

June 28 local-polish completion/readout. Job `11959283` completed cleanly: all
12 array tasks exited `0:0`. Final best is task 12, label
`pattern_i014_d12_p`, saved locally at
`output/model/intergen_age18_local_polish_2h_20260628/best.json`. A clean local
packet re-solve with `candidate_replacement_post_audit_v1`, `J=17`, `Nb=60`,
5 income states, `H_own=[2,4,6,8,10]`, `hR_max=6.0`, `max_iter_eq=10`, and
linear interpolation reproduces rank loss `11.7463453`, market residual
`7.32e-05`, and price `0.664872`. Full readout with every target row, every
free parameter/bound, and all extracted model moments:
`output/model/intergen_age18_local_polish_2h_20260628/packet/complete_readout.md`.
Classic diagnostic packet:
`output/model/intergen_age18_local_polish_2h_20260628/packet/`, including
`contact_sheet.png`, `first_look_policies_markets*.png`, wealth and
total-wealth densities, age profiles, tenure by age, owner thresholds/rungs,
room-bin fit, standard diagnostics, and deterministic policy overlay
`policy_overlays/deterministic_policy_overlay_age30_lowhigh_renter_ownerH4.png`.
The fit improvement is real, but the optimum sits on several search bounds:
`beta_annual` lower, `c_bar_0` upper, `h_bar_0` lower, `h_bar_jump` near upper,
`psi_child` upper, and `tenure_choice_kappa` lower (`0`, deterministic tenure).

June 28 positive-tenure-noise search launch. In response to the bound/corner
diagnosis, the active search bound for `tenure_choice_kappa` was changed from
`[0.000,0.080]` to `[0.005,0.080]` in
`code/model/intergen_housing_fertility/local_panel.py`; this constrains the
calibration search away from literal deterministic tenure choice but does not
remove the model's ability to solve explicitly specified `kappa=0` points.
Local compile and a tiny low-grid local-polish smoke passed. The updated code
was synced to the existing Torch scratch snapshot
`/scratch/td2248/projects/Fertility_Spring26_20260628_age18_12h`, along with
the prior deterministic-corner best seed at
`output/model/intergen_age18_local_polish_2h_best_20260628.json`. Remote
compile and target/bounds checks passed: 14 targets, 13 free parameters,
`tenure_choice_kappa` bound `[0.005,0.080]`. Full-grid Slurm preflights also
passed with empty stderr: global-DE job `11969650` and local-polish job
`11969670`, each exit `0:0`; clipping the prior seed to `kappa=0.005` gives
loss `13.5780780`, residual `2.39e-05`, price `0.663730`, `tfr=1.86778`,
young ownership `0.36872`, aggregate ownership `0.58986`, old ownership
`0.87450`, and young childless renter wealth/income `0.64386`.

A hard 4-hour mixed positive-kappa search is now running: global-DE job
`11969716` and true local-polish job `11969717`, each array `1-12%12`, wall
time `03:55:00`, 4G per task. Common settings:
`candidate_replacement_post_audit_v1`, `J=17`, `Nb=60`, `income_states=5`,
`INTERGEN_N_HOUSE=5`, `max_iter_eq=10`, `INTERGEN_MINUTES=225`, seeded from
the prior best. Global-DE uses population size `24` and up to `3000`
evaluations/task; local-polish uses up to `1200` evaluations/task, min step
`0.002`, and the existing step-grid split between Nelder-Mead and coordinate
pattern search. Result roots:
`code/cluster/results_intergen_housing_fertility_intergen_age18_poskappa_globalde_4h_20260628/`
and
`code/cluster/results_intergen_housing_fertility_intergen_age18_poskappa_local_polish_4h_20260628/`
under the same scratch copy. Monitor with
`squeue -j 11969716,11969717`. Early checkpoint scan after roughly 7--9 cases
per task already improved the constrained seed from `13.5780780` to
`13.5070457` in local-polish task 10 (`pattern_i000_d11_m`), with all best
records respecting `tenure_choice_kappa=0.005`.

June 28 young-wealth frontier diagnostic launch. A diagnostic target set
`candidate_replacement_post_audit_wealthstress_v1` was added with the same 14
target values as `candidate_replacement_post_audit_v1`, but with
`young_childless_renter_liquid_wealth_to_annual_gross_income_2535` weight raised
from `12.0` to `120.0`. This is not a production SMM target set; it is a
frontier run to learn which moments/parameters break when the optimizer is
forced to lower young childless renter liquid wealth. Local compile and target
metadata checks passed: 14 targets, 13 free parameters, wealth weight `120.0`,
positive `tenure_choice_kappa` bound `[0.005,0.080]`. Remote compile and exact
Slurm preflight passed with empty stderr: local-polish preflight job `11970841`,
exit `0:0`; the seed point has wealth-stress loss `36.8937570`, residual
`2.39e-05`, and young wealth `0.64386`.

The wealth-frontier search is running in parallel with the positive-kappa main
search. Seed is the current positive-kappa checkpoint best at launch time
(`loss=13.2051449` under the regular post-audit objective). Jobs:
global-DE `11970889` (`1-6%6`) and local-polish `11970890` (`1-12%12`), wall
time `02:05:00`, internal `INTERGEN_MINUTES=115`, target set
`candidate_replacement_post_audit_wealthstress_v1`, `J=17`, `Nb=60`, 5 income
states, `INTERGEN_N_HOUSE=5`, `max_iter_eq=10`, 1G per task. Because the
regular positive-kappa search is still running, this wealth-frontier run may
partially pend under `QOSMaxCpuPerUserLimit`; at launch all 6 global tasks and
4 local-polish tasks were running while local-polish tasks 5--12 were pending.
Monitor with `squeue -j 11970889,11970890`. Result roots:
`code/cluster/results_intergen_housing_fertility_intergen_age18_wealthstress_globalde_2h_20260628/`
and
`code/cluster/results_intergen_housing_fertility_intergen_age18_wealthstress_local_polish_2h_20260628/`
under the same scratch copy.

June 28 positive-kappa and wealth-frontier completion. Jobs `11969716`,
`11969717`, `11970889`, and `11970890` all completed cleanly with exit
`0:0`. The best regular positive-kappa record is local-polish task 12,
`pattern_i019_d08_p`, with loss `13.1713650`, residual `4.90e-05`, price
`0.665065`, and `tenure_choice_kappa=0.005`; it does not beat the earlier
deterministic-kappa local-polish best loss `11.7463453`. The best
wealth-frontier diagnostic record is local-polish task 9,
`pattern_i022_d11_m`, with wealth-stress loss `19.9825422`, residual
`1.25e-03`, price `0.666666`, and `tenure_choice_kappa=0.005`. It lowers
`young_childless_renter_liquid_wealth_to_annual_gross_income_2535` from about
`0.659` in the regular positive-kappa best to `0.314`, closer to the target
`0.179`, but worsens the lifecycle/room block. Full side-by-side target and
parameter tables are saved in
`output/model/intergen_policy_counterfactuals_current/wealthstress_diagnostic/wealthstress_vs_regular_completion_readout.md`.

June 28 wealth-frontier policy diagnostic. The one-market policy proof-of-
concept runner now accepts direct `best.json` records plus explicit `J`,
`n_house`, `income_states`, and `target_set` arguments. A fixed-theta policy
packet was run from the wealth-frontier best at
`output/model/intergen_policy_counterfactuals_current/wealthstress_diagnostic/`
with `J=17`, `Nb=60`, `n_house=5`, 5 income states, and `max_iter_eq=25`.
Compared with the wealth-frontier baseline, parent-targeted LTV relief to
`phi=0.95` raises the completed-fertility-equivalent `tfr` by only `0.0059`
and lowers childlessness by `0.0011`; applying the relief to all parent states
is numerically identical in this run. Universal `phi=0.95` lowers `tfr` by
`0.0445` while raising young ownership; `hR_max=5` lowers `tfr` by `0.0464`;
the property-tax and estate-tax diagnostics have near-zero fertility effects.
This is not a production policy result, but it confirms that the wealth-frontier
diagnostic improves the targeted parent-credit effect only modestly.

Follow-up fixed-price check: holding the wealth-frontier GE baseline price fixed
at `p=0.66666564`, parent-targeted `phi=0.95` raises the completed-fertility-
equivalent `tfr` by `0.00983` and lowers childlessness by `0.00135`; extending
eligibility to all parent states is again identical. A fixed-price parent-`phi`
sweep from `0.80` to `1.00` gives a maximum `tfr` increase of only `0.01153`
at zero down payment. Output:
`output/model/intergen_policy_counterfactuals_current/wealthstress_fixed_price_parent_credit/`.
Conclusion: GE price feedback is not the main reason parent credit has weak
fertility effects in this candidate; the household-level exposed fertility
margin is itself thin.

June 28 wide-bound overnight calibration launch. At user request, an additional
overnight frontier search was launched to see whether the current
`candidate_replacement_post_audit_v1` target system can improve further before
redesigning tenure segmentation. This is a search-box experiment, not a model
equation change. The active global/local search bounds in
`code/model/intergen_housing_fertility/local_panel.py` were widened from the
post-audit box to:
`beta_annual=[0.880,0.995]`, `alpha_cons=[0.400,0.950]`,
`c_bar_0=[0.080,1.280]`, `c_bar_n=[0.050,1.500]`,
`h_bar_0=[1.000,6.000]`, `h_bar_jump=[0.050,2.500]`,
`h_bar_n=[0.020,2.000]`, `psi_child=[0.000,0.350]`,
`kappa_fert=[1.000,12.000]`, `tenure_choice_kappa=[0.000,0.120]`,
`chi=[0.400,2.400]`, `theta0=[0.000,2.000]`, and
`theta_n=[0.000,1.500]`. Local compile passed. A local full-grid two-evaluation
smoke of `local-polish` from the current deterministic best reproduced loss
`11.7463453` with residual `7.32e-05`; the first coordinate probe solved
cleanly. A fresh Torch scratch snapshot was synced to
`/scratch/td2248/projects/Fertility_Spring26_20260628_widebounds_12h`.
Remote compile and target checks passed: 14 targets, 14 weights, 13 free
parameters. Exact Slurm preflights passed with empty stderr: global-DE job
`11987103` and local-polish job `11987104`, each one full-grid evaluation,
both exit `0:0`, reproducing the seed loss `11.7463453`.

The full chain is three dependent 3:55 waves on `cpu_short`, 24 concurrent
tasks per active wave: 12 global-DE plus 12 local-polish. Common settings:
`INTERGEN_TARGET_SET=candidate_replacement_post_audit_v1`, `J=17`, `Nb=60`,
5 income states, `INTERGEN_N_HOUSE=5`, `max_iter_eq=10`, and seed theta
`output/model/intergen_age18_local_polish_2h_20260628/best.json` under the
scratch copy. Global-DE uses `3000` max evals/task, population size `28`,
mutation `0.90`, crossover `0.70`, and `INTERGEN_MINUTES=225`. Local-polish
uses `1200` max evals/task, `INTERGEN_MINUTES=225`, and min step `0.0015`.
Job chain: global-DE `11987205 -> 11987207 -> 11987209`; local-polish
`11987206 -> 11987208 -> 11987210`. Result roots under the scratch copy are
`code/cluster/results_intergen_housing_fertility_intergen_age18_widebounds_globalde_12h_w{1,2,3}_20260628/`
and
`code/cluster/results_intergen_housing_fertility_intergen_age18_widebounds_local_polish_12h_w{1,2,3}_20260628/`.
Monitor with `squeue -j 11987205,11987206,11987207,11987208,11987209,11987210`.

June 24 code-audit fixes changed the live mechanics before any new calibration:
`use_pti_constraint=False` by default, optional PTI uses actual transaction
debt, and `chi` is now an owner utility/service premium on residual owner
housing services after physical family-space needs, not a multiplier that
reduces the physical room floor.

June 25 tenure-segmentation correction: the live diagnostic convention now caps
renter housing at `hR_max=6.0` while keeping the owner ladder
`H_own=[2,4,6,8,10]`. This keeps strong product-support separation while
allowing renters to reach the first large-room threshold; the upper 8/10-room
family-sized rungs remain owner-only in this diagnostic convention. Treat
pre-correction diagnostic plots with `hR_max=8.0` as useful pathology evidence,
not as the current model convention.

June 25 (late) audit + measurement fixes. A code/diagnostics audit found and
fixed three issues in the live model; all committed to `main`
(`7219f64`, `8f97ed6`, `8ade95f`). Full handoff:
`docs/model/intergen_calibration_handoff_20260625.md`.
(1) `compute_markov_statistics` collapsed the income dimension before applying
nonlinear (threshold/median) operators to the renter policy (Jensen error), so
`prime30_55_childless_renter_share_rooms_ge6` read `0.013` instead of the true
`0.124` (≈ target `0.138`); the renter median and cap shares were also corrupted.
Fixed via `markov_renter_room_moments`. (2) `housing_increment_0to1` was a raw
12-year (3-period) post−pre window with no control, reading `1.467` (mostly
lifecycle drift) vs target `0.664`; redefined as a controlled difference-in-
differences (birth cohort minus no-birth control) at a configurable
`housing_event_horizon` (default `0` = birth period ≈ 3 years), now `0.477`.
(3) added a configurable wealth grid and a real `interp_method` switch
(`linear` default / `monotone_cubic`); the default solve is unchanged. The
"childless" room moments were checked and are NOT a bug (model `current_child_bin_dt`
matches the ACS `ni==0` "childless-in-household" definition), so
`owner_share_rooms_ge6` (`0.964` vs `0.596`) is a genuine economic miss, not
measurement. After the fixes, the current diagnostic point's aggregate loss is
`29.85` (down from `38.05`; verified re-solve, residual `1.5e-5`). The dominant
remaining misses are all economic: young ownership (`own_rate_2534≈0` vs `0.341`),
aggregate ownership (`0.371` vs `0.575`), old-age ownership (`0.916` vs `0.764`),
owner large-room share, and young liquid wealth. IMPORTANT: any prior cluster
run or incumbent that weighted renter room-share/median or `housing_increment`
was optimizing against corrupted moments and must be re-scored with the fixed
stats before reuse.

June 25 21:20 EDT calibration launch on Torch. Local `main` was pulled
(`git pull --ff-only`, already up to date), then a clean committed scratch copy
was created at
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib` from commit
`47e5542` (`Seed intergen calibration from warm-start theta`). The required
fix commits are in the scratch history: `7219f64`, `8f97ed6`, `8ade95f`.
Warm-start file copied to
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib/output/model/intergen_room_distribution_current_best_20260623/summary.json`.
Remote verification passed:
`PYTHONPATH=$PWD python -m compileall -q intergen_housing_fertility tools/build_intergen_mechanics_packet.py tools/collect_intergen_panel_results.py`
and
`PYTHONPATH=$PWD python -m intergen_housing_fertility.cli smoke --quiet`.
Exact Slurm preflight also passed on `cpu_short`: global-DE job `11810090`
and local-panel job `11810091`, each one full-grid/full-target warm-start
evaluation, both exit `0:0`, stderr empty, corrected loss `29.8504933`,
market residual `4.29e-05`. A single 12-hour CPU job was rejected by Torch
(`partition 'cpu_short' is not valid for this job`; same for `cpu_prem`), so
the run was launched as three dependent `03:55:00` waves. An initial 6G-per-task
chain (`11810113`--`11810118`) was canceled after Torch admitted only 20 tasks
because of `QOSMaxMemoryPerUser`; the final chain uses 4G per task and the
first wave is running all 24 concurrent tasks (12 global-DE array tasks plus 12
seeded local-panel array tasks). Target/grid:
`candidate_replacement_young_old_roomgap_v1`,
`J=16`, `Nb=60`, `income_states=5`, `INTERGEN_N_HOUSE=5` (therefore
`H_own=[2,4,6,8,10]` under `base_overrides`), `hR_max=6.0`,
`max_iter_eq=10`, `interp_method=linear` default, `use_pti_constraint=False`.
Global-DE jobs: `11810145 -> 11810147 -> 11810149`; local-panel jobs:
`11810146 -> 11810148 -> 11810150`. First wave was running at status-write
time; waves 2 and 3 are dependency-held. Result roots are:
`code/cluster/results_intergen_housing_fertility_intergen_fixedstats_seeded_globalde_12h_mem4_w1_20260625/`,
`..._globalde_12h_mem4_w2_20260625/`, `..._globalde_12h_mem4_w3_20260625/`,
`..._panel_12h_mem4_w1_20260625/`, `..._panel_12h_mem4_w2_20260625/`, and
`..._panel_12h_mem4_w3_20260625/` under the scratch copy. Monitor with
`squeue -j 11810145,11810146,11810147,11810148,11810149,11810150` and logs
under `code/cluster/logs/slurm_ihf_de_<job>_<task>.out` and
`code/cluster/logs/slurm_ihf_2hr_<job>_<task>.out`. Collect after completion
from the scratch copy with:
`python tools/collect_intergen_panel_results.py --results-dir ../cluster/results_intergen_housing_fertility_<RUN_TAG>`.

June 26 morning completion/readout. The full three-wave Torch chain completed
cleanly: jobs `11810145`--`11810150` all exited `0:0`, Slurm stderr files were
empty, and `squeue` showed no remaining tasks. The useful branch was global-DE;
the seeded local-panel branch never beat the warm start. Best saved global-DE
record: `de_g044_i022` from wave 3 task 10,
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib/code/cluster/results_intergen_housing_fertility_intergen_fixedstats_seeded_globalde_12h_mem4_w3_20260625/task_10/best.json`,
saved rank loss `16.9807648`, market residual `5.43e-05`. A local reproducible
inspection folder was created at
`output/model/intergen_fixedstats_overnight_review_20260626/`. Re-solving the
best candidate through `tools/build_intergen_mechanics_packet.py` with the live
fixed-stats grid (`J=16`, `Nb=60`, 5 Markov income states,
`H_own=[2,4,6,8,10]`, `hR_max=6.0`, `max_iter_eq=10`, `interp_method=linear`)
wrote
`output/model/intergen_fixedstats_overnight_review_20260626/best_de_g044_i022_packet/`
with re-solved rank loss `17.071` and market residual `6.31e-06`. Shortlist
comparison artifacts are
`output/model/intergen_fixedstats_overnight_review_20260626/shortlist_comparison.md`
and `shortlist_key_moments.png`. The complete calibration report with every
target moment, target value, model moment, weight, loss contribution, and all
estimated parameters with search bounds is
`output/model/intergen_fixedstats_overnight_review_20260626/best_complete_calibration_report.md`.
Main readout: aggregate ownership and the
first-child housing response improved materially, but the old-age ownership
absorbing-margin problem remains (`old_age_own_rate≈0.978` vs target `0.764`),
young ownership is still too low (`own_rate_2534≈0.124`--`0.127` in the
re-solve/best record vs `0.341`), `housing_increment_1to2` is too low
(`≈0.16` vs `0.488`), renter large-room share is too low (`≈0.054` vs `0.138`),
and owner large-room share remains too high (`≈0.737` vs `0.596`). Visual
diagnostics in the packet show small market residuals and no obvious blank or
exploding policy object; the remaining failures look economic/discrete-choice
rather than an immediate numerical crash.

June 26 housing-grid / smoothing audit. The dedicated audit note is
`docs/model/intergen_housing_smoothing_audit_20260626.md`. Two read-only tools
were added: `code/model/tools/audit_intergen_active_tenure_values.py` and
`code/model/tools/run_intergen_kappa_smoothing_audit.py`. The active tenure
value-gap diagnostic exactly recovers stored tenure probabilities from saved
policies and next-period values (max error `2.65e-08` on the dense-grid GE
cache). Findings: mid-wealth wiggles around `b≈4.2` are genuine close tenure
and owner-rung margins (`rent` vs `H4`, and `H8` vs `H6` for parents), but the
large high-mass consumption drop at `b_entry_fixed≈0.146514` is not a simple
owner/renter near-indifference issue; owner branches are infeasible or far below
the renter value at the near-zero states. A dense-grid fixed-theta kappa sweep
over `tenure_choice_kappa ∈ {0,0.005,0.01,0.02,0.05}` shows that smoother
tenure choice does not remove the entry-node drop: the drop stays near `-0.34`
through `kappa=0.02` and is still `-0.308` at `kappa=0.05`. At `kappa=0.05`,
economics move materially (old-age ownership about `0.900`, renter large-room
share about `0.111`, owner large-room share about `0.636`), so it is not a
harmless numerical smoothing fix. Next housing audit should focus on liquid-grid
placement and the entry/liquid wealth atom before changing fertility mechanics.

June 26 calibration-plumbing update after the smoothing audit:
`tenure_choice_kappa` is now part of the local/global calibration search vector
in `code/model/intergen_housing_fertility/local_panel.py`, with bounds
`[0.000, 0.080]` and a backward-compatible warm-start default of `0.01` for
older seed theta files. A new explicit 14-hard-moment target set was added as
`candidate_replacement_roomgap_14moment_v1` in
`code/model/intergen_housing_fertility/calibration.py`. It implements the
June 23 14-moment roomgap rescore convention: keep the replacement wealth,
fertility, ownership, and room moments; keep `own_rate`, `own_rate_2534`,
`old_age_own_rate`, and `own_family_gap`; replace separate owner/renter room
levels with `prime30_55_childless_owner_minus_renter_mean_rooms`; and demote
the old diagnostic-only moments `prime30_55_childless_owner_mean_rooms`,
`prime30_55_childless_renter_share_rooms_ge6`, and
`old_age_parent_childless_gap`. The target count is now exactly 14 against the
14-parameter search vector including `tenure_choice_kappa`.

June 26 14:25 EDT one-hour Torch smoke calibration launch. Local `main` was
pushed at commit `20eacf9` (`Add intergen roomgap 14 moment target set`), and
the plain scratch snapshot at
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib` was refreshed from
`git archive HEAD`; `SYNC_COMMIT.txt` records full SHA
`20eacf9c2c4db7ce2d77402b23ea1592f7e18ccb`. Remote compile passed after the
snapshot refresh, and the remote target check returned 14 parameters / 14
moments / no weight mismatch for `candidate_replacement_roomgap_14moment_v1`.
Warm-start file:
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib/output/model/intergen_room_distribution_current_best_20260623/summary.json`
(top-level `theta`; missing `tenure_choice_kappa` defaults to `0.01`).
Two exact-loop Slurm preflights completed cleanly with empty stderr:
global-DE job `11863502` and local-panel job `11863503`, each one seeded
full-grid evaluation, rank loss `27.60915015410047`, market residual
`4.289731608961121e-05`. The actual smoke is a small 4+4 split, not a formal
production calibration: four global-DE array tasks and four local/random-panel
array tasks, each with `INTERGEN_MINUTES=55`, target set
`candidate_replacement_roomgap_14moment_v1`, `J=16`, `Nb=60`,
`income_states=5`, `INTERGEN_N_HOUSE=5`, `max_iter_eq=10`,
`interp_method=linear` default, `use_pti_constraint=False`, and
`tenure_choice_kappa` searched over `[0.000, 0.080]`. Slurm jobs:
global-DE `11865125`, local-panel `11865127`; both arrays were running and had
written per-task metadata at status-write time. Result roots:
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib/code/cluster/results_intergen_housing_fertility_intergen_roomgap14_kappa_globalde_smoke1h_20260626/`
and
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib/code/cluster/results_intergen_housing_fertility_intergen_roomgap14_kappa_panel_smoke1h_20260626/`.
Monitor with:
`squeue -j 11865125,11865127`. Collect each branch after completion from
`code/model` with:
`python tools/collect_intergen_panel_results.py --results-dir ../cluster/results_intergen_housing_fertility_intergen_roomgap14_kappa_globalde_smoke1h_20260626`
and the analogous panel result directory.

June 26 15:23 EDT solver-accuracy diagnostics launch. Two read-only tools were
added and pushed to `main`: `code/model/tools/audit_intergen_solver_accuracy.py`
and `code/model/tools/audit_intergen_savings_globality.py` (latest launch
commit `3a8d531`, full SHA
`3a8d531a4e870bd9323d739c1f8f57fc8b330e01`). The Torch scratch snapshot at
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib` was refreshed from
that commit; `SYNC_COMMIT_DIAGNOSTICS.txt` records the full SHA. Diagnostic
inputs copied to
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib/output/model/intergen_solver_accuracy_inputs_20260626/`:
`de_w3_task10_de_g044_i022_best.json` and
`best_de_g044_i022_solution_cache.pkl`. Remote compile/import preflight passed.
This is separate from the 14-moment smoke calibration above. Three `cpu_short`
jobs were launched in parallel, all running at status-write time:
grid convergence `11868769`, equilibrium/equivalence `11868770`, and
savings-globality `11868771`. Output root:
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib/output/model/intergen_solver_accuracy_20260626/`.
The grid job writes fixed-price and full-GE `Nb ∈ {60,120,240}` target/moment
drift tables and solution caches. The equilibrium/equivalence job writes
`max_iter_eq`/`scalar_market_refine` sensitivity, interior-renter Euler
residuals, shape/KFE checks, Python-vs-Numba small-grid array differences,
fast-stats-vs-full differences, and fixed-price linear-vs-monotone-cubic
differences. The savings job writes all-branch dense-grid `b'` globality checks
from the saved cache, separating raw branch failures from KFE-relevant failures
using `mass × tenure_probability`. Monitor with
`squeue -j 11868769,11868770,11868771`; logs are
`code/cluster/logs/slurm_ihf_acc_grid_<job>.out`,
`slurm_ihf_acc_eq_<job>.out`, and `slurm_ihf_acc_save_<job>.out`.

June 26 solver-accuracy diagnostics completion/readout. Jobs `11868769`,
`11868770`, and `11868771` all completed `0:0` with empty stderr, and outputs
were pulled locally to
`output/model/intergen_solver_accuracy_20260626/`. Main grid finding: increasing
`Nb` from 60 to 240 at the fixed source theta materially changes target moments.
Fixed-price loss falls `16.4486 -> 14.7749`; GE loss falls
`16.4486 -> 14.8566`. Fixed-price and GE drifts are very similar, so the
dominant measured error here is household/grid approximation rather than
price-feedback error. Key GE moment drifts from `Nb=60` to `Nb=240`: young
ownership `0.127 -> 0.192`, aggregate ownership `0.589 -> 0.617`, renter mean
rooms `3.888 -> 3.811`, young liquid wealth/income `0.405 -> 0.330`, TFR
`1.867 -> 1.895`, and `housing_increment_1to2` `0.158 -> 0.080`. The old-age
ownership pathology barely moves (`0.978 -> 0.977`), so that remains economic /
structural rather than a coarse-grid artifact. Equilibrium-quality check:
with scalar refinement on, `max_iter_eq ∈ {3,5,10,25}` gives small residuals
(`3.4e-06` to `5.4e-05`). With scalar refinement off, low iteration budgets
fail (`0.1265` residual at 3, `0.0331` at 5, both triggering the +100 penalty);
by 10 iterations the damped loop reaches `7.0e-05`. This confirms the audit:
the live one-market run is fine with Brent on, but low-iteration/no-refinement
runs are not production-safe. Savings-globality check: 1,287 feasible branch
rows, 297 raw dense-grid gaps above `5e-3`, but only 23 KFE-relevant failures
after weighting by `mass × tenure_probability`; weighted p95 gap is `0.00309`
overall, `0.00040` for renter branches, `0.0212` for `own_H4`, `0.00406` for
`own_H6`, and near zero for the upper owner rungs. Huge unweighted gaps mostly
come from zero-weight / penalty-cliff branches and should be read as a mask/
reporting hazard, not as mass-relevant optimizer failure. Shape/KFE checks:
total mass is conserved, grid-edge mass is zero, income/child transition row
sums are clean, but an approximate same-tenure continuation check finds small
positive mass near `-1e10` penalty interpolation (`~2e-4` at `Nb=60`, falling to
`~6e-7` at `Nb=240`). Python-vs-Numba small-grid moments are close
(e.g. ownership diff `-1.5e-4`, TFR diff `-3.3e-5`), while fixed-price
monotone-cubic interpolation moves some moments nontrivially
(`housing_increment_1to2 -0.041`, owner-renter room gap `-0.066` in the small
test), so PCHIP remains a sensitivity experiment rather than a default. Interior
renter Euler residuals were computed only for classified smooth states; they
are not yet a clean accuracy certificate and should be treated as a flag for
better branch-specific residual construction.

June 26 18:53 EDT TFR-equivalent 1.92 target smoke launch. A new target-set
variant was added and pushed on `main` at commit `a9c7e59`
(`Add intergen tfr192 target variant`):
`candidate_replacement_roomgap_14moment_tfr192_v1`. It is identical to
`candidate_replacement_roomgap_14moment_v1` except that the fertility target is
`tfr=1.92` instead of `1.70`; the old target set is unchanged. This is a
diagnostic target experiment motivated by the distinction between period TFR and
stationary completed fertility in the one-shot model. The target count remains
14 against the 14-parameter search vector, with no target/weight mismatch. Local
one-case smoke passed on the live grid. The Torch scratch copy at
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib` was refreshed from
full SHA `a9c7e5953fc2b242f3576028a44d05d7d3fd64dd`; the scratch marker is
`SYNC_COMMIT_TFR192.txt`. Seed theta is the current fixed-stats overnight best:
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib/output/model/intergen_tfr192_inputs_20260626/seed_de_g044_i022_best.json`
(copied from
`output/model/intergen_fixedstats_overnight_review_20260626/candidates/de_w3_task10_de_g044_i022_best.json`).
Remote target/import/compile checks passed under the cluster Anaconda Python.
Two exact-loop Slurm preflights also passed with empty stderr: global-DE job
`11881021` and local-panel job `11881022`, each one full-grid seeded evaluation
with rank loss `15.9486492`, market residual `5.43e-05`,
`TFR-equiv=1.86681`, and childlessness `0.23218`. The actual run is a small 4+4
split, not a production calibration: global-DE job `11881031` and local-panel
job `11881032`, each with `INTERGEN_MINUTES=55`, target set
`candidate_replacement_roomgap_14moment_tfr192_v1`, `J=16`, `Nb=60`,
`income_states=5`, `INTERGEN_N_HOUSE=5`, `max_iter_eq=10`,
`interp_method=linear` default, `use_pti_constraint=False`, and
`tenure_choice_kappa` searched over `[0.000,0.080]`. Result roots:
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib/code/cluster/results_intergen_housing_fertility_intergen_tfr192_globalde_smoke1h_20260626/`
and
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib/code/cluster/results_intergen_housing_fertility_intergen_tfr192_panel_smoke1h_20260626/`.
Monitor with `squeue -j 11881031,11881032`. After completion, collect from
`code/model` with
`python tools/collect_intergen_panel_results.py --results-dir ../cluster/results_intergen_housing_fertility_intergen_tfr192_globalde_smoke1h_20260626`
and the analogous panel directory. Next step is to re-solve the best TFR-1.92
candidate through `tools/build_intergen_mechanics_packet.py --run-policy-cases`
and compare fertility/completed-fertility deltas, not losses, against the old
baseline policy packet.

June 26 19:09 EDT TFR-equivalent 1.92 smoke completion/readout. Jobs
`11881031` and `11881032` completed cleanly (`0:0` for all eight array tasks);
all Slurm stderr files are empty. Local pulled review folder:
`output/model/intergen_tfr192_smoke_review_20260626/`. The run produced `960`
finite/ok records. Scalar-best under the new `tfr=1.92` target is still the old
fixed-stats incumbent: rank loss `15.94865`, `TFR-equiv=1.86681`, model
`E[n]=0.93340`, childlessness `0.23218`, ownership `0.58900`, young ownership
`0.12709`, old ownership `0.97830`, `H01=0.76281`, `H12=0.15820`. The search did
find high-fertility frontier points, but they are worse on the overall target
bundle. Lowest-loss candidate with `TFR-equiv >= 1.90`:
`output/model/intergen_tfr192_smoke_review_20260626/candidates/best_tfr_ge190.json`,
loss `31.13658`, `TFR-equiv=1.90673`, `E[n]=0.95337`, childlessness `0.20388`,
ownership `0.47594`, young ownership `0.11086`, old ownership `0.80918`,
`H01=0.85530`, `H12=0.24209`. Closest-to-1.92 candidate under loss `<80`:
`.../candidates/closest_tfr192_loss_lt80.json`, loss `47.04853`,
`TFR-equiv=1.92188`, childlessness `0.21810`, ownership `0.47236`, young
ownership `0.04792`, old ownership `0.92371`, `H01=0.21938`, `H12=0.33394`.
All target moments for selected candidates are in
`output/model/intergen_tfr192_smoke_review_20260626/selected_candidate_comparison.md`;
all selected parameters with bounds are in
`output/model/intergen_tfr192_smoke_review_20260626/selected_parameter_bounds.md`.
Standard policy cases were re-solved for `best_tfr_ge190` and
`closest_tfr192`; see
`output/model/intergen_tfr192_smoke_review_20260626/policy_readout.md`.
Main policy result: retargeting fertility upward alone does not restore a large
parent-credit channel. Parent LTV95 raises `TFR-equiv` by only `+0.00181` at
`best_tfr_ge190` and `+0.00240` at `closest_tfr192`, versus `+0.00118` at the
old baseline. The previously large `+0.047` partial-equilibrium effect remains
a rental-cap/space-support result (`hR_max=5`), not something generated by just
moving the fertility target to `1.92` under the current `hR_max=6` support.

June 26 21:34 EDT full overnight TFR-equivalent 1.92 launch. Goal: run a full
overnight search on `candidate_replacement_roomgap_14moment_tfr192_v1` after the
small smoke showed a high-fertility frontier but no scalar-best replacement.
The Torch scratch copy at
`/scratch/td2248/projects/Fertility_Spring26_20260625_calib` was refreshed from
local commit `d23d3badfd50e040776d2bd6965c3cf7d88cc32b`; marker file
`SYNC_COMMIT_TFR192_OVERNIGHT.txt`. Two seed files are used:
old fixed-stats incumbent
`output/model/intergen_tfr192_inputs_20260626/seed_de_g044_i022_best.json`
(`TFR-equiv=1.86681`, loss `15.94865` under the 1.92 target) and high-fertility
frontier seed
`output/model/intergen_tfr192_inputs_20260626/seed_best_tfr_ge190.json`
(`TFR-equiv=1.90673`, loss `31.13658`). Remote target/import/compile checks
passed under cluster Anaconda Python. Exact-loop Slurm preflights passed with
empty stderr: old-seed DE `11883059`, high-seed DE `11883060`, old-seed panel
`11883061`, high-seed panel `11883062`; each re-solved its seed with the
expected loss/residual. Overnight design: 24 concurrent tasks per wave, split
as 6 old-seed global-DE, 6 high-seed global-DE, 6 old-seed local/random-panel,
and 6 high-seed local/random-panel. Each task has `INTERGEN_MINUTES=225`,
`--time=03:55:00`, `--mem=4G`, `J=16`, `Nb=60`, `income_states=5`,
`INTERGEN_N_HOUSE=5`, `max_iter_eq=10`, `interp_method=linear` default,
`use_pti_constraint=False`, and `tenure_choice_kappa` searched over
`[0.000,0.080]`. Three dependent waves were submitted per branch:
old-seed DE `11883112 -> 11883113 -> 11883114`; high-seed DE
`11883115 -> 11883116 -> 11883117`; old-seed panel
`11883118 -> 11883119 -> 11883120`; high-seed panel
`11883121 -> 11883122 -> 11883123`. First wave started with all 24 tasks
running; later waves are dependency-held with `afterok` dependencies. First-wave
metadata/checkpoint files were confirmed under all four result roots. Monitor:
`squeue -j 11883112,11883113,11883114,11883115,11883116,11883117,11883118,11883119,11883120,11883121,11883122,11883123`.
Result roots under the scratch `code/cluster/` directory are
`results_intergen_housing_fertility_intergen_tfr192_overnight_oldseed_globalde_w{1,2,3}_20260626/`,
`..._highseed_globalde_w{1,2,3}_20260626/`,
`..._oldseed_panel_w{1,2,3}_20260626/`, and
`..._highseed_panel_w{1,2,3}_20260626/`. Morning collection should use
`tools/collect_intergen_panel_results.py` on each result root, then re-solve the
best scalar candidate and at least the best high-fertility frontier candidate
through `tools/build_intergen_mechanics_packet.py --run-policy-cases`. Policy
readout should report completed-fertility/TFR-equivalent and childlessness
deltas first, not calibration losses.

Current reference diagnostic point: global-DE diagnostic best from
`output/model/intergen_globalde_final_best_diagnostics/source_record.json`,
label `de_g008_i011`, stored loss `11.503191936648555` under the default
one-market owner ladder. The June 10 verified fixed-stack replication with
kernel clamps, Brent scalar refinement, and `max_iter_eq=3` gives loss
`11.601065514992962`, market residual `5.219e-05`, and zero material borrowing
floor violations. This point still has serious lifecycle ownership pathologies
and should not be presented as a calibrated benchmark.

June 17 exploratory cluster pulses use the same one-market/no-location
intergenerational strand with Markov income shocks, `J=16`, `Nb=60`,
`income_states=5`, `n_house=5`, and `max_iter_eq=3`. The default
`candidate_no_timing_v0` pulse fixed the owner-median room rung but selected a
low-ownership basin: best loss `13.631`, ownership `0.125`, family ownership
gap `0.023`, old-age ownership `0.623`, TFR `1.556`. A follow-up
`candidate_no_timing_ownheavy_v1` pulse raises ownership-related weights as a
basin-finding diagnostic, not a final SMM objective. It completed `960` cases
with `954` finite records and zero Slurm stderr; best recorded rank loss
`11.755`, ownership `0.441`, family ownership gap `0.122`, old-age ownership
`0.954`, TFR `1.590`, childless rate `0.272`, and housing user-cost share
`0.420`. Interpretation: the low-ownership failure is partly objective-ranking,
but the current basin still overstates old-age ownership and housing costs and
does not yet deliver a production quantitative calibration.

The same day, a two-hour `candidate_no_timing_refinement_v1` pulse restored
stronger penalties on housing user-cost share, old-age ownership, fertility,
wealth, and housing increments while keeping ownership pressure high. Slurm job
`11026855` completed `7,901` cases, with `7,882` finite records and zero stderr
bytes. The best recorded rank loss was `20.790`, with TFR `1.642`, childless
rate `0.285`, ownership `0.391`, family ownership gap `0.170`, old-age
ownership `0.920`, old-age parent-childless gap `0.071`, housing user-cost
share `0.368`, and owner median rooms `6.0`. Frontier slices point to a target
or measurement mismatch in the one-market model: no finite candidate in the run
jointly had ownership above `0.50`, old-age ownership below `0.85`, and housing
user-cost share below `0.32`; candidates with owner median rooms equal to `6`
typically carry too-high housing costs and old-age ownership. Do not launch
longer blind searches on this exact target bundle before auditing moment
definitions and deciding which moments are hard targets for this simplified
strand.

Overnight on June 17, an intergen frontier ensemble was launched to separate
target-definition failures from weighting failures. The ensemble uses four
diagnostic target sets added to `calibration.py`: core feasibility
(`candidate_no_timing_core_feasibility_v1`), cost test
(`candidate_no_timing_cost_test_v1`), old-age test
(`candidate_no_timing_oldage_test_v1`), and room-cost test
(`candidate_no_timing_roomcost_test_v1`). Each target set runs four chained
two-hour `cpu_short` waves, with `J=16`, `Nb=60`, `income_states=5`,
`n_house=5`, `max_iter_eq=3`, `pop_size=32`, and `INTERGEN_MINUTES=115`.
Slurm chains are: core `11038977 -> 11038978 -> 11038986 -> 11038987`; cost
`11038988 -> 11038989 -> 11038990 -> 11038991`; old-age
`11038992 -> 11038993 -> 11038995 -> 11038996`; room-cost
`11038997 -> 11038998 -> 11038999 -> 11039000`. Result tags are
`intergen_overnight_<label>_w<wave>_20260617`. The morning readout should
compare feasible frontier slices across target sets, not just the scalar best.

Morning June 18 readout: all overnight chains completed with zero stderr bytes
and exit code `0:0` in `sacct`. Across the four target sets, the core
feasibility problem is easy relative to the full hard target bundle: best core
loss `2.033`, with TFR `1.888`, childless `0.173`, ownership `0.487`, family
ownership gap `0.224`, housing increments `0.843/0.599`, renter median rooms
`4.04`, and owner median rooms `6.0`. Adding the cost target gives best loss
`7.264` but still leaves housing user-cost share at `0.359`. Adding old-age
targets gives best loss `4.800`, ownership `0.503`, old-age ownership `0.894`,
old-age parent-childless gap `0.113`, but housing user-cost share rises to
`0.504`. Dropping room medians and keeping cost plus old-age targets gives best
loss `5.895`, ownership `0.550`, housing user-cost share `0.272`, but owner
median rooms collapse to `4`. The frontier diagnostic is sharp: no target set
found a candidate with ownership above `0.50`, housing user-cost share below
`0.32`, and owner median rooms equal to `6`; and only the core-feasibility run
found even one candidate with ownership above `0.50`, old-age ownership below
`0.85`, and housing user-cost share below `0.32`, but that candidate had
TFR `2.851`, negative family ownership gap, negative second-child housing
increment, and owner median rooms `4`. Interpretation: the one-market model can
match the basic fertility/tenure/room block when old-age and cost-share moments
are demoted, but the measured cost-share, old-age ownership, and owner-room
targets are not jointly coherent as hard targets without further measurement
or mechanism audits. This
demotion was diagnostic only: any formal SMM target revision must preserve
identification by replacing unreachable moments with moments that discipline
the same parameter blocks, or by fixing the affected parameters externally.

Sensitivity/Jacobian audit, June 18: `code/model/tools/audit_intergen_sensitivity_jacobian.py`
was run on Torch as Slurm job `11077696` using the best
`core_feasibility_v1` and `roomcost_test_v1` frontier points. All 54
finite-difference solves completed successfully; outputs were pulled to
`output/model/intergen_sensitivity_jacobian_20260618/` and summarized in
`docs/model/intergen_sensitivity_jacobian_audit_20260618.md`. The core point
has full local rank 13 but is badly conditioned, with condition number
`2.69e4`; the room-cost point has rank 12 because the owner-median-room target
is locally flat. Interpretation: the target count is formally adequate, but
owner median rooms, old-age ownership, and aggregate housing user-cost share are
weak or non-smooth identifying objects. Replace them with identification-
preserving moments or fix the affected parameters externally before calling a
new target system a production SMM calibration.

Replacement-moment data map, June 18:
`docs/model/intergen_replacement_moment_data_map_20260618.md` records the
recommended empirical direction: ACS/MMS household heads are the primary source
for housing, tenure, rooms, and renter cost moments; AHS is a robustness/stock
source; PSID is the source for old-age wealth, completed-fertility, lifecycle
wealth, and down-payment/access moments. All proposed replacement targets still
need a data re-audit before entering the SMM objective. In particular, source
files, sample restrictions, weights, household unit, room/bedroom convention,
wealth definition, timing, and formula must be checked against the model object.
For old-age bequest moments, the recommended primary object is nonhousing or
liquid wealth; total net worth is a robustness check.

Candidate replacement values have now been extracted and summarized in
`docs/model/intergen_candidate_replacement_targets_20260618.md`. The key ACS
candidate values are childless owner mean rooms `6.224`, childless renter mean
rooms `3.805`, owner-renter room gap `2.419`, owner share with at least six
rooms `0.596`, renter share with at least six rooms `0.138`, and childless
renter median rent-to-income `0.240`. The key PSID old-age candidate values are
old nonhousing net worth to income `6.419` by mean and `2.230` by median, a
parent-minus-childless nonhousing net-worth-to-income mean gap of `1.007`, and
an old-age parent-minus-childless ownership gap of `0.083`. These values are
not final SMM targets until the data-target audit is complete.

A 72-case local panel using `candidate_replacement_v1` was run on June 18; see
`docs/model/intergen_candidate_replacement_v1_panel_20260618.md`. The panel
completed all submitted cases with `J=16`, `Nb=60`, `income_states=5`,
`n_house=5`, and `max_iter_eq=3`. Best rank loss improved from baseline
`48.883` to `34.009`, but the economic fit is not acceptable as a benchmark:
the best case has TFR `1.952` versus target `1.700`, childlessness `0.230`
versus `0.150`, ownership `0.750` versus `0.575`, old nonhousing
wealth-to-income `2.214` versus `6.419`, parent-childless old nonhousing wealth
gap `0.268` versus `1.007`, childless renter mean rooms `4.621` versus
`3.805`, and childless owner mean rooms `5.273` versus `6.224`. High-old-wealth
candidates exist but collapse ownership, while high-owner-room candidates also
put renters into too much space. Treat `candidate_replacement_v1` as a
diagnostic target set until the old-wealth data object and owner-renter room
separation mechanism are re-audited.

Follow-up cluster wave, launched June 18 at 15:25 EDT: see
`docs/model/intergen_replacement_cluster_wave_20260618.md`. Three 8-task
global-DE arrays are running on Torch: `11095203` for `candidate_replacement_v1`
(old nonhousing mean), `11095204` for `candidate_replacement_nh_median_v1`
(old nonhousing median matched to a model median statistic), and `11095205` for
`candidate_replacement_total_median_v1` (old total-wealth median plus total
parent-childless wealth gap). A separate Jacobian job, `11095206`, audits the
current best replacement point and the high-old-wealth/low-ownership basin under
`candidate_replacement_v1`. These runs preserve 13 target moments for 13 varied
internal parameters; they are a diagnostic comparison, not a production
calibration.

Readout from the same cluster wave: all jobs completed with exit code `0:0` and
the collected results are under
`output/model/cluster_pulls/intergen_replacement_cluster_wave_20260618/`.
The old-nonhousing-mean target set improved from local-panel loss `34.009` to
`26.553`, but still misses old nonhousing mean wealth badly (`2.590` versus
`6.419`) and has a small owner-renter room gap (`0.581`). The old-nonhousing
median variant has much better scalar loss, `9.093`, with ownership `0.509`,
TFR `1.882`, childlessness `0.248`, childless renter rooms `4.052`, childless
owner rooms `5.556`, and old nonhousing median wealth `1.601` versus target
`2.230`. The total-wealth median variant has loss `13.430`, matches old total
median wealth (`5.364` versus `5.264`) and produces a larger room gap, but
pushes ownership down to `0.393` and first-birth room growth up to `1.334`.
No run produced a convincing joint fit of moderate ownership, old wealth, and a
large owner-renter room gap. The Jacobian audit is full rank but ill-conditioned
at the high-ownership basin (rank `13`, condition `9.41e3`) and rank deficient
at the high-old-wealth basin (rank `12`, condition `4.46e8`). Old nonhousing
wealth loads primarily on `beta`, not `theta0`; the current old-wealth targets
therefore do not cleanly identify the intended bequest block.

For the intergen strand, use:

- package: `code/model/intergen_housing_fertility/`
- venv: `code/model/.venv/bin/python`
- proof-of-concept runner: `code/model/tools/run_intergen_policy_poc.py`
- parent-credit margin audit: `code/model/tools/audit_intergen_parent_credit_margin.py`
- pathology audit: `code/model/tools/audit_intergen_final_best_pathologies.py`

The older center-periphery `dt_cp_model` calibration status below remains
historical/live for that strand, but it is not sufficient orientation for the
June 2026 one-market intergenerational work.

The remaining sections below preserve the older discrete-time center-periphery
status history. Historical MATLAB status notes, handoffs, and planning memos
have been archived; use them only as background unless a future session
explicitly returns to that spatial strand.

## Active Codebase

The active model implementation is the Python port:

- model setup and targets:
  `code/model/dt_cp_model/parameters.py`
- equilibrium solver:
  `code/model/dt_cp_model/solver.py`
- SMM objective:
  `code/model/dt_cp_model/objective.py`
- direct-geometry calibration logic:
  `code/model/dt_cp_model/direct_calibration.py`
- local/cluster calibration entry point:
  `code/model/tools/calibrate_direct_geometry.py`
- result collector:
  `code/model/tools/collect_direct_geometry_results.py`
- active Torch launcher:
  `code/cluster/submit_python_direct_geometry_overnight.sh`

Archived MATLAB reference code:

- `calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/`
- `calibration_archive/legacy_matlab_2026-05-07/code_matlab/`

## Current Closure

Benchmark closure:

- `population_closure = outside_option_benchmark_normalized`
- The model is stationary but not fixed-population in counterfactuals.
- At the benchmark, the stationary scale identity is imposed mechanically with
  \(S=1\):
  \[
  S E_0(p)=q^E(p)\left[M+S B_0(p)\right],
  \qquad
  M=\frac{E_0(p)}{q^E(p)}-B_0(p).
  \]
- The empirical outside-margin normalization is an outside-origin entrant
  share, not \(q^E\) itself. Baseline uses
  \(s^{E,\mathrm{out}}=0.169\), based on ACS 2012--2023 cumulative
  across-CBSA arrival rates from ages 18--22 in the MMS metro sample.
- The code's required city-entry probability is then
  \(q^{E,*}=(1-s^{E,\mathrm{out}})/(B_0/E_0)\). Using the previous
  near-calibrated \(B_0/E_0\simeq 0.93\) gives \(q^{E,*}\simeq 0.89\).
  Because the benchmark equilibrium is invariant to \(q^{E,*}\), the exact
  \(q^{E,*}\) should be recomputed after the new best fit to hit
  \(s^{E,\mathrm{out}}\) exactly.
- `outside_entry_flow` and `outside_value` are benchmark accounting objects,
  not calibrated SMM parameters under the default direct-geometry setup.
- `implied_total_population` is not an SMM target under this closure.
- Counterfactuals should hold benchmark \(M\), \(\bar W^E\), and \(\kappa_E\)
  fixed, then let \(S=q^E(p)M/[E_0(p)-q^E(p)B_0(p)]\) move.

Finite-scale guardrail:

- Candidates with \(E_0(p)-\rho B_0(p)\le 0\), negative outside entry flow, or
  nonfinite scale diagnostics should be rejected before objective comparison.

## Current Target System

Ownership targets now use ACS/MMS household heads rather than person records,
so adult children living in parent-owned homes are not counted as owners. The
audit script is `code/data/mms_center_periphery/audit_ownership_targets.R`; its
diagnostic packet is written to
`code/data/mms_center_periphery/output_ownership_audit/`.

Current ownership targets in `parameters.py`:

| Moment | Target | Source |
|---|---:|---|
| `own_rate` | `0.575` | ACS/MMS heads, ages 30--55, DUE housing restrictions |
| `own_gradient` | `0.186` | ACS/MMS heads, ages 30--55, periphery minus center |
| `own_family_gap` | `0.168` | ACS/MMS heads, ages 30--55, new parents minus childless |
| `old_age_own_rate` | `0.764` | ACS/MMS heads, ages 65--75, DUE housing restrictions |
| `old_age_parent_childless_gap` | `0.070` | PSID completed-children target; ACS co-resident `NCHILD` is not the right object |

The lifecycle ownership slope, ages 65--75 minus ages 25--34, and the
prime-age childless renter/owner median room levels are diagnostic only. They
are reported when available but are not hard SMM targets.

Tenure segmentation now allows small owner units while preserving the renter
cap:

- benchmark owner ladder:
  `H_own = [2.0, 4.0, 6.0, 8.0, 9.5, 11.0]`
- fast owner ladder:
  `H_own = [2.0, 4.0, 6.0, 11.0]`
- renter cap remains `hR_max = 8.0`

## Current Overnight Reduced-Target Runs

Launched: `2026-05-27 00:13 EDT`

Purpose: recalibrate after dropping the lifecycle ownership slope and
prime-age childless renter/owner median room levels from the hard SMM target
system. Local `main` and Torch scratch active calibration files were checksum
matched at commit `27c87f3`.

Verified target counts before launch:

- base targets: `15`
- direct targets: `17`, including `inv_pop_share_C` and
  `inv_rent_ratio_C_over_P`
- direct parameters: `15`

The current Python worker has one optimizer family, a seeded random/global plus
local adaptive proposal search. To get algorithmic variation without changing
model code, the overnight launch uses two proposal regimes per room
specification:

| Room spec | Regime | Run tag | Slurm waves |
|---|---|---|---|
| `hR_max=8.0` baseline | default proposal mix | `py_direct_reduced_targets_hR8_default_overnight_20260527` | `9670117`, then `9670118` |
| `hR_max=8.0` baseline | global-heavy proposal mix | `py_direct_reduced_targets_hR8_globalheavy_overnight_20260527` | `9670119`, then `9670120` |
| `hR_max=6.0` diagnostic | default proposal mix | `py_direct_reduced_targets_hR6_default_overnight_20260527` | `9670121`, then `9670122` |
| `hR_max=6.0` diagnostic | global-heavy proposal mix | `py_direct_reduced_targets_hR6_globalheavy_overnight_20260527` | `9670123`, then `9670124` |

Update after launch: the first 8 tasks from each cold-start first wave were
left running. The pending cold-start tasks `9--16` and second waves were
canceled after confirming that the saved incumbents were much better under the
reduced target system than the generic seed-bank starts. Incumbent losses
computed from saved moments under the reduced target system:

- `hR_max=8.0` incumbent: `16.097`
- `hR_max=6.0` incumbent: `35.230`

Additional incumbent-seeded arrays:

| Room spec | Regime | Run tag | Slurm job |
|---|---|---|---|
| `hR_max=8.0` baseline | default proposal mix | `py_direct_reduced_targets_hR8_incumbent_default_overnight_20260527` | `9670313` |
| `hR_max=6.0` diagnostic | default proposal mix | `py_direct_reduced_targets_hR6_incumbent_default_overnight_20260527` | `9670315` |
| `hR_max=8.0` baseline | tight local proposal mix | `py_direct_reduced_targets_hR8_incumbent_tight_overnight_20260527` | `9670412` |
| `hR_max=6.0` diagnostic | tight local proposal mix | `py_direct_reduced_targets_hR6_incumbent_tight_overnight_20260527` | `9670413` |

The first global-heavy incumbent arrays, `9670314` and `9670316`, were canceled
after the default/global-heavy starts accumulated roughly 60 evaluations per
run without improving the incumbents. The replacement tight local arrays use
`DT_DIRECT_GLOBAL_PROB=0.02`, `DT_DIRECT_INITIAL_SCALE=0.05`,
`DT_DIRECT_MIN_SCALE=0.003`, `DT_DIRECT_SHRINK=0.60`, and
`DT_DIRECT_STALL_WINDOW=6`.

Two-hour checkpoint, `2026-05-27 02:35 EDT`:

- best `hR_max=8.0`: run
  `py_direct_reduced_targets_hR8_incumbent_tight_overnight_20260527`,
  worker `7`, evaluation `295`, loss `10.697`, \(TFR=1.734\), ownership
  `0.337`, center share `0.467`, rent ratio `1.087`
- best `hR_max=6.0`: run
  `py_direct_reduced_targets_hR6_incumbent_default_overnight_20260527`,
  worker `3`, evaluation `274`, loss `15.474`, \(TFR=1.930\), ownership
  `0.622`, center share `0.432`, rent ratio `1.180`
- Status: both incumbent-local specifications improved relative to the saved
  incumbents, so the default and tight incumbent arrays were left running.

Second checkpoint, `2026-05-27 04:36 EDT`:

- best `hR_max=8.0`: run
  `py_direct_reduced_targets_hR8_incumbent_default_overnight_20260527`,
  worker `3`, evaluation `444`, loss `9.907`, \(TFR=1.662\), ownership
  `0.334`, center share `0.457`, rent ratio `1.066`
- best `hR_max=6.0`: run
  `py_direct_reduced_targets_hR6_incumbent_default_overnight_20260527`,
  worker `3`, evaluation `433`, loss `15.392`, \(TFR=1.959\), ownership
  `0.623`, center share `0.434`, rent ratio `1.181`
- First-wave tasks finished for the incumbent arrays; tasks `9--16` started
  and were left running.

Third checkpoint, `2026-05-27 06:37 EDT`:

- best `hR_max=8.0`: run
  `py_direct_reduced_targets_hR8_incumbent_default_overnight_20260527`,
  worker `10`, evaluation `289`, loss `9.323`, \(TFR=1.700\), ownership
  `0.358`, center share `0.456`, rent ratio `1.076`
- best `hR_max=6.0`: run
  `py_direct_reduced_targets_hR6_incumbent_default_overnight_20260527`,
  worker `14`, evaluation `268`, loss `15.329`, \(TFR=2.075\), ownership
  `0.606`, center share `0.483`, rent ratio `1.144`
- Because `hR_max=6.0` had largely flattened and the tight run was worse,
  `9670413` was canceled and replaced with a micro-local refinement seeded
  from the live `hR_max=6.0` best: run
  `py_direct_reduced_targets_hR6_microbest_overnight_20260527`, Slurm job
  `9676802`, with `DT_DIRECT_GLOBAL_PROB=0.00`,
  `DT_DIRECT_INITIAL_SCALE=0.015`, `DT_DIRECT_MIN_SCALE=0.0008`,
  `DT_DIRECT_SHRINK=0.50`, and `DT_DIRECT_STALL_WINDOW=4`.

Micro-local checkpoint, `2026-05-27 07:09 EDT`:

- best `hR_max=8.0` unchanged from the third checkpoint: loss `9.323`
- best `hR_max=6.0`: run
  `py_direct_reduced_targets_hR6_microbest_overnight_20260527`, worker `3`,
  evaluation `51`, loss `14.747`, \(TFR=1.995\), ownership `0.611`,
  center share `0.477`, rent ratio `1.131`
- Status: micro-local refinement improved the `hR_max=6.0` best and was left
  running.

Morning checkpoint, `2026-05-27 08:10 EDT`:

- best `hR_max=8.0`: run
  `py_direct_reduced_targets_hR8_incumbent_default_overnight_20260527`,
  worker `10`, evaluation `289`, loss `9.323`, \(TFR=1.700\), ownership
  `0.358`, center share `0.456`, rent ratio `1.076`
- best `hR_max=6.0`: run
  `py_direct_reduced_targets_hR6_microbest_overnight_20260527`, worker `3`,
  evaluation `131`, loss `14.743`, \(TFR=1.997\), ownership `0.612`,
  center share `0.477`, rent ratio `1.131`
- The hR8 default/tight and hR6 default incumbent arrays finished cleanly.
  The hR6 micro-local job `9676802` was still running.
- Current best JSON records were mirrored locally under
  `output/model/reduced_target_overnight_20260527/records/`.

Live pull, `2026-05-27 08:59 EDT`:

- best `hR_max=8.0` unchanged: loss `9.323`, \(TFR=1.700\), ownership
  `0.358`, center share `0.456`, rent ratio `1.076`
- best `hR_max=6.0`: run
  `py_direct_reduced_targets_hR6_microbest_overnight_20260527`, worker `3`,
  evaluation `306`, loss `14.730`, \(TFR=1.988\), ownership `0.613`,
  center share `0.478`, rent ratio `1.132`
- New summary fit note:
  `latex/current_fit_reduced_target_diagnostics_20260527.pdf`

Afternoon fast pulses, `2026-05-27`:

- Results and comparison table:
  `output/model/fast_calibration_pulses_20260527/summary.csv`
- Best default hR8 pulse record:
  `output/model/fast_calibration_pulses_20260527/records/hR8_standard_best.json`
- Best diagnostic high-owner/old-age record:
  `output/model/fast_calibration_pulses_20260527/records/hR8_chiowner_best.json`
- Best owner-ladder support diagnostic:
  `output/model/fast_calibration_pulses_20260527/records/Hown_lowtop_best.json`

Main results:

| Case | Loss | TFR | Own | Old-own | Old gap | Pop C | Rent ratio |
|---|---:|---:|---:|---:|---:|---:|---:|
| hR8 morning | `9.323` | `1.700` | `0.358` | `0.510` | `0.072` | `0.456` | `1.076` |
| hR8 standard swarm | `8.860` | `1.688` | `0.363` | `0.510` | `0.074` | `0.466` | `1.086` |
| Hown low-top exact | `8.854` | `1.688` | `0.363` | `0.511` | `0.073` | `0.466` | `1.086` |
| hR8 chi-owner | `11.806` | `1.755` | `0.413` | `0.738` | `0.066` | `0.449` | `1.120` |
| hR7 | `99.602` | `1.705` | `0.537` | `0.735` | `-0.119` | `0.485` | `1.037` |
| hR6 final | `14.716` | `1.983` | `0.614` | `0.999` | `-0.001` | `0.477` | `1.131` |

Read:

- A short hR8 incumbent-local swarm improved the default hR8 loss from
  `9.323` to `8.860`, but did not repair ownership; prime-age ownership
  remains around `0.36` and old-age ownership around `0.51`.
- hR6 and hR7 can force ownership mechanically, but they break gradients,
  parent-childless gaps, and/or fertility. They are diagnostics, not current
  benchmark candidates.
- Raising the effective owner-housing margin through the chi-owner basin can
  bring old-age ownership close to target (`0.738` vs `0.764`) and preserves
  the old-age parent-childless gap, but it raises TFR and still leaves
  prime-age ownership too low.
- Finer owner ladders around the 6--8 room range do not by themselves fix the
  room collapse. The low-top ladder gives a tiny objective improvement but
  leaves the median owner room choice at `8.0`.
- Parent down-payment waiver diagnostics did not repair the hR7/hR8 tradeoff;
  they should remain diagnostic unless given a structural interpretation.

Rush continuation, `2026-05-27 17:30--17:51 EDT`:

- Results:
  `output/model/rush_calibration_20260527/summary.csv`
- Best high-ownership candidate:
  `output/model/rush_calibration_20260527/records/hR8_chiowner_fertgrid_local.json`
- Best lower-TFR high-ownership candidate:
  `output/model/rush_calibration_20260527/records/hR8_chiowner_tradegrid2_local.json`

| Case | Loss | TFR | Own | Old-own | Old gap | Pop C | Rent ratio |
|---|---:|---:|---:|---:|---:|---:|---:|
| hR8 low-loss continue | `8.860` | `1.688` | `0.363` | `0.510` | `0.074` | `0.466` | `1.086` |
| hR8 chi-owner fertgrid local | `10.588` | `1.748` | `0.420` | `0.748` | `0.076` | `0.471` | `1.151` |
| hR8 chi-owner tradegrid local | `10.846` | `1.729` | `0.427` | `0.745` | `0.073` | `0.468` | `1.147` |
| hR8 expanded-fertility exact | `18.444` | `1.761` | `0.417` | `0.752` | `0.039` | `0.468` | `1.147` |

Read:

- The chi-owner basin can fit the old-age ownership block much better than the
  low-loss hR8 basin. The best rush candidate hits old-age ownership
  `0.748` versus target `0.764` and old-age parent-childless gap `0.076`
  versus target `0.070`.
- Prime-age ownership remains too low (`0.420--0.427` versus target `0.575`),
  and TFR remains above target (`1.729--1.748` versus `1.700`).
- Expanding fertility-cost bounds (`psi_child`, `c_bar_n`, `kappa_fert`) did
  not immediately improve the fit; exact high-cost probes worsened the loss.
  The high-owner basin is not simply blocked by the old fertility-cost bounds.
- The current practical choice is between the low-loss hR8 benchmark
  (`8.860`, bad ownership lifecycle) and the chi-owner diagnostic benchmark
  (`10.588`, much better lifecycle ownership but worse TFR/housing increments).

Room-discipline rush, `2026-05-27 18:00 EDT`:

- Code hook added for diagnostic-only extra targets:
  `DT_DIRECT_EXTRA_TARGETS="moment=target:weight,..."`.
- Code hook added for warm-start records:
  `DT_DIRECT_WARM_START_JSON=/path/to/best.json`.
- Warm starts now fail loudly if outside the active bounds; an initial
  `scout`-bound launch clipped the current best records and was canceled.
- Active room-discipline branches use `bounds=global`, `hR_max=8.0`, and soft
  room-median targets:
  `prime_childless_renter_median_rooms=4.0`,
  `prime_childless_owner_median_rooms=6.0`.

Active Torch jobs:

| Job | Run tag | Workers | Purpose |
|---|---|---:|---|
| `9720451` | `py_direct_roomsoft_hR8_lowloss_global_20260527_rush` | `8` active after pruning | soft room targets from low-loss hR8 |
| `9720452` | `py_direct_roomsoft_hR8_chiowner_global_20260527_rush` | `8` active after pruning | soft room targets from chi-owner hR8 |
| `9720595` | `py_direct_roomsoft_denseH_chiowner_global_20260527_rush` | `8` | dense owner ladder diagnostic, `H_own=[2,4,5.5,7,8.5,10]` |
| `9722184` | `py_direct_roomhard_noh12_chiowner_global_20260527_rush` | `8` | high room-target weights, diagnostic `housing_increment_1to2` weight `0.5` |

First checkpoints:

| Branch | Evals | Best loss | TFR | Own | Old-own | Renter med. | Owner med. | Read |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| low-loss soft | `312` | `11.120` | `1.687` | `0.364` | `0.511` | `6.093` | `8.000` | soft room targets barely move the low-loss basin |
| chi-owner soft | `286` | `12.512` | `1.746` | `0.420` | `0.750` | `5.896` | `8.000` | good old-age block, owner median still stuck at `8` |
| dense owner ladder | `176` | `30.939` | `1.764` | `0.324` | `0.819` | `6.029` | `7.000` | support can lower owner median, but at large objective cost |
| high room weights | `32` | `14.325` | `1.748` | `0.420` | `0.748` | `5.880` | `8.000` | first eval only; waiting for local proposals |

Current read:

- Aggregate housing demand discipline is not enough to discipline room
  heterogeneity. The price/shifter block can match aggregate demand and rents
  while the deterministic housing choice still bunches owner households on
  high rungs.
- Reintroducing room medians as soft targets has so far added loss without
  moving the owner median off `8.0` in the benchmark ladder.
- A denser owner ladder can mechanically move the owner median to `7.0`, but
  the current branch damages ownership and housing-response moments badly.
- Low-housing-need seed probes were canceled after poor first checks: they
  reduced renter rooms slightly but kept owner median at `8.0` and worsened
  fertility, geography, and ownership gradients.

Room-distribution sprint update, `2026-05-27 20:20 EDT`:

- Added diagnostic room-distribution moments:
  `owner25_45_rooms_le6_share`, `owner25_45_rooms_7to8_share`,
  `owner25_45_rooms_ge9_share`, renter cap shares, and renter/owner mean rooms
  by current-child state.
- Added a default-off owner-size wedge:
  `owner_size_cost * p_i * max(H_own - owner_size_cost_ref,0)^owner_size_cost_power`.
- Added a default-off diagnostic tenure/product logit:
  `tenure_choice_kappa`. It smooths the discrete choice over rent plus owner
  rungs. Conditional renter housing remains the continuous FOC object; the
  logit does not discretize renter rooms.
- Hard-argmax room-weight/lower-needs searches were launched widely on Torch.
  At peak, `240` short-partition workers were running under the effective
  per-user CPU cap. Current best room-diagnostic hard-argmax candidates still
  show a bad tradeoff:
  - `py_direct_roomwedge_hbarscale085_c012_20260527_rush`: loss `59.710`,
    \(TFR=2.018\), own `0.764`, old-own `0.997`,
    owner bins `<=6=0.613`, `7-8=0.387`, `>=9=0.000`,
    \(pop_C=0.276\).
  - `py_direct_widenet_H_ladder_midtail_20260527`: loss `63.888`,
    \(TFR=1.771\), own `0.435`, old-own `0.934`,
    owner bins `<=6=0.020`, `7-8=0.847`, `>=9=0.133`.
- Local capped-GE diagnostics with `tenure_choice_kappa` are encouraging:

| `tenure_choice_kappa` | Loss | TFR | Own | Old-own | Owner `<=6` | Owner `7-8` | Owner `>=9` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| `0.03` | `62.819` | `1.759` | `0.640` | `0.875` | `0.000` | `0.999` | `0.001` |
| `0.08` | `40.388` | `1.760` | `0.603` | `0.823` | `0.012` | `0.982` | `0.006` |
| `0.15` | `33.503` | `1.747` | `0.553` | `0.688` | `0.143` | `0.822` | `0.035` |
| `0.25` | `62.367` | `1.728` | `0.532` | `0.655` | `0.389` | `0.492` | `0.119` |
| `0.35` | `72.811` | `1.690` | `0.568` | `0.680` | `0.472` | `0.419` | `0.109` |
| `0.50` | `78.242` | `1.634` | `0.621` | `0.713` | `0.534` | `0.360` | `0.106` |

- These local diagnostics used `max_iter_eq=10`, `owner_h_bar_scale=0.90`,
  `owner_size_cost=0.006`, `H_own=[2,4,5.5,7,8.5,10]`, and the chi-owner
  warm start. Treat levels as directional until full GE runs are collected.
  The key result is that a small dispersion parameter can create real room
  heterogeneity without collapsing TFR immediately.
- Torch SSH authentication stopped accepting noninteractive commands after
  the wide hard-argmax launch (`BatchMode` returns permission denied). The
  dispersion code is compiled locally but has not yet been synced/launched on
  Torch. Refresh SSH/Kerberos before launching the next cluster wave.

Launch settings:

- setup: `benchmark`
- bounds: `global`
- closure: `outside_option_benchmark_normalized`
- geography weight: `100`
- array size: `1-16%8` per run, four first-wave arrays active together
- chained waves: second wave submitted with `afterany` dependency on the first
  wave and the same run tag, so workers resume through existing `best.json` and
  `evaluations.jsonl`
- per-worker budget: `13,500` seconds per wave; Slurm wall time `03:55:00`
- default proposal: `DT_DIRECT_GLOBAL_PROB=0.12`,
  `DT_DIRECT_INITIAL_SCALE=0.18`
- global-heavy proposal: `DT_DIRECT_GLOBAL_PROB=0.30`,
  `DT_DIRECT_INITIAL_SCALE=0.30`

Pre-launch checks:

- local and Torch target-count checks both returned `15` base targets,
  `17` direct targets, and `15` direct parameters
- active model and cluster launcher checksums matched between local `main` and
  Torch scratch
- one-evaluation Slurm smoke jobs completed for both room specs:
  `py_direct_reduced_targets_hR8_smoke_20260527` and
  `py_direct_reduced_targets_hR6_smoke_20260527`
- the obsolete old-target array `9662249` was canceled before launch

## Latest Cluster Search

Active household-head ownership / small-owner-ladder outside-option search:

- Slurm job: `9662249`
- Run tag: `py_direct_outside_headown_smallown_global_4h_20260526`
- Results directory:
  `/scratch/td2248/projects/Fertility_Spring26/code/cluster/results_python_direct_geometry_py_direct_outside_headown_smallown_global_4h_20260526`
- Purpose: first broad recalibration after correcting ownership targets to
  household heads/reference persons and allowing small owner units while
  preserving the renter cap.
- setup: `benchmark`
- bounds: `global`
- workers: `40` on Torch `cpu_short`, submitted as `1-40%32`
- internal worker budget: `13,500` seconds; Slurm wall time `03:55:00`
- population closure: `outside_option_benchmark_normalized`
- owner ladder: `H_own=[2.0,4.0,6.0,8.0,9.5,11.0]`
- renter cap: `hR_max=8.0`
- Pre-launch checks:
  - local one-evaluation direct-worker smoke completed
  - `check_population_closure.py` passed
  - Torch two-worker Slurm smoke job `9662202` completed and wrote
    `config.json`, `evaluations.jsonl`, `best.json`, and `status.json`
- Initial production check: `32` worker directories wrote configs/status/best
  records within the first minute; tasks `33--40` were pending on the array
  concurrency limit.
- On 2026-05-26 evening, while this run was still active, the array throttle
  was lowered to `24` and the eight weakest active workers by current best loss
  were canceled to free CPU slots for the `hR_max=6` diagnostic below. Their
  partial JSON outputs remain in the results directory.
- Snapshot used for the current fit note, pulled at `2026-05-26 22:45 EDT`:
  worker `29`, evaluation `440`, loss `19.662`, \(TFR=1.772\),
  ownership `0.313`, childless-renter median rooms `5.974`, center share
  `0.449`, and rent ratio `1.064`. The run was still active when this
  snapshot was pulled.
- Current fit note:
  `latex/current_fit_slide_diagnostics_20260526.pdf`; generated source is
  `code/model/tools/make_direct_fit_slide_note.py`.

Parallel renter-cap diagnostic:

- Slurm smoke job: `9664816`
- Slurm diagnostic job: `9664883`
- Run tag: `py_direct_outside_headown_smallown_hr6_30m_20260526`
- Results directory:
  `/scratch/td2248/projects/Fertility_Spring26/code/cluster/results_python_direct_geometry_py_direct_outside_headown_smallown_hr6_30m_20260526`
- Purpose: test whether lowering only the renter cap from `hR_max=8.0` to
  `hR_max=6.0` helps the childless-renter room moment and prime-age ownership
  without treating it as a benchmark change.
- setup: `benchmark`
- bounds: `global`
- requested workers: `16` on Torch `cpu_short`, submitted as `1-16%16`
- active split at launch: `24` running workers remain on the `hR_max=8`
  production run and `8` workers started on the `hR_max=6` diagnostic; tasks
  `9--16` were pending on the user CPU limit.
- internal worker budget: `1,800` seconds; Slurm wall time `00:40:00`
- population closure: `outside_option_benchmark_normalized`
- owner ladder unchanged:
  `H_own=[2.0,4.0,6.0,8.0,9.5,11.0]`
- renter cap override: `hR_max=6.0`
- Smoke job `9664816` completed and verified written configs with
  `hR_max=6.0`.
- Final collection succeeded with `16` worker directories. Best diagnostic
  point: worker `12`, evaluation `59`, loss `55.106`, \(TFR=1.931\),
  ownership `0.629`, childless-renter median rooms `6.000`, old-age ownership
  `0.985`, center share `0.413`, and rent ratio `1.033`.
- Read: `hR_max=6` mechanically binds the childless-renter room moment but
  looks too distortionary relative to the live `hR_max=8` search: ownership is
  pushed high, late-life ownership is nearly universal, and the lifecycle
  ownership slope is too steep.

Corrected DUE-common-support diagnostic pulse:

- Slurm job: `9546318`
- Run tag: `py_direct_outside_commonH_owner0_tfrlt2_20m_20260525_003535`
- Results directory:
  `/scratch/td2248/projects/Fertility_Spring26/code/cluster/results_python_direct_geometry_py_direct_outside_commonH_owner0_tfrlt2_20m_20260525_003535`
- Purpose: diagnostic run with the owner ladder lowered to zero and rental
  access extended to the owner maximum:
  `H_own=[0.0, 2.2, 4.4, 6.6, 8.8, 11.0]`,
  `hR_max=max(H_own)=11.0`.
- Existing property tax `tau_H=0.01`, transaction cost `psi=0.06`, and
  financed share `phi=0.80` remained active.
- setup: `benchmark`
- bounds: `global`
- workers: `32` on Torch `cpu_short`
- internal worker budget: `1,200` seconds; Slurm wall time `00:30:00`
- population closure: `outside_option_benchmark_normalized`
- launcher city-entry probability: `q^{E,*}=0.89`
- hard TFR cap: candidates with \(TFR\ge 2.0\) receive loss `1e6`.
- Smoke job `9546293` completed first with
  `DT_DIRECT_H_OWN_MIN=0` and `DT_DIRECT_HR_MAX=owner_max`; configs verified
  the lowered owner ladder and `hR_max=11.0`.
- Final collection succeeded with `32` worker directories, `29--49`
  evaluations per worker, and empty `squeue` for job `9546318`.
- Partial best: worker `26`, evaluation `24`, loss `82.2035`,
  \(TFR=1.893\), childless rate `0.074`, mean age first birth `29.68`,
  ownership `0.194`, childless-renter median rooms `6.872`,
  childless-owner median rooms `6.6`, center share `0.437`, rent ratio `1.538`.
- Read: the corrected common-support/owner-floor-zero diagnostic still does
  not repair the small-renter-unit miss. Lowering the owner floor reduces the
  owner-room median but leaves childless renters around `6.9` rooms and
  ownership far below target.

Earlier common-support diagnostic with owner grid left fixed:

- Slurm job: `9543538`
- Run tag: `py_direct_outside_commonH_tfrlt2_global_6h_20260525_000050`
- Results directory:
  `/scratch/td2248/projects/Fertility_Spring26/code/cluster/results_python_direct_geometry_py_direct_outside_commonH_tfrlt2_global_6h_20260525_000050`
- Purpose: diagnostic run with perfect overlap between the renter support and
  owner support, leaving the existing owner grid fixed:
  `H_own=[4.0, 5.4, 6.8, 8.2, 9.6, 11.0]`,
  `hR_max=max(H_own)=11.0`.
- This does not add a new economic mechanism. The existing property tax
  `tau_H=0.01`, transaction cost `psi=0.06`, and financed share
  `phi=0.80` remain active.
- setup: `benchmark`
- bounds: `global`
- workers: `32` on Torch `cpu_short`
- wall time: `6:00:00`
- internal worker budget: `20,700` seconds
- population closure: `outside_option_benchmark_normalized`
- launcher city-entry probability: `q^{E,*}=0.89`
- hard TFR cap: candidates with \(TFR\ge 2.0\) receive loss `1e6`.
- Cluster smoke job `9543502` completed first with
  `DT_DIRECT_HR_MAX=owner_max`; configs verified `hR_max=11.0`.
- Stopped manually with `scancel` after about 20 minutes total runtime, by
  design. This was a diagnostic pulse, not a production calibration.
- Final partial collection succeeded with `32` worker directories, `15--36`
  evaluations per worker, and empty `squeue` for job `9543538`.
- Partial best: worker `7`, evaluation `13`, loss `83.5964`, \(TFR=1.902\),
  ownership `0.441`, childless-renter median rooms `6.700`,
  childless-owner median rooms `8.2`, center share `0.223`, rent ratio `1.143`.
- Read: common rental/owner support does not mechanically repair the
  small-rental-unit miss. Even with `hR_max=11.0`, the best diagnostic point
  keeps childless renters around `6.7` rooms and has much weaker ownership and
  geography than the previous `hR_max=8` outside-option search.

Active outside-option recalibration:

- Slurm job: `9516170`
- Run tag: `py_direct_outside_sout169_tfrlt2_global_6h_20260524_151347`
- Results directory:
  `/scratch/td2248/projects/Fertility_Spring26/code/cluster/results_python_direct_geometry_py_direct_outside_sout169_tfrlt2_global_6h_20260524_151347`
- setup: `benchmark`
- bounds: `global`
- workers: `32` on Torch `cpu_short`
- wall time: `6:00:00`
- internal worker budget: `20,700` seconds
- population closure: `outside_option_benchmark_normalized`
- empirical outside-origin entrant-share normalization:
  \(s^{E,\mathrm{out}}=0.169\)
- launcher city-entry probability: \(q^{E,*}=0.89\)
- hard TFR cap: candidates with \(TFR\ge 2.0\) receive loss `1e6` and cannot
  be selected as best.
- first collection succeeded with `32` worker directories; initial best was
  worker `1`, evaluation `1`, loss `353.731`, \(TFR=0.891\),
  ownership `0.121`, \(N=1.000\), \(pop_C=0.393\), rent ratio `1.276`.
- The prior uncapped run
  `py_direct_outside_sout169_global_6h_20260524_145525` was canceled after the
  best-so-far had \(TFR>2\); use the capped run above for live results.

The latest completed cluster search below was run under the prior
`renewal_valve_calibrated` closure. It is useful as a pre-outside-option
benchmark, but it is not a recalibration of the current live closure.

Results directory:

- `code/cluster/results_python_direct_geometry_py_direct_renewal_calibrated_global_12h_20260506`

Torch arrays:

- `8113159`, `8113160`, `8113161`
- `120` total tasks on `cpu_short`
- all tasks completed with exit code `0:0`

Search settings:

- setup: `benchmark`
- bounds: `global`
- parameters: `15`
- `geo_weight = 100`
- `renewal_retention = 1.0`
- `scale_target = 1.0` imposed mechanically
- per-stage time budget: `13,500` seconds

Best collected candidate:

- worker `37`, evaluation `1425`
- loss `8.1275`
- solve time `44.2` seconds
- GE accepted, but not strict: final equilibrium error `0.003006`,
  convergence reason `soft_tol_10x`
- prices: \(p_P=0.49766\), \(p_C=0.58848\)

Closure diagnostics:

| Diagnostic | Value |
|---|---:|
| `entry_per_unit_scale` | `0.0166667` |
| `mature_cityborn_per_unit_scale` | `0.0155334` |
| `outside_entry_flow` | `0.0011333` |
| `population_scale_denominator` | `0.0011333` |
| `implied_total_population` | `1.000` |
| `scale_factor` | `1.000` |
| `stationary_entry_relative_residual` | `0` |

## Target Fit

Targets come from `build_calibration_setup()` in the active Python
`parameters.py`. The model column below is from the latest completed
pre-ownership-correction renewal-valve run, so it is a historical benchmark
until a new outside-option calibration is collected. The last two rows are
direct-geometry targets disciplined by the geometry-weight block.

| Moment | Target | Model |
|---|---:|---:|
| `tfr` | `1.700` | `1.898` |
| `childless_rate` | `0.150` | `0.145` |
| `mean_age_first_birth` | `26.000` | `33.535` |
| `tfr_gradient` | `0.133` | `0.119` |
| `own_rate` | `0.575` | `0.643` |
| `own_gradient` | `0.186` | `0.139` |
| `own_family_gap` | `0.168` | `0.114` |
| `housing_increment_0to1` | `0.664` | `0.441` |
| `housing_increment_1to2` | `0.566` | `0.192` |
| `young_liquid_wealth_to_income` | `0.600` | `0.527` |
| `center_share_nonparents` | `0.494` | `0.405` |
| `center_share_newparents` | `0.416` | `0.382` |
| `migration_rate` | `0.032` | `0.035` |
| `old_age_own_rate` | `0.764` | `0.947` |
| `old_age_parent_childless_gap` | `0.070` | `0.062` |
| `inv_pop_share_C` | `0.450` | `0.441` |
| `inv_rent_ratio_C_over_P` | `1.140` | `1.182` |

Current read:

- Under the prior renewal-valve search, benchmark scale was imposed
  mechanically, the outside flow was positive, and the stationary entry
  residual was zero.
- Remaining misses are concentrated in first-birth timing, childless renter
  baseline rooms, the first- and second-birth housing responses, and old-age
  ownership.
- Do not compare this loss mechanically to pre-Python, pre-renewal-valve, or
  MATLAB-inversion losses.

## Reporting Rules

When reporting calibration results, always show targets next to model moments.
Do not report only losses or only model moments when target values are
available.

Do not treat `parity_progression_1to2` as a calibration target unless the
fertility architecture has been changed to support a sequential second-birth
hazard.

Use `CALIBRATION_STATUS.md` first. Historical notes are archived here:

- `docs/archive/calibration_docs_2026-05-07/`
- `docs/archive/root_notes_2026-05-07/`
- `calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/status_notes/`
