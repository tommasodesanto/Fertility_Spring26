# HANDOFF: update the paper draft with the July 16–17 decisions

Target file: `latex/intergenerational_housing_fertility_paper_draft.tex`.
Rules: minimal-change edits only — do not rewrite surrounding prose, do not
touch unrelated sections. Follow `docs/style/econ_writing_style_guide.md`.
The word "parity" is banned everywhere; write "number of children" or use
(n,s) notation. Do not compare losses across different target systems.

## 1. Current working calibration = M5 (use these numbers, no others)

Strict winner, 15 moments / 14 estimated parameters, Nb=120 (the project
standard), loss `9.044`, established-12 block `3.166` (the project's best),
equilibrium residual 1.5e-5, bit-identical tight repeats. Do not cite M3 or
M4 numbers as current anywhere.

Estimates: beta_annual 0.9912, alpha_cons 0.5912, c_bar_0 1.2596,
c_bar_n 0.3929, h_bar_0 0.3915, h_bar_jump 1.6021, h_bar_n 0.1644,
psi_child 0.1979, kappa_fert 2.0434, chi 1.1133, H0 8.6456,
theta0 0.3118, theta1 0.3973, tenure_choice_kappa 0.0100.
External restriction: theta_n = 0.

Target fit (target / model): tfr 1.918/2.034; childless 0.188/0.189;
own_rate 0.575/0.658; own_family_gap 0.168/0.219; housing_increment_0to1
0.664/0.681; renter mean rooms 3.805/3.963; owner share rooms>=6
0.596/0.760; young liquid wealth 0.179/0.328; owner-renter rooms
2.419/2.323; own_rate_2534 0.341/0.274; 3plus-vs-1to2 rooms 0.368/0.607;
estate median 76-84 6.501/6.431; nonhousing>=1x-income share 65-75
0.608/0.610; old_age_own_rate 0.764/0.954; aggregate rooms 5.780/5.673.

Caveats that MUST accompany any use of M5:
- theta1 is weakly identified; the identification audit is in progress. Do
  not describe theta1 as precisely identified; "estimated, with a flat
  objective direction" is the honest phrasing.
- Old-age ownership (0.954 vs 0.764) is a disclosed structural miss: the
  model has no late-life exit margin (no health/LTC expense risk, no
  maintenance burden), so owners hold housing until death. Make no claims
  that rest on the late-life ownership path or late-life portfolio
  composition (untargeted: estate p90/p50 1.75 vs 3.45; 76-84 nonhousing
  quantiles 0 vs PSID 2.4/10.3/26.3).
- tenure_choice_kappa estimated interior (0.0100) for the first time — its
  identifying moment (old-age ownership) is now in the target system. One
  sentence of economics: tenure choice has a small idiosyncratic-taste
  component; the ownership age profile identifies its scale.

## 2. Income process — REPLACE the circular sentence (currently ~lines 482–486)

Delete the claim that the annual AR(1) "preserv[es] the stationary
log-earnings dispersion of the annual process" (it is circular). Insert, in
the calibration section only (not the introduction):

> Household income risk follows an annual AR(1) with persistence 0.960,
> inside the 0.88–0.96 range PSID studies report (Card 1994;
> Hubbard–Skinner–Zeldes 1995; Heathcote–Storesletten–Violante 2010, as
> surveyed in Sommer and Sullivan 2018), sampled at the model's four-year
> frequency and discretized by a five-state Rouwenhorst chain. The
> innovation variance is set below econometric estimates of gross income
> risk. The reason is feasibility at the bottom of the distribution: the
> model requires every household to afford a subsistence consumption and
> housing bundle, and it contains no means-tested safety net, family
> transfers, or labor-supply margin. Under measured income dispersion, a
> household holding the lowest persistent income state could not afford the
> bundle in some states, which no feasible choice can resolve. The
> calibrated variance should therefore be read as the effective,
> post-insurance risk households bear rather than gross income risk
> (Blundell, Pistaferri and Preston 2008 estimate that roughly one third of
> permanent income shocks are insured; De Nardi 2004 similarly feeds
> transfer-inclusive income into a model without an explicit safety net).
> We treat the innovation variance as provisional; disciplining it
> externally — through a measured safety-net bound on the lowest income
> state or a preference specification without subsistence intercepts — is
> characterized in a companion memo and left to future revision.

(If the draft cites Guvenen–Smith 2014 elsewhere it may be added to the
insurance sentence; do not add it otherwise without verifying the citation.)

## 3. Bequest block — corrections and language

- Specification: child-blind De Nardi (2004) luxury warm glow over the
  estate W = max{b + pH, 0}, zero-normalized; theta0 and theta1 estimated;
  theta_n = 0 externally.
- theta_n = 0 citations (all verified from the papers): Hurd (1989,
  Econometrica) — elderly with children decumulate no differently;
  Kopczuk and Lupton (2007, REStud) — children affect motive incidence,
  not scale; De Nardi, French, Jones and McGee (2025) — children not a
  conditioning variable in retirement saving; Kvaerner (2023, RFS 36(8),
  3382–3422, "How Large Are Bequest Motives? Estimates Based on Health
  Shocks" — CHECK the bibliography: an earlier note carried a wrong title).
  Our own PSID family-size estate gap is 0.101 (bootstrap SE 0.563) —
  statistically zero, and a marital-composition cancellation (−0.886 among
  married, +0.884 among nonmarried households).
- The "children and bequests" defense, if a referee-facing sentence is
  wanted: children direct resources through housing while parents are
  alive — they raise required space, parents buy larger owner-occupied
  homes, transaction costs make those homes sticky into old age, and
  housing is the dominant component of estates — so the model generates
  larger estates for parents with no child-directed bequest taste.
- CORRECT anywhere they appear: Nakajima–Telyukova (2017, JF, Table I
  Panel B) bequest parameters are gamma = 20.534 and zeta = $7,619 per year
  in 2000 dollars (approximately one quarter of their retired homeowners'
  median after-tax income). The values gamma = 0.43 / zeta ≈ $19,600 that
  circulated in earlier project notes appear in NO version of that paper.

## 4. Target-system story (if/where the draft describes calibration targets)

- Hard wealth/ownership targets now: the 76–84 median total-estate-to-income
  ratio (6.501, person-bootstrap weight); the 65–75 share of households with
  nonhousing wealth of at least one year of income (0.608, person-bootstrap
  weight); the 65–75 ownership rate (0.764, ACS). Twelve established
  fertility/ownership/rooms/liquid-wealth moments unchanged.
- The estate p90/p50 ratio and the family-size estate gap are DIAGNOSTICS,
  not targets. Reasons, if stated: the PSID upper tail is dominated by
  business, IRA and other financial wealth with no model counterpart (top
  estate decile housing share 0.186 in PSID); the family-size gap is
  statistically zero and marital-composition-driven, and the model has no
  marital state.
- A nonhousing wealth MEDIAN is not targeted because the model's liquid
  wealth lives on a grid: median-type moments are locally flat (zero
  derivative), share-type moments move smoothly. One sentence at most; per
  the style rules, do not narrate numerical conventions beyond that.
- Entrant wealth: entrants at age 18 draw from the PSID net-worth-to-income
  distribution of childless renters aged 18–24 (previously an age 25–35
  object whose mean was also a target — that circularity is resolved; do
  not mention the old construction unless the draft already does).

## 5. Do NOT write into the draft

- Any suggestion that the income variance is estimated or externally
  disciplined (it is provisional, per the Section 2 note).
- The equivalence-scale / Stone–Geary revision (under consideration, not
  decided), the safety-net transfer mechanism (characterized in a memo,
  not implemented), or any M6 promises.
- Claims resting on the late-life ownership path, late-life portfolio
  composition, or the estate upper tail.
- Loss comparisons across target systems, or any M3/M4-era numbers.

## 6. File pointers (for your reference, not for citation in the draft)

- M5 report: `output/model/intergen_income_disciplined_recalibration_20260716/report/`
- Income decision memo: `docs/model/intergen_income_risk_feasibility_decision_memo_20260717.md`
- M4 audit (history of the corrections above): `docs/model/intergen_m4_calibration_audit_20260716.md`
- Bequest review: `docs/model/intergen_bequest_balance_sheet_fable_review_20260716.md`
