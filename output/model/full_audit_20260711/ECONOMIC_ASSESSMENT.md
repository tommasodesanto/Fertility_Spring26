# Economic assessment of the 14.780 current-bound candidate (lead, 2026-07-11)

Scope: interpretation of the five persistent misses, the at-bound parameters, and the
relaxed-housing diagnostic, at the July 10 combined specification (Rouwenhorst income,
normalized bequests, r=2%/\delta=1.1% annual, \eta=1.75, H0-identified rooms target,
15 moments / 14 parameters). Numbers are from the extracted report and the lead
reproduction battery. Not judged by loss alone.

## 1. The tenure-age tension (old-age ownership vs young ownership vs aggregate)

The model delivers own(25-34)=0.251 vs 0.341, own(all)=0.653 vs 0.575,
own(65-75)=0.948 vs 0.764. One deterministic tenure margin (tenure_choice_kappa=0)
with a uniform down payment (1-\phi)=0.20 and no old-age exit margin can only
generate an ownership age profile that rises late and never falls. Households buy
at or near the fertility deadline (the July 6 mechanics memo located the jump at
the last fertile age), overshooting old ownership by 18 points while undershooting
young ownership by 9. The old-age miss carries weight 160 and alone contributes
5.39 of 14.78 (36%). The economics are clear: nothing in the model makes an
80-year-old sell — no mortality risk before J, no health/LTC shocks, no downsizing
motive beyond the bequest calculus, and the bequest function now actively rewards
holding the gross estate b+pH. In the data the 65-75 rate is already past its peak.
This is the D11 missing-exit-margin diagnosis, now SHARPENED: in the June record
old_age_own_rate loaded on box-pinned beta/chi/theta0, so the miss could be read as
bound-limited; in the July 10 candidate beta (0.9577 annual) and theta0 (0.132) are
INTERIOR. The pull is no longer against the box — the model genuinely cannot lower
old ownership at any admissible preference vector. The 0.18 gap is the shadow price
of the missing mechanism, not of the bounds.

## 2. The wealth-composition tension (young liquid vs old liquid vs housing)

young_childless_renter liquid/income = 0.543 vs 0.179 (3.0x);
old nonhousing median/income = 1.029 vs 2.230 (0.46x). The model cannot
simultaneously (i) accumulate enough liquid wealth by 25-34 to clear the (1-\phi)pH
down payment on family-sized rungs (needed for own_rate_2534 = 0.341), (ii) show
young renters holding only 0.18 of annual income, and (iii) have retirees hold 2.2x
income in nonhousing wealth while 95% of them hold the housing asset. One
discount factor and one riskless asset serve all three. The estimator's compromise:
beta_annual=0.9577 (interior — patience up from the old 0.94 corner once bequests
switched on), c_bar_0 pinned at its 1.28 cap (max subsistence burn to damp young
saving), and the down-payment saving motive still forces the young ratio to 3x its
target. Note the direction change: the July 9 second opinion predicted the age-18-21
income fix would RELIEVE the young-wealth overshoot (from 0.353); instead the
combined candidate sits at 0.543. The income fix lowered young income (denominator
of the ratio) while the saving need (numerator driver: down payments at p=0.69,
family rungs H=6,8) is unchanged — the ratio worsened. The miss is structural:
either young households in the model need a reason not to save (higher discounting
of the young alone, entry debt, rental insurance) or the ownership block needs a
non-wealth channel into young ownership (credit access, parental transfers — the
D7 design). The old-liquid miss is the mirror image: old wealth is warehoused in
housing because nothing converts it back to the liquid asset before death.

## 3. The housing-quantity tension (renter rooms and the owner-renter gap)

renter mean rooms 4.281 vs 3.805 (+0.48, contrib 1.36); owner-renter gap 1.948 vs
2.419 (-0.47, contrib 2.66); owner share>=6 rooms 0.741 vs 0.596 (+0.15). The gap
miss is mostly a RENTER-side failure: renters consume too much housing. With the
renter cap at 6 and the calibrated floors (h_bar_0 at its lower bound 1.0), renter
demand piles up high in the 1-6 band; the historical cap-bunching diagnosis (~0.73
of the fertility derivative cap-pinned) lives in the same geometry. The search
wants h_bar_0 lower still (bound) and chi higher (bound at 1.15): both would widen
the owner-renter service wedge. chi is doing double duty — it is the tenure lever
(owner premium raises ownership everywhere, feeding the old-age overshoot) AND the
room-gap lever; the two roles conflict, which is the chi ridge the July 8 Jacobian
found. The relaxed-housing arm shows exactly this: chi 1.20 buys the room gap
(contrib 2.66 -> 0.64) and young tenure (own_2534 hits) but pushes aggregate and
old ownership further up and empties old liquid wealth. One parameter cannot price
owner services and ration tenure at the same time.

## 4. At-bound parameters — adjudication

- H0 = 9.9997 at the [1,10] upper bound: BOUND ERROR, not economics. At p_eq=0.6898
  supply is H0*(p*ucr/r_bar)^1.75 = 5.55; the rooms target 5.78 needs H0 ~ 10.4 at
  this price. The box mechanically caps the moment H0 exists to match. The bound is
  arbitrary (H0 units are supply-scale, not rooms); widen to ~[1,20]. Costless fix,
  re-polish required; expect the rooms contribution (0.31) to vanish and small GE
  reallocation elsewhere.
- c_bar_0 = 1.28 at the upper cap (0.32/yr x 4): serves as the young-wealth damper
  (Section 2). Whether the data admit a higher subsistence level is a calibration
  judgment; the preference-relaxed arm allowed 1.80 but its best was NOT extracted
  locally (blocked: cluster). Treat as weak identification + missing mechanism, not
  as a defensible interior optimum.
- chi = 1.15 at the upper bound: the bound is the binding constraint on the
  room-gap block (relaxed arm: chi -> 1.20 interior, loss -0.88). But the chi that
  fixes rooms breaks tenure (Section 3). The bound is not obviously wrong — chi>1.15
  drives ownership toward universal. The real issue is a missing tenure-rationing
  mechanism, not the box.
- h_bar_0 = 1.0 at the lower bound: same block as chi (relaxed arm: 0.92 interior).
  A childless housing floor below one room is hard to defend economically; the
  search pushing there says the renter-side quantity fit wants something the floor
  architecture cannot give.
- tenure_choice_kappa = 0 at the lower bound: user-imposed deterministic tenure
  (Frechet smoothing rejected). This is an external restriction sitting inside the
  search box; it should be REMOVED from the search vector and fixed at 0, which
  also repairs the parameter count honestly (13 free parameters, 15 moments).
- beta_annual = 0.9577 INTERIOR: for the first time the 0.94 floor is slack. The
  bequest activation (theta0 0.132, theta_n 0.768, both interior) resolved the old
  beta/theta0 corner ridge and D11's theta_n non-identification concern in its old
  form (theta0 is no longer ~0, so theta_n has a live gradient). Fresh local
  Jacobian at the candidate: see diagnostics/jacobian_current_bound.json.

## 5. The relaxed-housing diagnostic — informative, misleading, do not adopt

Contribution ledger (report summary): improves room gap -2.02, renter rooms -1.03,
young wealth -0.83, own_2534 -0.65; worsens old liquid median +1.59 (model value
0.379 vs 2.230 — retirees essentially liquid-broke), old own +0.93, own_rate +0.86,
tfr +0.56. Net -0.88. The scalar improvement rides on the weight asymmetry (room
gap w=12 vs old-wealth median w=0.8). Under any weighting that respects sampling
precision (the old-wealth median SE is not 4x the room-gap SE), the trade reverses.
Verdict: the relaxed arm correctly measures the shadow price of the chi/h_bar_0
bounds on the room block; as a candidate it is economically inadmissible
(near-universal ownership, empty old portfolios, TFR 2.10). Keep as diagnostic.

## 6. Identification (15 moments, 14 parameters)

Count is adequate only nominally. Effective concerns: (i) tenure_choice_kappa at
its externally-preferred corner contributes no identifying variation — treat as
fixed, 13 free; (ii) H0 at a mis-set bound removes its gradient direction;
(iii) the chi ridge (chi, h_bar_0, h_bar_jump, h_bar_n jointly move the room block
against the tenure block); (iv) the young-wealth moment disciplines c_bar_0 only
through the cap. The June 18 Jacobian (rank 13, cond 2.7e4) predates the combined
spec; the audit Jacobian at the new candidate (14 params incl. H0, tight
equilibrium) quantifies the current effective rank — see
diagnostics/jacobian_current_bound.json when complete. Weights remain ad hoc
(documented hand weights); code/data/moment_standard_errors/ (new, untracked)
suggests SE work exists but is not yet wired into the objective — until then,
influence is concentrated (old_age_own_rate alone carries 36% of the loss) and
cross-candidate comparisons inherit that arbitrariness.

## 7. What would move the economics (for discussion only — no changes made)

Ranked by expected value per unit of new complexity, given the ledger history
(D3/D5/D7/D11/D12 all prepared): (1) an old-age exit/downsizing margin (mortality
before J with estate realization, or a forced-sale/health shock) — attacks the
single largest contribution (5.39) and the old-liquid miss simultaneously;
(2) the D7 parental/down-payment transfer — attacks young ownership without
inflating young saving (its design memo exists; pi_T, tau_p external); (3) H0
bound widening — free; (4) SE-based weights — turns the scalar loss into an
interpretable statistic before any further search is bought. A renter-cap /
floor-architecture rework (D3/D5) remains the deeper fix for the renter-rooms
block but is a model change of a different order.
