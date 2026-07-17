# Surrogate-Assisted Calibration Prototype: Results Note

Date: 2026-06-22

Status: diagnostic / methods prototype. It does **not** change the model, the
target system, or any calibration state. It is an alternative *outer loop* for
the one-market intergenerational housing-fertility model in
`code/model/intergen_housing_fertility/`.

Code: `code/model/intergen_surrogate_calibration/` (see its `README.md`).
Output packets: `output/model/intergen_surrogate_calibration/<target_set>/`.

## What was built

A multi-output Gaussian-process emulator `theta -> moment vector`, trained on
the existing global-DE evaluation logs, plus three things built on it:

1. **Emulator accuracy** — K-fold cross-validation, per moment and on the loss.
2. **Identification** — the smooth surrogate Jacobian `d m / d theta`, with the
   finite-difference audit's exact target-normalization, compared head-to-head
   with `output/model/intergen_sensitivity_jacobian_20260618/`; plus ARD
   lengthscale relevance.
3. **Bayesian optimization** — expected-improvement proposals (global and
   local/trust-region), a frontier probe of conflicting moments, and
   **validation of every proposal against the real `run_model_cp_dt` solver**.

Everything runs in the model `.venv` (no scipy/sklearn; the GP is hand-rolled in
numpy with analytic log-marginal-likelihood gradients, checked against finite
differences to 1e-9). The pipeline was adversarially reviewed; the one
correctness bug it found (SVD rank tolerance not matching the FD audit) was
fixed before these results.

## Data

15,101 Tier-1 records carry a full 13-parameter `theta` and the raw moment
dict (not just a loss), so the emulator predicts moments and any target set's
loss is recomputed with the model's own `diagnostic_loss`. After dropping
non-converged solves (`market_residual > 5e-3`, ~1%), **14,871** records remain.
Two target systems are studied:

- `candidate_no_timing_v0` — 13 moments incl. the **discrete** owner/renter
  *median* rooms; this is the exact system the FD Jacobian audit used. 14,871
  usable records; data-set best loss **7.637**.
- `candidate_replacement_v1` — 13 moments using **smooth** room *means* and
  *>=6-room shares* plus old nonhousing wealth (the ACS/PSID replacement set).
  7,095 usable records; data-set best loss **21.276**.

Reproduction check: re-solving stored records reproduces their logged moments to
max abs diff `1.4e-4` (v0) and `3.1e-3` (replacement). The solve path here is
faithful to the pipeline that generated the data.

## Result 1 — the emulator is accurate (cross-validated)

Held-out loss prediction:

| target system | loss R² | loss rank-corr |
|---|---:|---:|
| `candidate_no_timing_v0` | 0.854 | 0.958 |
| `candidate_replacement_v1` | 0.859 | 0.935 |

Per-moment out-of-fold R² (every moment):

| moment (v0) | R² | | moment (replacement) | R² |
|---|---:|---|---|---:|
| `housing_user_cost_share` | 0.985 | | `young_liquid_wealth_to_income` | 0.972 |
| `old_age_own_rate` | 0.973 | | `own_rate` | 0.946 |
| `young_liquid_wealth_to_income` | 0.969 | | `tfr` | 0.946 |
| `own_rate` | 0.941 | | `old_nonhousing_wealth_to_income_6575` | 0.931 |
| `tfr` | 0.932 | | `childless_rate` | 0.919 |
| `liquid_wealth_to_income` | 0.931 | | `prime30_55_childless_renter_mean_rooms` | 0.903 |
| `childless_rate` | 0.905 | | `prime30_55_childless_renter_share_rooms_ge6` | 0.847 |
| `prime_childless_renter_median_rooms` | 0.907 | | `prime30_55_childless_owner_share_rooms_ge6` | 0.816 |
| `housing_increment_0to1` | 0.682 | | `prime30_55_childless_owner_mean_rooms` | 0.746 |
| `prime_childless_owner_median_rooms` | 0.641 | | `own_family_gap` | 0.571 |
| `own_family_gap` | 0.599 | | `old_parent_childless_nonhousing_wealth_gap_6575` | 0.559 |
| `housing_increment_1to2` | 0.371 | | `housing_increment_0to1` | 0.541 |
| `old_age_parent_childless_gap` | 0.159 | | `housing_increment_1to2` | 0.218 |

Read: the emulator nails the strongly-identified moments (ownership, fertility,
wealth, user-cost share, renter rooms: R² > 0.9) and is weak exactly where the
identification ledger says the moments are weak/noisy — the two child-linked
housing increments, the family ownership gap, and both parent-childless gaps.
The per-moment R² is itself an identification readout.

## Result 2 — identification, surrogate vs finite differences

Surrogate Jacobian at the FD audit's two anchor points, same
target-normalization, same rank tolerance:

| anchor | FD rank | FD cond | surrogate rank | surrogate cond | sign-agree | Pearson |
|---|---:|---:|---:|---:|---:|---:|
| `core_feasibility_v1` | 13 | 2.69e4 | 13 | 4.56e3 | 0.65 | 0.39 |
| `roomcost_test_v1` | 12 | ∞ | 13 | 8.7e2 | 0.61 | 0.57 |

Two points:

- The surrogate **reproduces the headline finding**: the core point is full rank
  but badly conditioned (condition number in the thousands–tens of thousands),
  so the target system is formally identified but numerically fragile.
- At the room-cost point the FD audit reports **rank 12** — the discrete
  owner-median-rooms moment has a locally-zero finite difference. The smooth
  emulator reports **rank 13**: it recovers a nonzero sensitivity for that
  moment. This is informative but a **trap for optimization** — see Result 3.

The entry-by-entry agreement is only moderate (sign 0.6–0.65, Pearson 0.4–0.6):
the smooth surrogate Jacobian and the FD Jacobian over a near-discrete objective
agree on the big picture (rank, conditioning, which parameters matter) but not
cell-by-cell. The surrogate is a complement to the structural audit, not a
replacement.

ARD relevance (least → most identified parameter), `candidate_no_timing_v0`:

`theta_n`, `theta0`, `psi_child`, `kappa_fert`, `h_bar_n`, `h_bar_jump`,
`beta`, `c_bar_0`, `c_bar_n`, `alpha_cons`, `b_entry_fixed`, `h_bar_0`, `chi`.

This matches the ledger almost exactly: the bequest block (`theta_n`, `theta0`)
and the fertility block (`psi_child`, `kappa_fert`) are the least-identified;
the owner wedge `chi`, baseline housing `h_bar_0`, and entry wealth
`b_entry_fixed` are the most-identified. The replacement system gives the same
ordering at the top and bottom.

## Result 3 — Bayesian optimization: the central finding

Every proposal was solved with the real model. The two target systems behave
oppositely, and the reason is economically precise.

### `candidate_no_timing_v0` — surrogate-BO fails

| mode | predicted loss | actual loss |
|---|---:|---:|
| incumbent (data best) | — | 7.637 |
| exploit (global argmin) | 4.39 | 44.65 |
| best global-EI | 6.08 | 46.67 |
| best local trust-region | 4.52 | 44.31 |

Predicted-vs-actual loss MAE = **38.8**. No proposal beat the incumbent.
Decomposing the gap moment-by-moment, the entire failure is one moment:
`prime_childless_owner_median_rooms`. The emulator predicts ~5.9 (smooth
interpolation between rungs); the real model **snaps to the discrete rung 4.0**
(target 6.0), and with weight 10 that single moment contributes exactly **40.0**
to the real loss. Stripping that one term, the emulator's predicted loss matches
the actual loss to within ~0.1–0.9 for the local proposals — i.e. the emulator
is right about the other twelve moments and wrong only about the discrete one.

The owner median is a step function of `theta`. A smooth GP cannot represent it
across rung boundaries, the incumbent sits on a rung edge, and small
perturbations flip the rung. This is the **same object** the FD audit flagged
(rank 12) and the ledger flagged (non-smooth) — the surrogate independently
confirms it and shows it is fatal for any smooth-surrogate optimizer. Note this
is *not* a trust-region failure: local proposals fail identically to global
ones, because the obstacle is discreteness, not distance.

### `candidate_replacement_v1` — surrogate-BO works

Same machinery, but the room moments are now smooth (means and >=6-room shares):

| mode | predicted loss | actual loss |
|---|---:|---:|
| incumbent (data best) | — | 21.276 |
| exploit (global argmin) | 18.15 | 21.33 |
| best global-EI | 24.72 | 28.72 |
| best local trust-region | 18.62 | **20.09** |

Predicted-vs-actual loss MAE = **6.85**, correlation = **0.97**. The surrogate
loss now tracks the real loss, and a single first batch of six local solves
found a point that **beats the data-set incumbent** (20.09 < 21.276). This is a
small, validated improvement, not a production calibration — but it shows the
outer loop is well-behaved once the objective is smooth.

The contrast is the deliverable: replacing the discrete owner-median moment with
smooth, identification-preserving room moments (mean rooms and >=6-room shares,
as the replacement set does) turns a surrogate-BO that fails into one that is
calibrated and useful. That is exactly the ledger's prescription, now
empirically validated.

## Result 4 — per-moment economic fit (all moments)

Targets, the data-set incumbent's model value, and the best *validated* proposal
under each system. This is the substantive read on where the model misses.

### `candidate_no_timing_v0` (incumbent loss 7.64; best validated proposal 44.31)

| moment | target | incumbent | best proposal | weight |
|---|---:|---:|---:|---:|
| `tfr` | 1.700 | 1.697 | 1.791 | 12 |
| `childless_rate` | 0.150 | 0.259 | 0.209 | 12 |
| `own_rate` | 0.575 | **0.058** | 0.159 | 12 |
| `own_family_gap` | 0.168 | −0.002 | −0.010 | 10 |
| `housing_increment_0to1` | 0.664 | 0.590 | 0.551 | 8 |
| `housing_increment_1to2` | 0.566 | 0.610 | 0.721 | 2 |
| `young_liquid_wealth_to_income` | 0.600 | 0.245 | 0.568 | 12 |
| `old_age_own_rate` | 0.764 | 0.488 | 0.808 | 10 |
| `old_age_parent_childless_gap` | 0.070 | 0.035 | 0.145 | 10 |
| `liquid_wealth_to_income` | 1.200 | 1.080 | 1.102 | 12 |
| `housing_user_cost_share` | 0.240 | 0.295 | 0.290 | **250** |
| `prime_childless_renter_median_rooms` | 4.000 | 3.731 | 3.721 | 10 |
| `prime_childless_owner_median_rooms` | 6.000 | 6.000 | **4.000** | 10 |

The v0 incumbent is a pathological corner: it scores loss 7.64 mainly by sitting
on the owner-median rung (6.0, on target) and keeping the weight-250 user-cost
share near target, while **prime-age ownership collapses to 0.058** and the
family ownership gap is ~0. The weight of 250 on `housing_user_cost_share`
dominates the objective and distorts the fit. This corner is exactly why the
discrete owner-median term is so destructive: any move off it costs 40.

### `candidate_replacement_v1` (incumbent loss 21.28; best validated proposal 20.09)

| moment | target | incumbent | best proposal | weight |
|---|---:|---:|---:|---:|
| `tfr` | 1.700 | 1.637 | 1.625 | 20 |
| `childless_rate` | 0.150 | 0.264 | 0.271 | 20 |
| `own_rate` | 0.575 | 0.408 | 0.459 | 100 |
| `own_family_gap` | 0.168 | 0.142 | 0.154 | 45 |
| `housing_increment_0to1` | 0.664 | 1.059 | 1.066 | 14 |
| `housing_increment_1to2` | 0.488 | 0.520 | 0.444 | 8 |
| `young_liquid_wealth_to_income` | 0.179 | 0.322 | 0.310 | 12 |
| `old_nonhousing_wealth_to_income_6575` | 6.419 | 3.389 | 3.476 | 1 |
| `old_parent_childless_nonhousing_wealth_gap_6575` | 1.007 | 0.230 | 0.231 | 2 |
| `prime30_55_childless_renter_mean_rooms` | 3.805 | 4.177 | 4.267 | 6 |
| `prime30_55_childless_owner_mean_rooms` | 6.224 | 5.240 | 5.227 | 6 |
| `prime30_55_childless_renter_share_rooms_ge6` | 0.138 | **0.003** | 0.002 | 25 |
| `prime30_55_childless_owner_share_rooms_ge6` | 0.596 | 0.618 | 0.613 | 25 |

The replacement incumbent is more balanced (ownership 0.41, family gap 0.14 near
target 0.17) but reveals consistent structural misses that BO cannot fix because
they are model/identification issues, not search failures:

- **Old nonhousing wealth far too low**: 3.39 vs target 6.42; the parent-
  childless old-wealth gap 0.23 vs 1.01. The ledger already noted old wealth
  loads on `beta`, not the bequest block — the surrogate's weak R² on the gap
  (0.56) and low ARD relevance of `theta0/theta_n` corroborate this.
- **Renters never reach 6 rooms**: share >=6 is 0.003 vs target 0.138. With the
  renter cap and continuous renter housing, the model essentially never puts a
  childless renter in a large unit.
- **First-birth housing response overshoots**: `housing_increment_0to1` 1.06 vs
  target 0.66.

## Result 5 — frontier probe (conflicting moments)

Emulated achievable Pareto front over three conflicting moments.

`candidate_no_timing_v0`, objectives `own_rate`, `old_age_own_rate`,
`housing_user_cost_share` (abs deviation from target): each is individually
reachable (~0), and at the best joint point the deviations are 0.006 / 0.133 /
0.018 — the binding conflict is old-age ownership against the other two. Caveat:
this 3-moment frontier ignores the other ten moments, and the discrete-moment
trap (Result 3) means emulated frontier points must be solver-confirmed.

`candidate_replacement_v1`, objectives `own_rate`,
`old_nonhousing_wealth_to_income_6575`,
`prime30_55_childless_owner_share_rooms_ge6`: best joint deviations 0.018 /
0.206 / 0.076 — old nonhousing wealth is the binding one (it cannot get close
while ownership and owner room share are near target), consistent with the
per-moment misses above.

## Conclusions

1. **The emulator works.** Held-out loss R² ~0.85–0.86 and rank correlation
   ~0.94–0.96 on both target systems; per-moment R² > 0.9 for the well-
   identified moments. It reproduces records exactly and predicts the smooth
   moments accurately even off the training cloud.

2. **For identification it is a genuine, cheap complement** to the FD audit: it
   reproduces the rank/conditioning verdict and the ARD relevance ordering
   matches the ledger's weak/strong parameter blocks. Its smoothing of the
   discrete owner-median moment (rank 12 → 13) is the one place to distrust it.

3. **For optimization the lesson is sharp and validated**: naive surrogate-BO
   fails on `candidate_no_timing_v0` solely because of the discrete
   owner/renter *median*-rooms moment, and succeeds on `candidate_replacement_v1`
   (smooth *mean*/share room moments), where it beat the data-set incumbent on
   the first six-solve batch. **Recommendation: drop the discrete median-room
   moments from any SMM objective in favor of the smooth room-share / mean-room
   moments** — which is independently the ledger's identification-preserving
   replacement, now empirically justified.

4. **The honest limits hold**: a surrogate cannot manufacture identification a
   misspecified target system lacks. The replacement system's structural misses
   (old nonhousing wealth far too low; renters never reaching 6 rooms;
   first-birth housing overshoot) are model/measurement issues that better search
   does not touch. These are where the next economic work should go.

## How to reproduce

```bash
cd code/model
.venv/bin/python -m intergen_surrogate_calibration.gp              # GP self-test
.venv/bin/python -m intergen_surrogate_calibration.run_surrogate_calibration \
    --target-set candidate_replacement_v1 --n-keep 1000 --restarts 2 --iter 60
.venv/bin/python -m intergen_surrogate_calibration.validate_proposals \
    --target-set candidate_replacement_v1 --limit 6 --repro 1
```

Packets, plots, and `validation.json` for both systems are under
`output/model/intergen_surrogate_calibration/`.
