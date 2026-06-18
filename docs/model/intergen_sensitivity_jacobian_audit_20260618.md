# One-Market Intergen Sensitivity Jacobian Audit

Date: 2026-06-18

## Run

This audit implements the finite-difference check proposed in
`docs/model/intergen_one_market_identification_ledger_20260618.md`.

- Script: `code/model/tools/audit_intergen_sensitivity_jacobian.py`
- Launcher: `code/cluster/submit_intergen_sensitivity_jacobian.sh`
- Slurm job: `11077696`
- Status: completed with exit code `0:0`
- Runtime: 7 minutes 4 seconds
- Solves: 54 finite-difference cases, all `ok`
- Output bundle: `output/model/intergen_sensitivity_jacobian_20260618/`

The audit perturbs the 13 internal parameters in
`candidate_no_timing_v0`:

\[
\Theta =
\{\beta,\alpha,b_0,\bar c_0,\bar c_n,\bar h_0,\bar h_{\mathrm{jump}},
\bar h_n,\psi_{\mathrm{child}},\kappa_n,\chi,\theta_0,\theta_n\}.
\]

It uses a 1 percent finite-difference step, `J=16`, `Nb=60`, five income
states, `n_house=5`, and `max_iter_eq=3`. The scaled Jacobian is

\[
\widetilde J_{mp}
=
\frac{\partial m_m}{\partial \theta_p}
\frac{\max\{|\theta_p|,1\}}{\max\{|m_m|,|\hat m_m|,1\}},
\]

except that \(\beta\) uses its own level as the parameter scale. This is a
conditioning diagnostic, not an elasticity.

## Baseline Points

Two economically different frontier points were audited:

| Point | Rank loss under `candidate_no_timing_v0` | Residual | Rank | Condition |
|---|---:|---:|---:|---:|
| `core_feasibility_v1` | 28.416 | \(6.36\times 10^{-6}\) | 13 | \(2.69\times 10^4\) |
| `roomcost_test_v1` | 58.034 | \(4.69\times 10^{-4}\) | 12 | infinite |

The core point is technically full rank, but the smallest singular value is
only \(6.7\times 10^{-4}\), so the system is nearly singular. The room-cost
point is locally rank deficient.

## Main Findings

1. The target system is not underidentified by raw count, but it is weakly
   identified in economically important blocks.

   The 13-by-13 target Jacobian has full rank at the core point, but it is
   badly conditioned. This means the count is formally adequate, while local
   identification is fragile.

2. The owner median room moment is not a reliable local identifying moment.

   At the room-cost point, `prime_childless_owner_median_rooms` has zero local
   derivative with respect to every internal parameter, so it contributes no
   local rank. At the core point it behaves like a discrete threshold: the
   median is 6 at the baseline, but small perturbations can jump it to 4. This
   is a non-smooth statistic, not a stable local moment.

   Identification-preserving replacement: use owner room-bin shares, owner
   rung shares, or the owner-renter room gap conditional on childless prime-age
   households.

3. The bequest parameters are not cleanly identified by the current old-age
   ownership moments.

   At the core point, the column norms for `theta0` and `theta_n` are tiny
   relative to the other parameters. At the room-cost point they move housing
   increments and ownership margins more than the old-age ownership moments.
   Meanwhile `old_age_own_rate` responds mostly to `chi`, `alpha_cons`, and
   `beta`, not to the bequest block.

   Identification-preserving replacement: add old liquid wealth, old
   parent-childless wealth gaps, downsizing/retention, or old owner room-bin
   moments before treating `theta0` and `theta_n` as internally calibrated.

4. The housing cost-share moment mainly identifies a combined
   \((\alpha,\chi)\) housing-demand/tenure block.

   In both audited points, `housing_user_cost_share` responds most to
   `alpha_cons` and `chi`. It is not enough by itself to separate the Stone-
   Geary housing share from the owner/renter service wedge.

   Identification-preserving replacement or supplement: tenure-specific
   housing cost shares, renter rent-to-income, imputed owner user-cost-to-income,
   and room distributions by tenure.

5. The fertility block is entangled with housing and wealth parameters.

   `tfr` and `childless_rate` are moved by `c_bar_n`, `alpha_cons`, `c_bar_0`,
   `h_bar_0`, and sometimes `beta`, not only by `psi_child` and
   `kappa_fert`. The strongest collinearity at the room-cost point is
   `psi_child` versus `kappa_fert` with cosine \(-0.984\), and both are also
   highly collinear with `theta_n`.

   Identification-preserving supplement: use coherent parity shares or another
   fertility-composition moment if `psi_child` and `kappa_fert` remain internal.

## Strongest Local Collinearities

At `core_feasibility_v1`:

| Parameter pair | Cosine |
|---|---:|
| `h_bar_jump`, `h_bar_n` | 0.965 |
| `h_bar_0`, `h_bar_jump` | 0.895 |
| `h_bar_0`, `h_bar_n` | 0.895 |
| `b_entry_fixed`, `h_bar_n` | -0.869 |
| `c_bar_n`, `psi_child` | -0.821 |

At `roomcost_test_v1`:

| Parameter pair | Cosine |
|---|---:|
| `psi_child`, `kappa_fert` | -0.984 |
| `psi_child`, `theta_n` | 0.975 |
| `kappa_fert`, `theta_n` | -0.946 |
| `kappa_fert`, `theta0` | -0.900 |
| `beta`, `kappa_fert` | 0.892 |

## Interpretation For Calibration Design

Do not drop weak moments without replacement. The audit says the formal
13-parameter / 13-moment count is not the whole problem. The replacement logic
should be:

| Weak current object | Problem | Replacement preserving identification |
|---|---|---|
| `prime_childless_owner_median_rooms` | discrete, locally flat or threshold-jumpy | owner room-bin/rung shares or owner-renter room gap |
| `old_age_own_rate` | mostly moved by tenure and saving parameters | old wealth, old retention/downsizing, old room-bin moments |
| `old_age_parent_childless_gap` | not enough to isolate `theta_n` | parent-childless old wealth or housing-retention gap |
| `housing_user_cost_share` | entangled \((\alpha,\chi)\) object | tenure-specific cost shares and room distributions |
| `tfr` / `childless_rate` alone | weakly separate child utility, child costs, and dispersion | parity-composition moment or externally fixed fertility-shape parameter |

The next calibration step should therefore be a target-measurement audit and
replacement-moment design, not another blind global search on the current
13 hard moments.
