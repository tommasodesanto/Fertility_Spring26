# Parallel nearby calibration starts

## Specification and verified source facts

The author requested more parallel computation for tonight's preliminary
calibration. This packet adds nine independent searches to the three already
running in `entry_calibration_round1_v1`; those existing jobs are untouched.
There are four starts per specification/grid, twelve tasks in total.

| Specification | Grid | New tasks | Unsecured limit |
|---|---|---:|---:|
| Option 1: five empirical wealth/income-bin means times current annual entrant earnings | 120 assets × 9 incomes | 3 | 0.25 |
| Option 3: negative bin means censored to zero, positive means rescaled to preserve the original mean | 120 × 9 | 3 | 0 |
| Option 3: same nonnegative wealth rule | 160 × 15 | 3 | 0 |

The unsecured limit is in economy-wide mean annual gross working-earnings units.
The positive-node multiplier is 0.3632385158888715; this transforms five bin
means, not raw survey records. The common annual real rate remains 2% for saving
and borrowing. Corrected homeowner collateral/full-sale repayment, mortality
repayment, earnings, taxes, pensions, targets and weights are unchanged from the
pilots and round 1. Relative to September 28 block0506, the experimental economic
changes remain corrected credit accounting, the entrant wealth/borrowing rules,
and birth-renewal price/population-housing closure. The grid change is numerical.
Option 3 is provisionally preferred pending the author's discussion with Corina;
these results are experimental and do not adopt a production reference.

House price clears birth renewal and population clears absolute housing supply.
H0 and psi_child stay fixed; nine free parameters face ten scored targets.
All 14 target and 31 parameter/bound rows plus the same 17 PNG diagnostics must
be reported. The actual authenticated 160×15 bundle supplies the fine grid;
there is no interpolation of coarse inputs into a nominally fine grid.

## Starting vectors and evidence scope

Each specification gets three nearby starting vectors around its pilot's
verified selected point. In parameter order
`beta_annual, chi, first_birth_fixed_cost, kappa_fert,
kappa_fert_continuation, theta0, delta_alpha_jump,
child_benefit_curvature, tenure_choice_kappa`, define

\[
h_j=\max\{10^{-6},\min[0.02\max(|p_j|,0.01),0.005(U_j-L_j)]\},
\qquad p_j^{(s)}=\operatorname{clip}(p_j+3h_jd_j^{(s)},L_j,U_j).
\]

Directions are `s1=[1,-1,-1,1,-1,1,-1,1,-1]`, `s2=-s1`, and
`s3=[-1,-1,1,1,1,-1,1,-1,1]`. These are modest local perturbations,
mostly about 6%, with annual discount-factor changes of 0.00075. They add
limited starting-point diversity, not a global parameter search. Both option 3
grids have exactly the same starting vectors for each suffix.

`plan.json` records parent coordinates, finite-difference scales, requested and
effective offsets, bounds, parent plan identity and verified pilot receipt hashes.
**The two recorded parent repeats verify the pilot parent only. The perturbed
starts have zero native verification repeats until their own jobs run.**
`parent_seed_parameter_table` is the parent's full table, not a table of the
perturbed proposal. Each new job must pass a fresh native baseline and its
independent repeat before search; failure stops the job with no fallback.

`prepare_multistart.py` deterministically creates the plan and pins. The three
base lane configurations are retained in that plan solely for unchanged test
compatibility. `dispatch_lanes` contains exactly nine new lane names; the base
three must not be submitted again.

## Clock, search and verification

All jobs share **September 30, 10:38:17 p.m. New York** as their absolute hard
end (`2026-10-01T02:38:17Z`, epoch 1790822297). These new jobs receive only the
remaining time; they do not receive another four hours from their later start.
The launcher passes the common epoch to the unchanged runner. One CPU, one
numerical thread and 24 GiB per task. Native evaluation caps remain 80 GE for
each coarse lane and 36 for each fine lane, 20 lifecycle evaluations per GE,
300 seconds per lifecycle/observer stage, and eight/four search rounds.

`runner.py`, `inputs.py`, `phase_b_pilot.py` and the original 13-test suite are
byte-identical to round 1. The iterative Jacobian/ridge/bounded Gauss–Newton
algorithm, accounting gates, exact-repeat comparisons, target fingerprint and
budget reserves are unchanged. Each selected result still needs two independent
full-GE repeats and the unchanged native diagnostics. Incomplete derivatives,
rank deficiency, native failures or insufficient verification must remain
visible. No restarts, budget extensions, economics changes or relaxed gates.

The local suite runs the 13 original tests over all twelve configurations and
four added tests for unique offsets, bound clipping, parent-proof labels,
matching coarse/fine starts, source identity and the common deadline. These
are zero-model-solve tests. Every submitted extra lane must also pass the exact
cluster launcher preflight before its native solve. Mock search success does
not certify native feasibility, fit or convergence.

```sh
cd output/model/fixed_reference_economics_20260928/entry_calibration_multistart_v1
../../../../code/model/.venv/bin/python -m unittest -v test_round1 test_multistart
```

Preparation and source freeze: see `preparation_verification.json`.
Deployment is owned by `deployment/`; launch/preflight receipts will be linked
here once the lead submits. The monitoring/collection contract remains compact
14/31 tables, source and repeat receipts, and standard PNGs, without large model
arrays or unsolicited PDFs. Keep failures and verified results separate when
selecting across starts.

## Launch

Additional array **18900753** submitted at18:54:25NewYork, alongside original
18899847. All nine preparatory jobs18900620 passed their exact CLI loops
(80mockGE in coarse lanes,36in fine lanes,two repeats,zero native lifecycle
solves). Shared17tests passed. Both arrays share22:38:17NewYork hard end.
See `launch.json`, `lead_verification.json` and `deployment/preflight_receipts/`.
The monitor covers all12chains; no automatic restart, extension or adoption.

All nine additional tasks started18:56:27NewYork. Together with the original
three,12single-core searches run across cs605/cs608/cs612/cs614.
