# Spec Audit B: Bequest utility normalization and terminal paths

Auditor: bequest specialist subagent. Date 2026-07-11. Object: dirty working tree at HEAD fe092ca
(the tree the July 10 overnight snapshot was taken from).

Intended display (paper `latex/intergenerational_housing_fertility_paper_draft.tex:504-513`):

    B(W^B,n) = theta0*(1+theta_n*n) * [ (theta1+max{W^B,0})^{1-sigma} - theta1^{1-sigma} ] / (1-sigma)
    W^B = b + P*h  (gross of mortgage; "no forced sale, no transaction cost, and no tax at death")
    sigma=1: log(theta1+max{W^B,0}) - log(theta1);  B(0,n)=0.

Implementation: `code/model/intergen_housing_fertility/solver.py:5604-5620` (`bequest_utility_vec`).
NOTE: the July 9 second-opinion citation "solver.py:5067-5076" is stale in this tree; those lines
are now moment-accounting code. The function moved to 5604-5620.

## Verdict

The normalized bequest specification used by the July 10 overnight is implemented correctly on
every terminal path of the active evaluation chain, and the seven-test unit file passes. The one
material problem is a staleness trap: `normalize_bequest_utility` defaults to `False`
(`parameters.py:66`) and none of the intergen diagnostic/policy tools nor most CLI subcommands set
it, so any post-hoc re-solve of a combined-spec record through those entry points silently reverts
to the old unnormalized bequest (level shift theta0*(1+theta_n*n)/theta1 = 13.2 utils childless,
33.4 utils at n=2, at the live candidate theta).

## 1. Formula, term by term (VERIFIED)

`solver.py:5604-5620`:

    b_gross = np.maximum(b, 0.0)                                   # max{W^B,0}
    taxable = np.maximum(b_gross - estate_tax_exemption, 0.0)
    b = np.maximum(b_gross - estate_tax_rate * taxable, 0.0)       # tax BEFORE utility; defaults 0
    scale = P.theta0 * max(1 + P.theta_n * nk, 0)
    if abs(P.sigma - 1) < 1e-6:  utility = log(theta1+b) [- log(theta1) if normalize]
    utility = (theta1+b)**(1-sigma)/(1-sigma) [- theta1**(1-sigma)/(1-sigma) if normalize]
    return scale * utility

Numeric check `scripts/check_bequest_formula.py` (all PASS):
- exact match to the paper display at sigma=2, theta1=0.01, for n in {0,1,2}, W in [-2, 20];
- B(0,n)=0 and B(W<0,n)=0 exactly under the flag;
- weakly increasing in W (flat only in the clip region W<0), strictly increasing for W>=0;
- marginal ratio n=2 vs n=0 equals 1+2*theta_n exactly;
- old-vs-new gap equals theta0*(1+theta_n*n)/theta1 exactly (sigma=2).

sigma: fixed at 2.0 (`parameters.py:79`), NOT searched (`production_profile.py:39-53` has only
theta0 in [0,2], theta_n in [0,1.5] among bequest params). The sigma=1 branch is an unreachable
guard in production; its formula is nonetheless correct (unit test passes). theta1=0.01 fixed
(paper line 952), not searched, so theta1^{1-sigma} is finite and theta1+W >= 0.01 > 0 after the
clip — no negative-base or singular cases. The `max(1+theta_n*nk, 0)` floor is inactive
(theta_n bound is [0,1.5]).

## 2. Gross estate and estate tax (VERIFIED)

Both Bellman variants construct the terminal array identically:
- markov-income variant `solver.py:2063-2070`: `hv = p_hat[i] * P.H_own[ten - 1] if ten > 0 else 0.0`
  and `Vbq[:, ten, i, nn, cs] = bequest_utility_vec(b_grid + hv, nk, P)`;
- core variant `solver.py:2431-2438`: identical lines.

So W^B = b + p*H at full market value: gross of mortgage debt (b_grid is negative for levered
owners, floor -phi*p*H per `solver.py:2059-2061`), no (1-psi) sale haircut (contrast `heq` at
2053/2418, used only for living sales), matching the paper's "no forced sale, no transaction
cost". Within the borrowing limit W^B >= (1-phi)p*H > 0 for owners, so the clip binds only on
infeasible grid nodes. Renters: ten=0 gives hv=0, estate = max(b,0), correct. Estate tax:
machinery exists but `estate_tax_rate = 0.0`, `estate_tax_exemption = 0.0` defaults
(`parameters.py:67-68`); no estate_tax override anywhere in `tmp/overnight_combined_20260710/`
(rg: zero hits), consistent with the paper's "no tax at death". Ordering (tax -> utility) is
correct and unit-tested.

## 3. Family-size scaling: encoded index (VERIFIED, matches paper intent)

`nk = get_completed_fertility(nn, cs, P)` (`solver.py:2069/2437`, def at 5593-5601): cs=0 -> 0,
cs=K+1 -> 1, cs=K+2 -> 2, else nn. So bequest scaling uses the ENCODED index n in {0,1,2}, not the
doubled CPS count. The paper explicitly intends this: "Preferences, child needs, bequests, and
demographic flows are evaluated in index units" with completed fertility reported as 2E[n]
(paper lines 600-607). No mismatch. Caveat for interpretation only: theta_n multiplies index
units, so per-child statements must halve it.

## 4. Terminal paths (VERIFIED)

- Warm glow enters as the j=J-1 continuation: `Vnr = Vbq` at `solver.py:2095` (markov) and
  `solver.py:2471` (core); it is discounted by beta inside eval_renter/eval_owner like any
  continuation (timing convention, consistent across tenures).
- No death before J: no survival/mortality code anywhere in solver.py (rg: zero hits).
- Child aging at the terminal boundary: `Vc = apply_child_aging(Vbq, ...)` is applied, but
  completed fertility nk is invariant under both the deterministic aging map
  (`solver.py:5420-5432`) and the stochastic `Pi_child` (matured mass routed nn=1 -> K+1,
  nn>=2 -> K+2, `parameters.py:489-494`), so the bequest scale is unaffected. Check 6 PASS.
- Off-grid: Vbq is exact at grid nodes and weakly monotone; linear (and pchip) interpolation in
  the golden-section b' search preserves monotonicity. theta1+W < 0 cannot occur (clip).

## 5. Overnight evaluation chain uses the normalized form everywhere (VERIFIED)

- wave 1 `tmp/overnight_combined_20260710/run_overnight_combined.py:155` kwarg
  `normalize_bequest_utility=True` AND `:163` in `fixed_model_overrides`;
- wave 2 `run_wave2_lateral.py:74,141,149` same double plumbing;
- `local_panel.py:578` `profile_extra_overrides["normalize_bequest_utility"] = bool(...)` then
  `:579` `.update(fixed_model_overrides or {})` — the True in fixed overrides wins even if the
  kwarg were left at default;
- `run_local_panel_case` merge `local_panel.py:848-853`: extra_overrides applied over
  base_overrides; theta (searched params only) applied last cannot clobber the flag.
The unnormalized branch is reachable ONLY when the flag is False/absent — as intended.

## 6. Tests

`code/model/intergen_housing_fertility/tests/test_bequest_normalization.py`: 7/7 PASS via
`python -m unittest` (pytest is NOT installed in code/model/.venv — ran unittest instead).
Tests cover B(0,n)=0, strict monotonicity, marginal family-size ratio, sigma=1 branch, CRRA
branch, legacy branch, and estate-tax-before-utility ordering.

`tools/compare_intergen_bequest_normalization.py` is the fixed-theta A/B + bounded-polish tool;
it sets the flag explicitly on both arms (`compare_intergen_combined_specification.py:233,256`
pattern; A/B seed at output/model/intergen_bequest_normalization_ab_20260710). Not re-run
(battery constraint); code inspected only.

## FINDINGS

### F1 (MAJOR, staleness trap): default-False flag silently reverts re-solves to the old bequest

- `parameters.py:66` `P.normalize_bequest_utility = False`; `solver.py:5611`
  `getattr(P, "normalize_bequest_utility", False)`.
- Intergen tools that re-solve records and NEVER set the flag (nor q/delta4/eta_supply/H0/
  Rouwenhorst of the combined fixed spec):
  `tools/audit_intergen_sensitivity_jacobian.py:316-320`,
  `tools/audit_intergen_final_best_pathologies.py:146-148`,
  `tools/audit_intergen_solver_accuracy.py:206-208`,
  `tools/audit_intergen_parent_credit_margin.py`, `tools/audit_intergen_mechanism_grid.py`,
  `tools/compare_intergen_income_processes.py:203-208`,
  `tools/run_intergen_code_repair_comparison.py`, `tools/run_intergen_policy_poc.py:193-203`.
  All build `overrides = {**base_overrides(...), **income_process_overrides(...), **theta}`;
  `base_overrides` (`calibration.py:997`) has no flag.
- CLI: only `local-polish` exposes `--normalize-bequest-utility` (store_true, default False,
  `cli.py:158,266`); `smoke/solve/diagnostics/local-panel/global-de-panel` cannot set it at all.
- Economic size at the live candidate (theta0=0.13183, theta_n=0.76792, theta1=0.01): reverting
  shifts terminal utility by -theta0*(1+theta_n*n)/theta1 = -13.18 (n=0), -23.30 (n=1),
  -33.43 (n=2) utils — a family-size-DEPENDENT level drop that changes fertility and saving
  policies, so replayed moments will not match the overnight record and any old-vs-new bequest
  comparison run through these tools is contaminated.
- Decisive test: solve the combined-spec theta once through
  `audit_intergen_sensitivity_jacobian.solve_theta` and once through the overnight overrides at
  Nb<=40 and diff moments; or simply assert
  `getattr(P, "normalize_bequest_utility")` after each tool's override merge.

### F2 (MINOR): per-record artifact under-specifies the model

`run_local_panel_case` return dict (`local_panel.py:879-894`) stores theta/moments/loss but not
the resolved overrides or the flag; only the run-level meta (`local_panel.py:621-624`) records
`normalize_bequest_utility` and `fixed_model_overrides`. `best.json` alone cannot be replayed
correctly — F1 tools take exactly this record as input.

### F3 (MINOR): stale line provenance

The July 9 second-opinion citation solver.py:5067-5076 for the old unnormalized form no longer
points at bequest code in this tree (now moment accounting); the function is at 5604-5620.
Line-based provenance from July 9 must be remapped before being quoted in the audit report.

### F4 (SUSPECTED, identification economics): bequest motive saturates at tiny estates

With theta1=0.01 and sigma=2, normalized B is bounded above by theta0*(1+theta_n*n)/theta1 and
reaches 90% of that bound at W^B=0.09 and 98% at W^B=0.5 (model wealth units are order 1). The
motive is effectively extensive-margin ("leave any positive estate"): marginal utility at W=0 is
10^4*scale but only ~scale at W=1. theta_n is then identified mostly by a discrete level jump in
parents' terminal value rather than by differential saving slopes, making (theta0, theta_n)
near-substitutable with any parenthood terminal-value shifter. Not a code bug; hand to the
identification adjudication (task 9). Decisive test: local Jacobian of target moments w.r.t.
theta1 in {0.01, 0.1, 0.5} at fixed theta0*(...)/theta1.

## Blocked

Cluster-side confirmation that the snapshot's solver.py:5604-5620 is byte-identical to this tree
(SSH unavailable) — covered by the provenance task's hash manifest if it captured solver.py.
