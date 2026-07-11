# Finance and Housing-Supply Audit (spec audit C)

Audit of the July 10 overnight combined-spec evaluation chain, dirty tree at
HEAD fe092ca. All numeric checks in
`output/model/full_audit_20260711/scripts/finance_supply_checks.py`, output in
`output/model/full_audit_20260711/spec_audit/finance_supply_checks.json`.
Tiny verification solve: Nb=32, J=17, full production+combined overrides,
current_bound_best theta (residual 7.9e-5, strict).

## 1. Defaults of q and delta (what runs when the fixed spec is NOT passed)

`parameters.py:80-81`:

```python
P.q = (1.04 ** P.period_years) - 1.0
P.delta = 1.0 - (1.0 - 0.02) ** P.period_years
```

Defaults are annual r=4%, annual delta=2%, compounded to 4-year objects:
q_default=0.1698585600, delta_default=0.0776318400, user_cost_rate_default=
0.2874904 (with tau_H=0.04). The combined fixed spec overrides these to
q=0.0824321600 (annual 2%) and delta=0.0432793094 (annual 1.1%),
user_cost_rate=0.1657114694. So the default gross period return is 1.16986 vs
1.08243 under the overnight spec (q ratio 2.06x, delta ratio 1.79x), and the
default user-cost rate is 1.73x the overnight one.

The "old values" are therefore NOT per-period 0.02*4 or 0.011*4; they are the
pre-revision annual 4%/2% compounded. Nothing in `production_profile.py`
carries q/delta/eta_supply/H0 (`production_profile_overrides()`,
production_profile.py:123-167, has no finance keys), so ANY tool that applies
only the production profile and not the combined `fixed_model_overrides`
silently runs annual 4%/2%. Verified by probe: profile-only overrides give
q=0.169859, delta=0.077632 (`silent_default_drift_probe` in the JSON). The
combined-spec drivers (tmp/overnight_combined_20260710/run_overnight_combined.py:159-164,
run_wave2_lateral.py:69-76, tools/compare_intergen_combined_specification.py:97-98,
tools/run_intergen_combined_recalibration.py, tools/build_intergen_mechanics_packet.py:421-422)
do pass the new values; older audit/plot tools generally do not.

The prompt-sheet value "delta4=0.04328066" is slightly wrong: exact
1-(0.989)^4 = 0.04327930935900004 (asserted to 1e-14 in
tools/compare_intergen_combined_specification.py:102). The code uses the exact
value.

`apply_overrides` recomputes `P.user_cost_rate = P.q + P.delta + P.tau_H` and
`P.R_gross = 1 + P.q` after every override merge (parameters.py:332-333), and
`run_model_cp_dt` recomputes them again (solver.py:905-906), so overridden
q/delta propagate everywhere; there is no stale user-cost path.

## 2. Every use of q

- Liquid return: `Rv = Rg * b + yj` with `Rg = P.R_gross` (solver.py:2007,
  2111; income-type path solver.py:2360, 2479). b runs over the full grid
  including negative (mortgage) positions, so (1+q) applies to assets and debt
  symmetrically: the mortgage rate equals the deposit rate. VERIFIED.
- User cost / rental rate: `r = P.user_cost_rate * p` (solver.py:1252, 1508,
  1589) with user_cost_rate = q+delta+tau_H. Matches paper eq
  (quant_user_cost) q=(r+delta+tau^p)P. VERIFIED.
- PTI cap: `q = max(float(getattr(P, "q", 0.0)), 0.0)`,
  `max_debt_pti = allowed_debt_payment / q` (solver.py:5634, 5647) — inactive
  (use_pti_constraint=False in the production profile).
- No other q use in the intergen package. The deceased's terminal wealth is
  bp + pH without a final-period return (see 3).

## 3. Every use of delta

- Owner per-period cost: `ocst[i, ten] = (P.delta + P.tau_H) * p_hat[i] * hs +
  extra_size_cost` (solver.py:2056, 2421) — maintenance + property tax on
  value, interest cost implicit through R_gross on the (possibly negative) b.
  Matches paper "Owners pay maintenance and property taxes through
  (delta+tau^p)PH_k". VERIFIED.
- Rental company: rent r=(q+delta+tau_H)p covers depreciation. VERIFIED.
- Supply curve: via user_cost in H_s = H0*(uc/r_bar)^eta. VERIFIED.
- There is no physical depreciation of the stock in the KFE: rungs are fixed
  and maintenance implicitly restores the unit. Consistent convention.
- Compounding note (minor): tau_H = 0.01*period_years = 0.04 is a LINEAR 4x
  annual aggregation while q and delta are compounded. Effect on
  user_cost_rate ~0.4% (0.04 vs 1-0.99^4=0.0394). MINOR consistency point.

## 4. Supply side and market clearing

Functional form (solver.py:5484 and identically 1471, 5518, 5568):

```python
supply = P.H0 * (user_cost / P.r_bar) ** P.xi_supply   # user_cost = user_cost_rate * p
```

i.e. H_s = H0 * ((q+delta+tau_H) p / r_bar)^eta, r_bar=0.16 (the reference
rental rate qbar in paper eq quant_supply), eta=xi_supply=1.75 (apply_overrides
copies eta_supply into xi_supply, parameters.py:282-284). eta_supply is a
per-location ARRAY remnant, indexed [i]; I=1 throughout — consistent, no
misindexing.

Units: population is normalized (N_target=1; housing_demand_normalizer,
solver.py:134-137), so demand and supply are in mean rooms per household.
H0 is the supply of rooms per household at user_cost = r_bar. Searching H0 in
[1,10] against a 5.78-room target is dimensionally sensible, BUT see finding
F2: at the candidate equilibrium price the ratio (uc/r_bar)^1.75 = 0.555, so
even H0=10 supplies only 5.55 rooms; H0 needed is 10.41 > bound.

Market clearing: ONE market, one price p per location. Renter rooms (continuous
hR at rent (q+delta+tau_H)p) and owner rooms (discrete H rungs at asset price
p) are pooled: `renter_demand + owner_demand` vs supply (solver.py:4566-4575,
5477-5499). Rental and owner housing are perfect substitutes on the supply
side; the rental rate is tied to p by zero-profit. The damped loop iterates on
the supply-inverse price target
`p_target[i] = P.r_bar[i] * (hd / P.H0[i]) ** (1.0 / P.xi_supply[i]) / P.user_cost_rate`
(solver.py:1025), and the scalar Brent refine drives the residual metric
`abs(excess) / max(abs(supply),1e-12)` (solver.py:1100-1104) below tol_eq=1e-4;
`strict_converged` uses that metric. VERIFIED numerically: the reported
`housing_supply` equals the formula exactly, and the supply identity holds on
all four saved repro records to within the reported market residual (~3e-5).

## 5. Owner menu and renter cap in force

- production_profile.py:17-18: `PRODUCTION_H_OWN = np.array([2.0, 4.0, 6.0,
  8.0, 10.0])`, `PRODUCTION_RENTER_CAP = 6.0`; production_profile_overrides()
  returns both (lines 127-128). base_overrides also sets
  `"H_own": np.linspace(2.0, 10.0, int(n_house))` = [2,4,6,8,10] at n_house=5
  and `"hR_max": 6.0` (calibration.py:1013-1014). Merge order in
  run_local_panel_case (local_panel.py:848-853): {base, **extra(profile+fixed),
  **income, **theta} — no later key touches them.
- Tiny solve confirms on P: H_own=[2,4,6,8,10], hR_max=6.0, and
  max(hR_pol[:,0,...]) = 6.0 exactly (the cap binds for some states); the
  owner slots of hR_pol are all zero.
- The reported candidates were produced through this exact chain (lead
  reproduced current_bound_best bit-consistently), so the menu/cap were in
  force in the overnight results. VERIFIED.

## 6. Down payment, debt floor, (1-phi)

- phi default 0.80*ones(n_parity) (parameters.py:86); nothing in the combined
  chain overrides it (verified P.phi=[0.8,0.8,0.8]).
- Down payment and debt limit (solver.py:2424-2429, comment included):

```python
# phi is the financed share, so the down-payment
# threshold is (1 - phi) * hcost and the borrowing
# limit is -phi * hcost.
dp_arr[i, ten, nn, cs] = (1 - phi_ncs) * hcost[i, ten]
bmo[i, ten, nn, cs] = -phi_ncs * hcost[i, ten]
```

- Buyer feasibility (kernels.py:646-649): renter->owner requires
  `bg_b >= dpn` and post-purchase `bab = bg_b - hc >= bmn` (equivalent
  conditions). Owner->owner: `dpc = dpn - sp` with sp=(1-psi)pH sale proceeds
  (kernels.py:651-655). No code path uses phi itself as the down-payment
  share. VERIFIED.
- Uniform owner debt floor: inside the owner Bellman the savings search is
  bounded below by `lo = bf` where bf = bmo = -phi p H (kernels.py:920-924)
  for ALL active owners, stayers included. Numerically: min bp_pol per rung
  equals -phi*p*H to 1e-15 with no violation. Since the state b next period
  equals the chosen bp, no owner ever begins a period below the floor at a
  stationary price. VERIFIED, but see F4: the PAPER's incumbent carve-out
  `b' >= min{b, -phi P H_k}` (eq quant_incumbent_ltv) is NOT implemented; at
  stationary prices it is vacuous, in counterfactuals with rising p it is not.
- Renter no-borrowing: `lo = 0.0` (kernels.py:803), matching paper b'>=0.
  VERIFIED.
- Terminal/bequest wealth: `Vbq = bequest_utility_vec(b_grid + hv)` with
  `hv = p_hat[i] * P.H_own[ten - 1]` (solver.py:2431-2438, 2063-2070): GROSS
  house value, no (1-psi) liquidation cost, no estate tax (defaults 0), no
  final-period return on b. The paper says W^B is "after liquidation costs and
  taxes" — mismatch, cross-flagged to spec-audit B (known "gross bequests").

## 7. Rooms moment counterpart

`extract_moments` (calibration.py:1090-1092):

```python
"aggregate_mean_occupied_rooms_18_85": float(getattr(sol, "aggregate_housing_demand", np.nan)) / max(total_mass, 1e-12),
```

aggregate_housing_demand = sum of renter rooms (mass-weighted hR_pol at
tenure-0 states) + owner rung rooms (mass times H_own) over ALL ages j=0..16
(ages 18-85, J=17 x 4y from 18) and all states (solver.py:5477-5494),
normalized by N_target; total_mass = sum(g_current) = 1. The distribution used
is g_current = realize_current_cross_section(...) — the POST-decision
cross-section by realized current location/tenure at post-transaction wealth
(solver.py:4014-4029), because the production profile sets
use_postdecision_current_distribution=True (production_profile.py:166).
Timing is stamped "post_housing_choice" (solver.py:4089-4093). So the moment
is post-decision occupied rooms per household, consistent with the July 9
timing repair, and it equals demand (=supply at convergence). Identity
verified exactly in the tiny solve (`rooms_equals_demand_over_mass: true`).

## 8. Owner service flow (kernels.py ~900-940)

Implemented object (kernels.py:915-918, identically solver.py:2183, 2572,
2669):

```python
ht_c = hsv - owner_h_bar_scale * hbc
...
ht_c = owner_service_premium * ht_c        # chi * (H - hbar)
```

with owner_service_premium = chi (solver.py:2012) and owner_h_bar_scale = 1.0
(default, not overridden in the combined chain). The code implements
chi*(H - hbar). NO code comment or docstring claims chi*H - hbar. The PAPER
does: latex/intergenerational_housing_fertility_full_draft.tex:509 has utility
over (h - hbar(n,s)) and :521 defines h = chi_O H_k for owners, which composes
to chi*H - hbar. See finding F1.

## Findings

### F1 (MAJOR, paper/model mismatch): owner housing surplus is chi*(H-hbar) in code, chi*H-hbar in the paper
- Intended (paper): u over (chi_O H_k - hbar(n,s)); full_draft.tex:509+521.
- Implemented: chi*(H_k - hbar(n,s)); kernels.py:915-918, solver.py:2183.
- Consequence: at the candidate (chi=1.15, hbar_0=1.0, jump=1.476,
  h_bar_n=0.984 -> hbar=3.46 for a one-child family), rung H=4 gives paper
  surplus 1.15*4-3.46=1.14 vs code 1.15*(4-3.46)=0.62 — nearly 2x; the code
  version also makes small rungs infeasible whenever H<hbar even though
  chi*H>hbar. Reads directly on the ownership-by-family-size moments and on
  the interpretation of chi as an ownership premium.
- Decisive test: none needed in code (unambiguous reading); fix is a paper
  edit or a one-line kernel change followed by re-calibration.

### F2 (MAJOR, calibration interpretation): H0 sits at its upper search bound; the rooms target is unreachable within [1,10]
- current_bound_best has H0=9.999672 (bound 10.0), p_eq=0.6898; at that price
  (uc/r_bar)^1.75 = 0.555, so max supply = 5.55 rooms vs target 5.78; the H0
  needed at that price is 10.41. rooms gap -0.229, loss contribution
  6*0.229^2 = 0.314 of the 14.78 total. The housing_relaxed candidate is
  interior (H0=8.56) but also under-target (-0.195) because its bound-driven
  chi/hbar block moved demand.
- Consequence: the calibration is a constrained optimum, not an interior SMM
  point; the rooms moment cannot discipline H0 (its designated parameter), so
  effective parameter count and moment count both drop by ~1 and the residual
  rooms misfit spills into other moments. The [1,10] bound is binding
  because r_bar=0.16 (the supply pivot) is inherited from the OLD finance
  regime (old p_init=r_bar/uc=0.556; new 0.966): the same demand now clears at
  uc*p well below r_bar, deflating supply for given H0.
- Decisive test: rerun the candidate with H0 bound widened (e.g. [1,14]) or
  equivalently re-pivot r_bar to the new user_cost_rate; if the optimizer
  moves H0 above 10 and the rooms gap closes, the bound was the binding
  distortion. (Cluster-side; not run here.)

### F3 (MINOR, reproducibility hazard): package defaults keep the old finance parameters
- parameters.py:80-81 defaults are annual 4%/2% (q=0.169859, delta=0.077632,
  uc=0.287490). production_profile_overrides() does not carry q/delta/
  eta_supply/H0, so any diagnostic/plot/audit tool that applies only the
  production profile silently evaluates the OLD finance regime (uc 1.73x the
  overnight one). The combined drivers and the three combined-spec tools do
  pass the overrides. Recommend hoisting the combined finance block into the
  production profile (or a named spec) before any further tool runs.

### F4 (MINOR, paper/code): incumbent-owner no-delever provision not implemented
- Paper eq quant_incumbent_ltv: b' >= min{b, -phi P H_k}. Code enforces the
  uniform floor b' >= -phi p H for stayers too (kernels.py:920-924). Vacuous
  at stationary prices (the state can never start below the floor), NOT
  vacuous in counterfactuals where p rises. Either implement the carve-out
  before price-raising policy experiments or fix the paper sentence.

### F5 (MINOR, convention): tau_H aggregated linearly, q and delta compounded
- tau_H = 0.01*4 = 0.04 vs compounded 1-(0.99)^4 = 0.0394; ~0.4% of the
  user-cost rate. Document or harmonize.

### Cross-flag to audit B: bequest base is gross (b + pH), paper says net of liquidation costs
- solver.py:2434 `hv = p_hat[i] * P.H_own[ten - 1]`; no (1-psi), no final
  return on b. Known project convention ("gross bequests") but inconsistent
  with full_draft.tex:530.

## Verified-correct summary

1. q/delta/eta_supply/H0/H_own/hR_max/phi all land on P through the exact
   overnight merge order (tiny solve, exact equality).
2. (1+q) symmetric on assets and debt; mortgage rate = deposit rate.
3. Down payment (1-phi)pH and borrowing limit -phi pH, with the correct
   owner-to-owner generalization dp - (1-psi)pH_old; no path uses phi as the
   DP share.
4. Uniform owner debt floor binds exactly at -phi p H per rung; renter b'>=0.
5. Rent = (q+delta+tau_H)p; owner flow cost (delta+tau_H)pH; supply
   H0*(uc/r_bar)^1.75; residual metric |excess|/supply with tol 1e-4; supply
   identity holds on all saved repro records to ~the reported residual.
6. Rooms moment = post-decision occupied rooms per household over ages 18-85,
   equals aggregate demand / total mass exactly.
7. The dirty-tree `.reshape(-1)` on the H0 override (parameters.py:299) is
   load-bearing for the scalar searched H0; without it the run would crash —
   the overnight snapshot required the dirty tree, consistent with the lead's
   provenance note.
