"""Fixed-price natural borrowing experiment for authenticated block0506 only.

Economic solvency is derived from worst-case lifetime resources, independently
of utility values. Boolean backward reachability supplies the conservative
grid representation. Its extra restriction is reported, never called an
economic borrowing limit. Execute only on Torch through the pinned harness.
"""
from pathlib import Path
import hashlib
import inspect
import json
import sys
import numpy as np

API_VERSION = 'block0506_credit_adapter_v1'
MODE = 'natural_solvency_boolean_grid_v1'
CHANGES = ['Remove artificial renter, purchaser and incumbent-owner borrowing limits; retain lifetime repayment and net-estate solvency']


def require(ok, text):
    if not ok:
        raise RuntimeError(text)


def write(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def reachable(grid, ok, points):
    """Both positive-weight destination nodes must be feasible."""
    idx = np.clip(np.searchsorted(grid, points, side='right')-1, 0, len(grid)-2)
    w = (points-grid[idx])/(grid[idx+1]-grid[idx])
    inside = (points >= grid[0]) & (points <= grid[-1])
    return inside & ((w >= 1) | ok[idx]) & ((w <= 0) | ok[idx+1])


def construct(model, P, grid, price):
    sd = model.precompute_shared(P, grid)
    require(P.I == 1 and P.preference_spec == 'eqscale', 'Only authenticated one-market eqscale')
    require(np.all(sd.cb_flat == 0) and np.all(sd.hb_flat == 0) and np.all(sd.gb_flat == 0), 'Nonzero floors/transfers')
    require(not model.child_earnings_penalty_active(P) and not model.rental_wedge_active(P), 'Additional spending mechanism')
    require(float(getattr(P, 'owner_size_cost', 0)) == 0 and not model.parent_age_maturation_active(P), 'Additional household mechanism')
    require(not P.use_pti_constraint and not np.any(sd.birth_dp) and not np.any(sd.birth_entry_grant), 'Purchase mechanism differs')
    z, _, pi = model.income_transition_values(P)
    require(np.all(pi > 0), 'This specialization requires every income transition reachable')
    require(P.R_gross > 1 and P.delta + P.tau_H >= 0 and 0 <= P.psi < 1, 'Liquidation dominance assumptions')
    cost = np.r_[0., float(price)*np.asarray(P.H_own)]
    sale = (1-float(P.psi))*cost
    oc = (P.delta+P.tau_H)*cost
    y = np.array([[model.income_at_state(P, 0, j, float(zz)) for zz in z] for j in range(P.J)])
    require(np.all(y > 0), 'Nonpositive earnings/pension')
    survival = np.r_[np.asarray(P.survival_probs), 0.]
    require(survival.shape == (P.J,), 'Survival clock')
    human = np.zeros(P.J)
    minimum = np.zeros((P.J, len(z)))
    floors = np.zeros((P.J, len(cost)))
    pre = np.zeros((P.J, len(z), len(cost), len(grid)), dtype=bool)
    stage = np.zeros_like(pre)
    for j in range(P.J-1, -1, -1):
        candidates = [0.] if survival[j] < 1 else []
        if survival[j] > 0:
            candidates.append(float(minimum[j+1].max()))
        human[j] = max(candidates)
        minimum[j] = (human[j]-y[j])/P.R_gross
        for h in range(len(cost)):
            candidates = [-sale[h]] if survival[j] < 1 else []
            if survival[j] > 0:
                ok = pre[j+1, :, h].all(axis=0)
                require(ok.any() and not np.any(ok[:-1] & ~ok[1:]), 'Nonmonotone/empty Boolean continuation support')
                candidates.append(float(grid[np.flatnonzero(ok)[0]]))
            floors[j,h] = max(candidates)
            require(floors[j,h] >= human[j]-sale[h]-1e-10, 'Discrete floor looser than economic solvency')
        for zz in range(len(z)):
            for new in range(len(cost)):
                stage[j,zz,new] = P.R_gross*grid + y[j,zz] - oc[new] >= floors[j,new]+1e-6
            for old in range(len(cost)):
                for new in range(len(cost)):
                    x = grid if new == old else grid+sale[old]-cost[new]
                    pre[j,zz,old] |= reachable(grid, stage[j,zz,new], x)
    return dict(cost=cost, sale=sale, oc=oc, income=y, survival=survival,
                human=human, minimum=minimum, floors=floors, pre=pre, stage=stage)


def install(*, prepared, reference, P, grid, output, overlay_root, credit_mode):
    require(credit_mode == MODE, 'Unexpected credit mode')
    model = prepared.rt['model']
    require(not getattr(model, '_block0506_credit_installed', False), 'Already installed')
    price = float(reference['solution'].p_eq[0])
    limits = construct(model, P, grid, price)
    P.native_due_stayer_credit = False
    P._credit_limits = limits
    P._credit_price = price
    P._credit_saving_checks = dict(calls=0, minimum_feasible_value=0.)
    original_savings = model._savings_stage
    original_renter = model.full_renter_block_kernel
    original_owner = model.full_owner_block_kernel

    def savings(Vc, P, b_grid, SD, ctx, r_hat, j, z_value, s_next, D_next, renter_floor, stay_floor=False):
        require(not stay_floor, 'Natural rule must apply equally to incumbent and purchaser')
        require(np.array_equal(ctx.current_prices, np.array([P._credit_price])), 'Credit limits built at a different price')
        zz = int(np.flatnonzero(np.asarray(P.z_grid) == z_value)[0])
        q = P._credit_limits
        def renter(*args):
            a = list(args)
            require(len(a) == 34, 'Frozen renter signature differs')
            a[-1] = np.full(SD.nc, q['floors'][j,0])
            return original_renter(*a)
        owner_counter = [0]
        def owner(*args):
            a = list(args)
            owner_counter[0] += 1
            h = owner_counter[0]
            a[12] = np.full(SD.nc, q['floors'][j,h])
            a[21] = 0.; a[22] = 0.
            require(not a[-2], 'DUE owner floor was not removed')
            return original_owner(*a)
        model.full_renter_block_kernel = renter
        model.full_owner_block_kernel = owner
        try:
            result = original_savings(Vc,P,b_grid,SD,ctx,r_hat,j,z_value,s_next,D_next,renter_floor)
        finally:
            model.full_renter_block_kernel = original_renter
            model.full_owner_block_kernel = original_owner
        require(owner_counter[0] == len(q['cost'])-1, 'Unexpected owner branch count')
        values = result[0]
        # This mask comes from resources and backward Boolean support, never V.
        ok = q['stage'][j,zz].T[:,:,None,None,None]
        ok = np.broadcast_to(ok, values.shape)
        require(np.isfinite(values[ok]).all() and np.all(values[ok] > -1e9), 'Finite utility conflicts with native feasibility classifier')
        values[~ok] = -1e10
        P._credit_saving_checks['calls'] += 1
        P._credit_saving_checks['minimum_feasible_value'] = min(P._credit_saving_checks['minimum_feasible_value'],float(values[ok].min()))
        return result

    original_tenure = model._tenure_location_stage
    def tenure(Vd,P,b_grid,SD,ctx,dp_choice,Vd_stay=None,bmo_purchase=None):
        return original_tenure(Vd,P,b_grid,SD,ctx,np.full_like(dp_choice,-np.inf),Vd_stay,
                               np.full_like(ctx.bmo,b_grid[0]))

    # The native logit must not mix feasible and infeasible transaction nodes.
    # Both node indicators are independently established in savings() above.
    kernel = model.tenure_logit_kernel
    function = getattr(kernel, 'py_func', kernel)
    source = inspect.getsource(function)
    require(source.count(', False, transaction_support)') == 5, 'Frozen transaction kernel differs')
    source = source.replace(', False, transaction_support)', ', True, transaction_support)')
    destination = Path(output)/'credit_tenure_kernel.generated.py'
    destination.write_text(source)
    namespace = dict(function.__globals__)
    exec(compile(source,str(destination),'exec'),namespace)
    model.tenure_logit_kernel = namespace['tenure_logit_kernel']
    model._savings_stage = savings
    model._tenure_location_stage = tenure
    model._block0506_credit_installed = True
    rows = []
    for j in range(P.J):
        rows.append(dict(age=float(P.age_start+j*P.period_years),
            economic_renter_limit=float(limits['human'][j]),
            numerical_renter_limit=float(limits['floors'][j,0]),
            maximum_grid_tightening=float(np.max(limits['floors'][j]-(limits['human'][j]-limits['sale'])))))
    write(Path(output)/'borrowing_limits.json', dict(reference_label=reference_label(),rows=rows,
        economic_formula="L_j=max(death floor 0, max reachable future (L_{j+1}-y_{j+1})/R); owner limit L_j-net liquidation value",
        numerical_method='Boolean backward reachability on retained grid; no utility-cutoff derivation',
        grid_refinement_verified=False))
    return dict(status='installed',credit_mode=MODE,economic_changes=CHANGES,
        generated_kernel_sha256=hashlib.sha256(source.encode()).hexdigest(),
        finite_grid_limitation='Conservative support can tighten economic limits; refinement not yet verified')


def reference_label():
    return '2007 stationary reference — block0506, September 28 verified export'


def audit_purchase_accounting(ev,P,sd,grid,model):
    from e5f_solvency_credit_benchmark import audit_purchase_accounting as base
    require(not P.native_due_stayer_credit, 'Audit requires common saving policy for stayers and buyers')
    out = base(ev,P,sd,grid,model)
    out['artificial_credit_limits_enforced'] = False
    out['native_value_cutoff_limitation'] = False
    out['audit_id'] = MODE
    return out


def audit_solvency(packet,prepared,output,*,stationary):
    P, ev, grid = packet['parameters'],packet['evaluation'],packet['b_grid']
    q = P._credit_limits
    g = ev.g_current
    bprime = ev.policy.bp_pol
    badmass = deathmass = tightmass = 0.
    minslack = float('inf')
    for j in range(P.J):
        for zz in range(len(P.z_grid)):
            for h in range(len(q['cost'])):
                mass = g[:,h,0,j,zz]
                bp = bprime[:,h,0,j,zz]
                economic_floor = q['human'][j]-q['sale'][h]
                bad = bp < economic_floor-1e-9
                if q['survival'][j] > 0:
                    ok = q['pre'][j+1,:,h].all(axis=0)
                    bad |= ~reachable(grid,ok,bp)
                badmass += float(mass[bad].sum())
                if q['survival'][j] < 1:
                    deathmass += float(mass[bp < -q['sale'][h]-1e-9].sum())
                positive = mass > 0
                if positive.any():
                    minslack = min(minslack,float((bp-economic_floor)[positive].min()))
                if q['floors'][j,h] > economic_floor+1e-9:
                    tightmass += float(mass[np.abs(bp-q['floors'][j,h]) < 1e-8].sum())
    initial = np.asarray(P.fixed_reference_entry_conditional)
    badentry = sum(float(initial[:,zz][~q['pre'][0,zz,0]].sum()) for zz in range(initial.shape[1]))
    require(float(getattr(P, '_entry_censored_mass', 0.)) == 0., 'Entrants were relocated to a numerical frontier')
    require(badmass <= 2e-10 and deathmass <= 2e-10 and badentry == 0., 'Natural-credit forward/entry solvency failed')
    require(P._credit_saving_checks['calls'] == P.J*len(P.z_grid), 'Not all saving cells independently checked')
    result = dict(status='passed',lifetime_solvency_and_repayment_retained=True,
        occupied_unreachable_continuation_mass=badmass,negative_estate_exposure_mass=deathmass,
        infeasible_fixed_entry_conditional_mass=badentry,
        occupied_mass_at_conservative_grid_floor=tightmass,
        minimum_occupied_slack_above_economic_limit=minslack,
        grid_refinement_verified=False,feasibility_classifier_checks=P._credit_saving_checks,
        finite_grid=dict(nodes=len(grid),minimum=float(grid[0]),maximum=float(grid[-1]),
            conservative_grid_floor_binding_mass=tightmass,refinement_verified=False))
    write(Path(output)/'credit_solvency.json',result)
    return result
