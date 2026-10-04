"""Zero-solve financing census at baseline childless prebirth exposure.

Reads cached arrays only. All sums preserve signed inversion roundoff.
"""
from __future__ import annotations
import csv, importlib.util, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
HELPER = HERE.parent / 'diagnostics/extract_common_states.py'
spec = importlib.util.spec_from_file_location('cached_birth_exposure', HELPER)
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
BASE = helper.DEFAULT / 'phi_080'
RELAXED = helper.DEFAULT / 'phi_100'

def main():
    p, a, pre, _, ages, fec, inversion = helper.load(BASE)
    p1, a1, _, _, _, _, inversion1 = helper.load(RELAXED)
    changes = [k for k in p if not k.startswith('_') and p[k] != p1[k]]
    assert changes == ['phi']
    assert p['native_purchase_income'] and not p['use_pti_constraint']
    assert np.all(np.array(p.get('owner_ltv_multipliers',[1.]))==1)
    assert not p['birth_dp_grant'] and not p['birth_entry_grant'] and not p['parent_dp_waiver']
    assert p['I'] == 1 and np.array_equal(a['b_grid'], a1['b_grid']) and np.array_equal(a['p_eq'], a1['p_eq'])
    extra = ('hR_pol', 'tenure_probs', 'shared.h_bar', 'g', 'g_stay_distribution', 'bp_pol', 'bp_pol_stay')
    with np.load(BASE/'solution_arrays.npz', allow_pickle=False) as z:
        a.update({k:z[k].copy() for k in extra})
    with np.load(RELAXED/'solution_arrays.npz', allow_pickle=False) as z:
        a1.update({k:z[k].copy() for k in extra})
    b, house = a['b_grid'], np.array(p['H_own'])
    price, R = float(a['p_eq'][0]), float(p['R_gross'])
    floor = float(a['shared.h_bar'][1,1])*float(p['owner_h_bar_scale'])
    suitable = house > floor
    assert bool(p['child_room_floor']) and suitable.tolist() == [False, True, True, True, True]
    G = pre[:,:,0,:,:,0,0].copy()  # b, origin tenure, age, income state
    fertile = (np.arange(p['J']) >= p['A_f_start']-1) & (np.arange(p['J']) < p['A_f_end'])
    G[:,:,~fertile,:] = 0
    attempt = a['fert_probs'][:,:,0,:,:,1]
    # Derivative of successful first births w.r.t. success-minus-wait utility.
    W = G*fec[None,None,:,None]**2*attempt*(1-attempt)/p['kappa_fert']
    sale = np.r_[0., (1-p['psi'])*price*house]
    income = np.array(p['income'])[0,:,None]*a['type_values'][None,:]+p['property_tax_lump_sum_transfer']
    cash = R*b[:,None,None,None]+sale[None,:,None,None]
    resources = {'soft': cash+income[None,None,:,:], 'before_current_income': np.broadcast_to(cash,G.shape)}
    groups = {'all':np.ones(G.shape,bool), 'beginning_renter':np.broadcast_to(np.arange(G.shape[1])[None,:,None,None]==0,G.shape), 'beginning_owner':np.broadcast_to(np.arange(G.shape[1])[None,:,None,None]>0,G.shape)}
    # An owner already in product h is a stayer, not a new purchaser of h.
    purchasing = np.arange(6)[None,:,None,None,None] != np.arange(1,6)[None,None,None,None,:]
    rows = []
    def row(kind, group, metric, weight, numerator, denominator, timing='', phi='', product=''):
        rows.append(dict(kind=kind,group=group,timing=timing,phi=phi,product_rooms=product,weight=weight,metric=metric,numerator=float(numerator),denominator=float(denominator),value=float(numerator/denominator) if denominator else None))
    for group, mask in groups.items():
        for weight, w in (('mass',G),('success_gap_susceptibility',W)):
            w = w*mask; den = w.sum()
            row('exposure',group,'total',weight,den,1)
            for timing, res in resources.items():
                afford0 = (res[...,None] >= (1-.8)*price*house)&purchasing
                parent0 = np.any(afford0&suitable,axis=-1)
                for phi in (.8,.95,1.):
                    afford = (res[...,None] >= (1-phi)*price*house)&purchasing
                    parent = np.any(afford&suitable,axis=-1)
                    childless = np.any(afford,axis=-1)
                    for metric, event in (('any_parent_suitable_purchase',parent),('only_unsuitable_small_purchase',childless&~parent),('no_purchase_product',~childless),('newly_parent_suitable_vs_same_timing_phi080',parent&~parent0)):
                        row('affordability',group,metric,weight,np.sum(w*event),den,timing,phi)
                    for k,h in enumerate(house):
                        row('product_affordability',group,'purchaser_screen_pass',weight,np.sum(w*afford[...,k]),den,timing,phi,h)
            current_parent = np.r_[False,suitable][None,:,None,None]
            row('exposure',group,'current_home_parent_suitable',weight,np.sum(w*current_parent),den)
            income_dependent = np.any(((resources['soft'][...,None]>=(1-.8)*price*house)&purchasing)&suitable,axis=-1)&~np.any(((resources['before_current_income'][...,None]>=(1-.8)*price*house)&purchasing)&suitable,axis=-1)
            row('affordability',group,'phi080_parent_purchase_requires_current_income',weight,np.sum(w*income_dependent),den,'soft',.8)
    housing_checks = []
    common_support = {}
    for n,m in ((0,0),(1,1)):
        common_support[n,m] = (a['tenure_probs'][:,:,0,:,:,n,m,:].sum(axis=-1)>0)&(a1['tenure_probs'][:,:,0,:,:,n,m,:].sum(axis=-1)>0)
    for arm, arrays in (('phi080',a),('phi100',a1)):
        for n,m,branch in ((0,0,'childless_wait'),(1,1,'one_child_at_home')):
            # hR is conditional on becoming a renter; owners' net sale enters
            # its wealth argument once through b+S/R, as in the tenure map.
            renter_h = np.empty(G.shape)
            cap = np.empty(G.shape)
            for ten in range(G.shape[1]):
                x = np.clip(b+sale[ten]/R,b[0],b[-1])
                for j in range(p['J']):
                    for z in range(len(a['type_values'])):
                        raw = arrays['hR_pol'][:,0,0,j,z,n,m]
                        renter_h[:,ten,j,z] = np.interp(x,b,raw)
                        cap[:,ten,j,z] = np.interp(x,b,(raw>=p['hR_max']-1e-8).astype(float))
            tp_raw = np.asarray(arrays['tenure_probs'][:,:,0,:,:,n,m,:],dtype=float)
            psum = tp_raw.sum(axis=-1,keepdims=True)
            tp = np.divide(tp_raw,psum,out=np.zeros_like(tp_raw),where=psum>0)
            pr = tp[...,0]
            h = pr*renter_h+np.sum(tp[...,1:]*house,axis=-1)
            hf = float(arrays['shared.h_bar'][n,m])
            dead = tp.sum(axis=-1)==0
            housing_checks.append(dict(arm=arm,branch=branch,dead_common_exposure=float(np.sum(G*dead)),live_tenure_probability_sum_max_error=float(np.max(np.abs(tp.sum(axis=-1)[(G>0)&~dead]-1))),conditional_renter_below_floor_exposure=float(np.sum(G*pr*(renter_h<hf-1e-8))),conditional_renter_above_cap_exposure=float(np.sum(G*pr*(renter_h>p['hR_max']+1e-8)))))
            for group,mask in groups.items():
                for weight,w in (('mass',G),('success_gap_susceptibility',W)):
                    w=w*mask
                    row('housing_'+branch,group,'common_menu_support_share',weight,np.sum(w*common_support[n,m]),w.sum(),phi=.8 if arm=='phi080' else 1.)
                    w=w*common_support[n,m]; den=w.sum(); renter_den=np.sum(w*pr)
                    for metric,num,d in (('expected_physical_rooms',np.sum(w*h),den),('renter_choice_probability',renter_den,den),('conditional_renter_rooms',np.sum(w*pr*renter_h),renter_den),('renter_cap_given_renter_choice',np.sum(w*pr*cap),renter_den),('renter_cap_joint_exposure',np.sum(w*pr*cap),den)):
                        row('housing_'+branch,group,metric,weight,num,d,phi=.8 if arm=='phi080' else 1.)
                    for k,hh in enumerate(house):
                        row('housing_'+branch,group,'owner_product_choice_probability',weight,np.sum(w*tp[...,k+1]),den,phi=.8 if arm=='phi080' else 1.,product=hh)
    # Realized baseline owner saving: new purchasers use collateral floor;
    # DUE stayers use min(b,collateral), plus net-estate solvency when death is possible.
    debt_checks = []
    for group in ('all_ages_all_family_states','fertile_childless'):
        total = binding = stay_total = stay_binding = new_total = new_binding = 0.
        for ten in range(1,6):
            bf = -.8*price*house[ten-1]
            for j in range(p['J']):
                if group=='fertile_childless' and not fertile[j]: continue
                g = a['g'][:,ten,0,j,:,:,:]
                gs = a['g_stay_distribution'][:,ten,0,j,:,:,:]
                bp = a['bp_pol'][:,ten,0,j,:,:,:]
                bps = a['bp_pol_stay'][:,ten,0,j,:,:,:]
                if group=='fertile_childless':
                    g,gs,bp,bps = (x[:, :, 0, 0] for x in (g,gs,bp,bps))
                gn = g-gs
                if gn.min() < -1e-12: raise ValueError('Stayer mass exceeds realized owner mass')
                floor_new = max(bf,b[0])
                death_possible = j==p['J']-1 or (p['use_age_survival'] and p['survival_probs'][j]<1)
                floor_stay = np.maximum(np.minimum(b,bf),-sale[ten] if death_possible else -np.inf)
                floor_stay = np.maximum(floor_stay,b[0]).reshape((len(b),)+(1,)*(g.ndim-1))
                bn = np.abs(bp-floor_new)<=1e-6; bs=np.abs(bps-floor_stay)<=1e-6
                total+=g.sum(); stay_total+=gs.sum();new_total+=gn.sum()
                binding+=np.sum(gn*bn+gs*bs);stay_binding+=np.sum(gs*bs);new_binding+=np.sum(gn*bn)
        for metric,num,den in (('ending_owner_floor_binding_share',binding,total),('stayer_floor_binding_share',stay_binding,stay_total),('new_owner_floor_binding_share',new_binding,new_total)):
            row('realized_baseline_owner_debt',group,metric,'mass',num,den,phi=.8)
        debt_checks.append(dict(group=group,owner_mass=float(total),new_owner_mass=float(new_total),stayer_mass=float(stay_total),floor_tolerance=1e-6))
    with (HERE/'census.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    median22 = int(np.searchsorted(np.cumsum(G[:,:,1,:].sum(axis=(0,1))),G[:,:,1,:].sum()/2))
    assert median22==4
    sources = [dict(path=str(path),sha256=helper.sha(path)) for path in (Path(__file__).resolve(),HELPER,Path('code/model/production/engine/household.py').resolve(),Path('code/model/production/engine/kernels.py').resolve(),Path('code/model/production/engine/shared.py').resolve(),Path('code/model/production/engine/distribution.py').resolve(),BASE/'executed_P.json',BASE/'solution_arrays.npz',RELAXED/'executed_P.json',RELAXED/'solution_arrays.npz')]
    receipt = dict(reference='post-interest soft chain13, no Estate A; phi080 fixed-price baseline exposure',zero_solves=True,changes=changes,price=price,R_gross=R,H_own=house.tolist(),parent_floor=floor,owner_service_premium=p['chi'],renter_cap=p['hR_max'],susceptibility='G_pre*pi_age^2*a*(1-a)/kappa_fert; derivative of first-birth flow w.r.t. successful-birth-minus-wait cardinal utility',timing='soft: Rb+y+S >= (1-phi)Q; diagnostic cash: Rb+S >= (1-phi)Q',income='period after-tax P.income[i,j]*z+zero property-tax transfer; distinct from annual gross earnings',worked_example=dict(age=22,z=float(a['type_values'][4]),income=float(income[1,4]),four_room_purchase_price=float(price*4),four_room_downpayment=float((1-.8)*price*4),minimum_b_renter_soft=float(((1-.8)*price*4-income[1,4])/R),minimum_b_renter_before_income=float((1-.8)*price*4/R)),stayers='same owner product excluded from purchase census; current-home suitability recorded separately',inversion_checks=[inversion,inversion1],housing_checks=housing_checks,debt_checks=debt_checks,sources=sources,rows=len(rows))
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'rows':len(rows),'exposure':float(G.sum()),'susceptibility':float(W.sum()),'housing_checks':housing_checks},indent=2))

if __name__ == '__main__': main()
