from pathlib import Path
import sys,json,time
import numpy as np
root=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907j')
sys.path[:0]=[str(root/'code/model'),str(root/'code/model/tools')]
import run_e5f_independent_numerical_audit as audit
import run_e5f_joint_overnight_case as a
source=root/'output/model/joint_nested_overnight/policy_loop_smoke/supply-plus-20/date_2023/dated_state.pkl.gz'
packet=audit.load_checkpoint(source);P=packet['parameters'];e=packet['evaluation'];p=e.policy;bg=packet['b_grid'];sd=packet['shared']
out=Path(__file__).resolve().parent;rows=[]
for age in (58,62):
 j=int(round((age-P.age_start)/P.da));badb=int(np.argmin(abs(bg+.25581395348837255)));n=c=3;z=0
 idx=(badb,0,0,j,z,n,c);contributions=[]
 for ten in range(e.g_post_fertility.shape[1]):
  for b in range(bg.size):
   mass=float(e.g_post_fertility[b,ten,0,j,z,n,c]);prob=float(p.tenure_probs[b,ten,0,j,z,n,c,0])
   if mass*prob<=0:continue
   # With one market there is no location transaction.
   kl=b;wl=0.0
   for bl,lw in ((kl,1-wl),(kl+1,wl)):
    if lw<=0:continue
    kt=int(p.maps.tmx_idx[0,ten,0,n,c,bl]);wt=float(p.maps.tmx_wt[0,ten,0,n,c,bl])
    for bt,tw in ((kt,1-wt),(kt+1,wt)):
     if bt!=badb or tw<=0:continue
     contributions.append(dict(origin_b=b,wealth=float(bg[b]),tenure=ten,post_fertility_mass=mass,renter_probability=prob,location_weight=lw,transaction_lower_node=kt,transaction_weight_upper=wt,destination_weight=tw,mass=mass*prob*lw*tw,
        origin_joint_value=float(p.V[b,ten,0,j,z,n,c]),origin_joint_probabilities=p.joint_choice.probabilities[b,ten,0,j,z,n,c].tolist(),
        destination_joint_values=[float(p.V[kt,0,0,j,z,n,c]),float(p.V[kt+1,0,0,j,z,n,c])],destination_consumptions=[float(p.c_pol[kt,0,0,j,z,n,c]),float(p.c_pol[kt+1,0,0,j,z,n,c])],destination_housing=[float(p.hR_pol[kt,0,0,j,z,n,c]),float(p.hR_pol[kt+1,0,0,j,z,n,c])]))
 total=sum(r['mass'] for r in contributions);expected=float(e.g_current[idx]);assert abs(total-expected)<1e-24,(total,expected)
 rows.append(dict(age=age,current_index=list(idx),current_mass=expected,current_value=float(p.V[idx]),current_consumption=float(p.c_pol[idx]),current_housing=float(p.hR_pol[idx]),current_saving=float(p.bp_pol[idx]),contributions=contributions))
proof=dict(checkpoint=str(source),checkpoint_sha256=a.digest(source),parameters=dict(price=p.price.tolist(),house_grid=P.H_own.tolist(),sale_cost=P.psi,financed_share=P.phi if np.isscalar(P.phi) else np.asarray(P.phi).tolist()),rows=rows,feasibility_projection_mass=e.feasibility_projection_mass,status='original_current_mass_exactly_traced',production_promoted=False)
a.write_json(out/'original_mass_trace.json',proof);print(json.dumps(proof),flush=True)
