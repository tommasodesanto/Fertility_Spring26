"""Diagnose zero occupied saving responses without changing the compiled fixture."""
import json
from pathlib import Path
import numpy as np
from test_e5f_social_security_compiled import CompiledSocialSecurityTests as T

T.setUpClass()
rows=[]
for changed_name, control_name in [('future_pension','fixed_tax'),('future_tax','fixed_pension')]:
    changed, control=T.results[changed_name],T.results[control_name]
    left,right=changed['dated'][0]['evaluation'],control['dated'][0]['evaluation']
    g=right.g_current
    occupied=g>1e-12
    bp=right.policy.bp_pol
    row=dict(case=changed_name,parameters={k:float(getattr(T.parameters,k)) for k in
        ('beta','R_gross','pension','tau_pay','b_min','b_max','lambda_d')},
        occupied_mass=float(g[occupied].sum()),
        at_zero_saving_mass=float(g[(bp<=1e-10)&occupied].sum()),
        positive_saving_mass=float(g[(bp>1e-10)&occupied].sum()),
        occupied_min_saving=float(bp[occupied].min()),occupied_max_saving=float(bp[occupied].max()),
        current_value_change_occupied=float(np.max(np.abs(changed['path'].values[0]-control['path'].values[0])[T.initial_g>1e-12])),
        current_value_change_all=float(np.max(np.abs(changed['path'].values[0]-control['path'].values[0]))),
        future_value_wealth_slope_change=float(np.max(np.abs(np.diff(changed['path'].values[1]-control['path'].values[1],axis=0)))))
    for field in ('bp_pol','c_pol','hR_pol','Pi_pol','fert2_probs'):
        if not hasattr(left.policy,field): continue
        a,b=np.asarray(getattr(left.policy,field)),np.asarray(getattr(right.policy,field))
        diff=np.abs(a-b)
        key=field+'_max_change_all'
        row[key]=float(diff.max())
        if diff.shape==g.shape:
            row[field+'_max_change_occupied']=float(diff[occupied].max())
            index=np.unravel_index(int(diff.argmax()),diff.shape)
            row[field+'_max_change_state']=list(map(int,index))
            row[field+'_control_at_max_change']=float(b[index])
            row[field+'_changed_at_max_change']=float(a[index])
            row[field+'_mass_at_max_change']=float(g[index])
    row['ages']=[]
    for age in range(int(T.parameters.J)):
        gm=g[:,:,:,age,:,:,:]
        bm=bp[:,:,:,age,:,:,:]
        dm=np.abs(left.policy.bp_pol-right.policy.bp_pol)[:,:,:,age,:,:,:]
        vm=np.abs(changed['path'].values[0]-control['path'].values[0])[:,:,:,age,:,:,:]
        occupied_age=gm>1e-12
        row['ages'].append(dict(age_index=age,mass=float(gm.sum()),
            zero_saving_mass=float(gm[bm<=1e-10].sum()),
            max_saving_change_all=float(dm.max()),
            max_value_change_all=float(vm.max()),
            max_saving_change_occupied=float(dm[occupied_age].max()) if occupied_age.any() else None))
    rows.append(row)
packet=dict(source_commit='a654219c',scope='same six conditional paths; inspect saving corners and off-support responses',rows=rows)
Path('output/anticipation_diagnostic.json').write_text(json.dumps(packet,indent=2,allow_nan=False)+'\n')
print(json.dumps(packet,allow_nan=False),flush=True)
