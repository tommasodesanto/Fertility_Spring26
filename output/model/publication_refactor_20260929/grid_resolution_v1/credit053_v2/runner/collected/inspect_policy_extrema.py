"""Read-only remote saved-policy extrema/occupancy; no numerical model solve."""
import json,shlex,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
script=r'''
import json,time,numpy as np
from pathlib import Path
root=Path('/scratch/td2248/projects/grid_resolution_credit053_v2/results/full');start=time.monotonic()
paths=[root/n/'phase_b_ge/selected_root/common_support_policies.npz' for n in ['control_160x15','proposal_120x9']]
gpaths=[root/n/'phase_b_ge/selected_repeat/stage/solution_arrays.npz' for n in ['control_160x15','proposal_120x9']]
with np.load(paths[0]) as a,np.load(paths[1]) as b,np.load(gpaths[0]) as ga,np.load(gpaths[1]) as gb:
    wealth=b['b_grid'];income=b['z_grid'];oldincome=a['z_grid'];idx=np.searchsorted(a['b_grid'],wealth);assert np.array_equal(a['b_grid'][idx],wealth)
    hi=np.clip(np.searchsorted(oldincome,income),1,len(oldincome)-1);lo=hi-1;w=((income-oldincome[lo])/(oldincome[hi]-oldincome[lo])).reshape(1,1,1,1,-1,1,1)
    av=a['V'];bv=b['V'];vlo=av[idx][:,:,:,:,lo];vhi=av[idx][:,:,:,:,hi]
    child=np.arange(bv.shape[6])[None,:]<=np.arange(bv.shape[5])[:,None];child=child.reshape((1,1,1,1,1)+child.shape)
    common=(vlo>-1e9)&(vhi>-1e9)&(bv>-1e9)&child
    ag=ga['g_beginning_distribution'];bg=gb['g_beginning_distribution'];assert ag.shape==av.shape and bg.shape==bv.shape
    oldlo=ag[idx][:,:,:,:,lo];oldhi=ag[idx][:,:,:,:,hi];occupied=common&(bg>0);bothoccupied=occupied&((oldlo>0)|(oldhi>0))
    result=dict(axis_order=['wealth','tenure','location','age','income','children','child_state'],distribution_source='native solution.g_beginning_distribution from selected exact-repeat stage: beginning/pre-current-tenure distribution, raw cell probability, not density. solution.g is postdecision current mass and is not used for these occupancy weights',policy_source='selected_root snapshot, each arm selected own renewal price',common_state_count=int(common.sum()),proposal_occupied_common_state_count=int(occupied.sum()),both_occupied_common_state_count=int(bothoccupied.sum()),proposal_beginning_g_total=float(bg.sum()),control_beginning_g_total=float(ag.sum()),fields={})
    def cell_info(key,gap,mask):
        k=np.unravel_index(np.argmax(np.where(mask,abs(gap),-1)),gap.shape);bi,t,l,age,z,n,m=map(int,k)
        oldk=(int(idx[bi]),t,l,age,int(lo[z]),n,m);oldkh=(int(idx[bi]),t,l,age,int(hi[z]),n,m)
        aa=a[key];bb=b[key];mapped=float((1-w[0,0,0,0,z,0,0])*aa[oldk]+w[0,0,0,0,z,0,0]*aa[oldkh])
        return dict(indices=list(map(int,k)),wealth=float(wealth[bi]),income_multiplier=float(income[z]),control_income_neighbors=[float(oldincome[lo[z]]),float(oldincome[hi[z]])],control_income_weight_high=float(w[0,0,0,0,z,0,0]),absolute_gap=float(abs(gap[k])),proposed_value=float(bb[k]),mapped_control_value=mapped,control_value_low=float(aa[oldk]),control_value_high=float(aa[oldkh]),proposal_g_mass=float(bg[k]),control_g_mass_low=float(ag[oldk]),control_g_mass_high=float(ag[oldkh]),proposal_V=float(bv[k]),control_V_low=float(av[oldk]),control_V_high=float(av[oldkh]))
    for key in ['V','bp_pol','c_pol','hR_pol','owner_choice_probability']:
        value=a[key];mapped=(1-w)*value[idx][:,:,:,:,lo]+w*value[idx][:,:,:,:,hi];gap=b[key]-mapped
        result['fields'][key]=dict(common_maximum=cell_info(key,gap,common),proposal_occupied_maximum=cell_info(key,gap,occupied),both_occupied_maximum=cell_info(key,gap,bothoccupied),proposal_g_weighted_mean_abs=float(np.sum(abs(gap[occupied])*bg[occupied])/np.sum(bg[occupied])),proposal_g_common_covered_mass=float(np.sum(bg[occupied])))
    result['remote_readonly_seconds']=time.monotonic()-start
    print(json.dumps(result,indent=2))
'''
result=subprocess.run(['ssh','torch','/share/apps/anaconda3/2025.06/bin/python -c '+shlex.quote(script)],capture_output=True,text=True,check=True)
packet=json.loads(result.stdout);(HERE/'policy_extrema_occupancy.json').write_text(json.dumps(packet,indent=2)+'\n')
print(json.dumps({k:dict(max=row['common_maximum']['absolute_gap'],max_proposal_mass=row['common_maximum']['proposal_g_mass'],occupied_max=row['proposal_occupied_maximum']['absolute_gap'],weighted_mean=row['proposal_g_weighted_mean_abs']) for k,row in packet['fields'].items()},indent=2))
