"""Small cached pair table: realized rooms/ownership and susceptible birth gap."""
from pathlib import Path
import argparse,json
import numpy as np
from extract_common_states import load,DEFAULT

def main():
    ap=argparse.ArgumentParser();ap.add_argument('alternative',type=Path);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    p,a,g,h,ages,fec,checks=load(DEFAULT/'phi_080')
    q,c,k,l,_,_,other=load(args.alternative)
    a0=a['fert_probs'][...,1];a1=c['fert_probs'][...,1];pw0=a['fert_probs'][...,0];pw1=c['fert_probs'][...,0];kappa=p['kappa_fert']
    pi=fec[None,None,None,:,None];w=g[...,0,0];fertile=np.arange(p['J'])<p['A_f_end']
    score=w*pi*pi*a0*pw0/kappa;score[:,:,:,~fertile,:]=0
    valid=(a0>0)&(pw0>0)&(a1>0)&(pw1>0)&(pi>0)
    with np.errstate(divide='ignore',invalid='ignore'):
        dD=kappa*(np.log(a1)-np.log(pw1)-np.log(a0)+np.log(pw0))/pi
    gap={}
    for label,mask in [('all',np.ones(w.shape,dtype=bool)),('inherited_renter',np.broadcast_to(np.arange(w.shape[1])[None,:,None,None,None]==0,w.shape)),('b_nonpositive',np.broadcast_to(a['b_grid'][:,None,None,None,None]<=0,w.shape))]:
        total=float(score[mask].sum());supported=valid&mask
        gap[label]=dict(susceptibility_sum=total,supported_susceptibility=float(score[supported].sum()),unsupported_susceptibility=float(score[mask&~valid].sum()),weighted_success_gap_delta=float(np.sum(score[supported]*dD[supported])/score[supported].sum()),first_birth_flow_linear_change=float(np.sum(score[supported]*dD[supported])),definition='baseline G*pi^2*a*(1-a)/kappa weights; derivative of realized first births with respect to successful-birth-minus-wait utility')
    rows=[]
    for label,path,P in [('baseline',DEFAULT/'phi_080',p),('alternative',args.alternative,q)]:
        with np.load(path/'solution_arrays.npz',allow_pickle=False) as z:
            mass=z['g'];rent=z['hR_pol'];t=mass.sum(axis=(0,2,4,5,6));age_mass=t.sum(0)
            owner_mass=t[1:].sum(0);rooms=(np.asarray(P['H_own'])[:,None]*t[1:]).sum(0)+(mass[:,0]*rent[:,0]).sum(axis=(0,1,3,4,5))
        for j in range(7):rows.append(dict(case=label,age=int(ages[j]),mass=float(age_mass[j]),ownership=float(owner_mass[j]/age_mass[j]),rooms=float(rooms[j]/age_mass[j])))
        js=[0,1,2];rows.append(dict(case=label,age='18-29 pooled',mass=float(age_mass[js].sum()),ownership=float(owner_mass[js].sum()/age_mass[js].sum()),rooms=float(rooms[js].sum()/age_mass[js].sum())))
    result=dict(age_table=rows,first_birth_gap=gap,checks=[checks,other])
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
