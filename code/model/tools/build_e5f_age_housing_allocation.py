"""Extract saved allocation by age, then render supplemental seminar panels."""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,pickle,sys
import numpy as np

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--checkpoint',type=Path);ap.add_argument('--source-root',type=Path)
    ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    if a.checkpoint:
        sys.path[:0]=[str(a.source_root/'code/model/tools'),str(a.source_root/'code/model')]
        with gzip.open(a.checkpoint,'rb') as stream:q=pickle.load(stream)
        P=q['parameters'];e=q['evaluation'];g=e.g_current;rows=[]
        for j in range(P.J):
            mass=0.;rooms=0.;owners=0.
            for z in range(P.Nz):
                for tenure in range(1+P.n_house):
                    for n in range(P.n_parity):
                        for d in range(P.n_child_states):
                            w=g[:,tenure,0,j,z,n,d];m=float(w.sum());mass+=m
                            if tenure:owners+=m;rooms+=m*float(P.H_own[tenure-1])
                            else:rooms+=float(np.sum(w*e.policy.hR_pol[:,tenure,0,j,z,n,d]))
            rows.append(dict(age=float(P.age_start+j*P.da),households=mass,rooms=rooms,owners=owners,rooms_per_head=rooms/mass,ownership=owners/mass))
        payload=dict(rows=rows,checkpoint=str(a.checkpoint),checkpoint_sha256=hashlib.sha256(a.checkpoint.read_bytes()).hexdigest(),
            model='Selected provisional overnight initial stationary economy; not a transition or policy effect',
            units='Occupied physical rooms, uncapped; actual post-choice household distribution',
            income_states=int(P.Nz),permanent_levels_enabled=bool(P.permanent_income_levels_enabled),
            income_group_index=P.permanent_income_group_index.tolist(),child_state_mode=P.child_state_mode)
        (a.output/'allocation.json').write_text(json.dumps(payload,indent=2)+'\n');return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d=json.loads((a.output/'allocation.json').read_text());rows=d['rows'];age=np.array([r['age'] for r in rows]);mass=np.array([r['households'] for r in rows]);rooms=np.array([r['rooms'] for r in rows])
    with (a.output/'allocation.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    fig,ax=plt.subplots(1,2,figsize=(11.5,4.7))
    ax[0].plot(age,rooms/mass,'o-',color='#145b83');ax[0].set(xlabel='Household-head age',ylabel='Physical rooms per household',title='Housing occupied across the life cycle')
    groups=[(18,34),(34,50),(50,66),(66,86)];labels=['18–33','34–49','50–65','66–85'];hh=[];h=[]
    for lo,hi in groups:
        mask=(age>=lo)&(age<hi);hh.append(100*mass[mask].sum()/mass.sum());h.append(100*rooms[mask].sum()/rooms.sum())
    x=np.arange(4);ax[1].bar(x-.18,hh,.36,label='Share of households',color='#9daeb7');ax[1].bar(x+.18,h,.36,label='Share of occupied rooms',color='#145b83')
    ax[1].set(xticks=x,xticklabels=labels,xlabel='Household-head age group',ylabel='Percent of total',title='How the housing stock is allocated');ax[1].legend(frameon=False,fontsize=9)
    for p in ax:p.spines[['top','right']].set_visible(False);p.grid(axis='y',alpha=.15)
    fig.suptitle('Housing allocation by age',fontsize=16)
    fig.text(.5,.035,'Provisional initial stationary calibration. Distribution-weighted physical rooms; no causal or policy comparison.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.08,1,.94));fig.savefig(a.output/'age_housing_allocation.png',dpi=180);fig.savefig(a.output/'age_housing_allocation.pdf');plt.close(fig)
    print(json.dumps(dict(mean_rooms=float(rooms.sum()/mass.sum()),group_household_shares=hh,group_room_shares=h)))

if __name__=='__main__':main()
