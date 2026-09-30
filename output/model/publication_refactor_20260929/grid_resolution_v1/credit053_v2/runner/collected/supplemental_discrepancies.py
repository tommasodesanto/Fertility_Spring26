"""Supplemental read-only discrepancy summary/plots; no model calls or source edits."""
from pathlib import Path
import base64,io,json,shlex,subprocess,os,sys
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',MPLBACKEND='Agg')
HERE=Path(__file__).resolve().parent
SCRIPT=r'''
import numpy as np,json,io,base64
from pathlib import Path
root=Path('/scratch/td2248/projects/grid_resolution_credit053_v2/results/full')
paths=[root/a/'phase_b_ge/selected_root/common_support_policies.npz' for a in ['control_160x15','proposal_120x9']]
gpaths=[root/a/'phase_b_ge/selected_repeat/stage/solution_arrays.npz' for a in ['control_160x15','proposal_120x9']]
with np.load(paths[0]) as a,np.load(paths[1]) as b,np.load(gpaths[1]) as gg:
    x=b['b_grid'];z=b['z_grid'];oz=a['z_grid'];idx=np.searchsorted(a['b_grid'],x);assert np.array_equal(a['b_grid'][idx],x)
    hi=np.clip(np.searchsorted(oz,z),1,len(oz)-1);lo=hi-1;t=((z-oz[lo])/(oz[hi]-oz[lo])).reshape(1,1,1,1,-1,1,1)
    V=b['V'];vlo=a['V'][idx][:,:,:,:,lo];vhi=a['V'][idx][:,:,:,:,hi];g=gg['g_beginning_distribution']
    child=(np.arange(V.shape[6])[None,:]<=np.arange(V.shape[5])[:,None]).reshape((1,1,1,1,1,V.shape[5],V.shape[6]))
    common=(vlo>-1e9)&(vhi>-1e9)&(V>-1e9)&child
    weights=np.where(common,g,0.);total=float(g.sum());covered=float(weights.sum())
    result=dict(weight_definition='Proposal native g_beginning_distribution: post-fertility, pre-current-tenure household mass. All strictly positive weights retained; no minimum occupancy cutoff.',comparison='Exact common wealth-node subset; linear interpolation across old income neighbors; both old stencil V and new V must exceed -1e9; m>n excluded. Each arm at own renewal-clearing price.',total_mass=total,common_covered_mass=covered,common_coverage_share=covered/total,fields={},group_summaries={},fertility_policy='Not analyzed: saved beginning mass is post-fertility and cannot weight fertility policies as pre-fertility mass.',units={'owner_choice_probability':'probability (percentage points in plots)','hR_pol':'conditional renter housing services, model room units','bp_pol':'conditional saving, model wealth units','c_pol':'conditional consumption, model wealth units','V':'conditional value-function units'},price_control=json.load(open(root/'control_160x15/phase_b_ge/selected_root/closure.json'))['price'],price_proposal=json.load(open(root/'proposal_120x9/phase_b_ge/selected_root/closure.json'))['price'])
    slices={'wealth':x,'income':z};cdf={};wealth_edges=np.array([-12,0,.5,1,2,5,10,20,50,100,300,1000,3000.00001]);owner_gap=None;housing_gap=None
    for key,thresholds in [('owner_choice_probability',[.01,.05,.10]),('hR_pol',[.1,.5,1.]),('bp_pol',[]),('c_pol',[]),('V',[])]:
        old=a[key];new=b[key];mapped=(1-t)*old[idx][:,:,:,:,lo]+t*old[idx][:,:,:,:,hi];gap=new-mapped
        w=weights.copy()
        if key=='hR_pol':w[:,1:]=0 # conditional renter margin only
        positive=w>0;values=abs(gap[positive]);mass=w[positive];den=float(mass.sum());order=np.argsort(values);cumulative=np.cumsum(mass[order])/den
        quantiles={str(q):float(values[order][min(np.searchsorted(cumulative,q),len(order)-1)]) for q in [.5,.9,.95,.99,.999]}
        quantile_grid=np.linspace(0,1,101);quantile_values=values[order][np.minimum(np.searchsorted(cumulative,quantile_grid),len(order)-1)]
        cdf[key]=dict(probability=quantile_grid.tolist(),absolute_gap=quantile_values.tolist())
        field=dict(conditioning_mass=den,mean_abs=float(np.sum(values*mass)/den),quantiles=quantiles,thresholds=[dict(threshold=v,mass=float(w[abs(gap)>v].sum()),conditional_share=float(w[abs(gap)>v].sum()/den),share_of_total_households=float(w[abs(gap)>v].sum()/total)) for v in thresholds])
        coordinate=np.unravel_index(np.argmax(np.where(positive,abs(gap),-1)),gap.shape)
        bi,ten,loc,age,zi,n,m=map(int,coordinate)
        field['largest_positive_mass_state']=dict(indices=list(map(int,coordinate)),wealth=float(x[bi]),income_multiplier=float(z[zi]),proposed_value=float(new[coordinate]),mapped_control_value=float(mapped[coordinate]),signed_gap=float(gap[coordinate]),post_fertility_pre_tenure_mass=float(g[coordinate]),children_ever_born=n,children_at_home=m)
        result['fields'][key]=field
        groups={}
        for label,axis in [('age',3),('children_ever_born',5),('beginning_tenure',1)]:
            rows=[]
            for j in range(gap.shape[axis]):
                sl=[slice(None)]*gap.ndim;sl[axis]=j;sl=tuple(sl);ww=w[sl];mm=float(ww.sum());rows.append(dict(index=j,mass=mm,mean_abs=float((abs(gap[sl])*ww).sum()/mm) if mm else None))
            groups[label]=rows
        rows=[]
        for left,right in zip(wealth_edges[:-1],wealth_edges[1:]):
            pick=(x>=left)&(x<right);ww=w[pick];mm=float(ww.sum());rows.append(dict(left=float(left),right=float(right),mass=mm,mean_abs=float((abs(gap[pick])*ww).sum()/mm) if mm else None))
        groups['wealth_bins']=rows;result['group_summaries'][key]=groups
        if key=='owner_choice_probability':owner_gap=gap;owner_mapped=mapped;owner_new=new;max_owner=coordinate
        if key=='hR_pol':housing_gap=gap;housing_mapped=mapped;housing_new=new;max_housing=coordinate
    # Two explicit conditional slices: occupied maximum owner gap and age30, income closest1, childless renter.
    benchmark=(0,0,0,3,int(np.argmin(abs(z-1))),0,0)
    for name,coordinate in [('largest_occupied_owner_gap',max_owner),('benchmark_age30_middle_income_childless_renter',benchmark),('largest_occupied_housing_gap',max_housing)]:
        bi,ten,loc,age,zi,n,m=map(int,coordinate);sel=(slice(None),ten,loc,age,zi,n,m);valid=common[sel]
        slices[name+'_indices']=np.array([ten,loc,age,zi,n,m]);slices[name+'_valid']=valid
        slices[name+'_owner_control']=np.where(valid,owner_mapped[sel],np.nan);slices[name+'_owner_proposal']=np.where(valid,owner_new[sel],np.nan);slices[name+'_mass']=g[sel]
        own=a['owner_choice_probability'];slices[name+'_owner_low_income']=np.where(valid,own[idx,ten,loc,age,lo[zi],n,m],np.nan);slices[name+'_owner_high_income']=np.where(valid,own[idx,ten,loc,age,hi[zi],n,m],np.nan)
        # Renter housing is always the conditional renter branch, preserving other state coordinates.
        rent_sel=(slice(None),0,loc,age,zi,n,m);rent_valid=common[rent_sel]
        slices[name+'_housing_control']=np.where(rent_valid,housing_mapped[rent_sel],np.nan);slices[name+'_housing_proposal']=np.where(rent_valid,housing_new[rent_sel],np.nan);slices[name+'_renter_mass']=g[rent_sel]
        result.setdefault('slices',{})[name]=dict(tenure=ten,location=loc,age_index=age,age_years=18+4*age,income_index=zi,income_multiplier=float(z[zi]),children_ever_born=n,children_at_home=m,owner_slice_mass=float(g[sel].sum()),renter_housing_slice_mass=float(g[rent_sel].sum()))
    result['weighted_quantile_grid']=cdf
    buf=io.BytesIO();np.savez_compressed(buf,**slices)
    print(json.dumps(dict(summary=result,slices_npz=base64.b64encode(buf.getvalue()).decode())))
'''
def main():
    if '--render-only' in sys.argv:
        data=json.loads((HERE/'supplemental_discrepancies.json').read_text())
    else:
        response=subprocess.run(['ssh','torch','/share/apps/anaconda3/2025.06/bin/python -c '+shlex.quote(SCRIPT)],capture_output=True,text=True,check=True)
        packet=json.loads(response.stdout);data=packet['summary'];(HERE/'supplemental_discrepancies.json').write_text(json.dumps(data,indent=2)+'\n');(HERE/'supplemental_slices.npz').write_bytes(base64.b64decode(packet['slices_npz']))
    import csv,numpy as np,matplotlib.pyplot as plt
    with (HERE/'supplemental_discrepancies.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['field','conditioning_mass','mean_absolute_gap','p95','p99']);writer.writeheader()
        for key,row in data['fields'].items():writer.writerow(dict(field=key,conditioning_mass=row['conditioning_mass'],mean_absolute_gap=row['mean_abs'],p95=row['quantiles']['0.95'],p99=row['quantiles']['0.99']))
    fig,axes=plt.subplots(3,2,figsize=(12,11));fig.suptitle('Supplemental grid discrepancies: D = 0.53, each grid at its own GE price',fontsize=13)
    bins=data['group_summaries']['owner_choice_probability']['wealth_bins'];labels=[f"{r['left']:g}–{min(r['right'],3000):g}" for r in bins];vals=[100*r['mean_abs'] if r['mean_abs'] is not None else np.nan for r in bins]
    axes[0,0].bar(range(len(bins)),vals);axes[0,0].set_xticks(range(len(bins)),labels,rotation=50,ha='right',fontsize=8);axes[0,0].set_ylabel('Mean absolute ownership gap (pp)');axes[0,0].set_title('Proposal mass weighted, by wealth bin')
    for j,row in enumerate(bins):
        if row['mean_abs'] is not None:axes[0,0].text(j,vals[j],f"{100*row['mass']/data['total_mass']:.2g}%",ha='center',va='bottom',fontsize=7)
    axes[0,0].margins(y=.25)
    row=data['weighted_quantile_grid']['owner_choice_probability'];cdf_x=np.array(row['probability'])*100;cdf_y=np.array(row['absolute_gap'])*100;keep=(cdf_x>=50)&(cdf_x<=99);axes[0,1].plot(cdf_x[keep],cdf_y[keep])
    axes[0,1].set_xlim(50,99);axes[0,1].set_ylim(0,110*data['fields']['owner_choice_probability']['quantiles']['0.99']);axes[0,1].set_xlabel('Weighted percentile');axes[0,1].set_ylabel('Ownership absolute gap (percentage points)');axes[0,1].set_title('Median through 99th percentile')
    with np.load(HERE/'supplemental_slices.npz') as s:
        x=s['wealth']
        for r,name in enumerate(['largest_occupied_owner_gap','benchmark_age30_middle_income_childless_renter'],start=1):
            for col,kind,ylabel in [(0,'owner','Owner choice probability'),(1,'housing','Conditional renter housing (rooms)')]:
                slice_name='largest_occupied_housing_gap' if r==2 and kind=='housing' else name
                info=data['slices'][slice_name];title=f"Age {info['age_years']}, z={info['income_multiplier']:.3g}; children ever born/home {info['children_ever_born']}/{info['children_at_home']}"
                ax=axes[r,col];mass=s[slice_name+('_mass' if kind=='owner' else '_renter_mass')];ax.plot(x,s[slice_name+'_'+kind+'_control'],label='160×15 interpolated',lw=1.6);ax.plot(x,s[slice_name+'_'+kind+'_proposal'],label='120×9',lw=1.6)
                if kind=='owner' and r==1:
                    ax.plot(x,s[slice_name+'_owner_low_income'],':',color='tab:blue',alpha=.5,label='Fine lower-income neighbor');ax.plot(x,s[slice_name+'_owner_high_income'],'--',color='tab:blue',alpha=.5,label='Fine upper-income neighbor')
                ax.set_xlim((-5,5) if slice_name=='largest_occupied_housing_gap' else ((-.1,.7) if kind=='owner' and r==1 else (-1,20)));
                if kind=='owner' and r==1:title+=' (wealth zoom)'
                ax.set_ylabel(ylabel);ax.set_xlabel('Wealth (model units)');ax.set_title(title+f"; slice mass={mass.sum():.3g}",fontsize=9);ax.legend(fontsize=7)
                twin=ax.twinx();twin.fill_between(x,0,mass,color='grey',alpha=.15);twin.set_ylabel('Beginning cell mass',fontsize=8);twin.tick_params(labelsize=7)
    fig.text(.02,.012,'Weight: post-fertility, pre-tenure proposal mass. Housing/saving/consumption are conditional policies; these are not realized-policy accuracy measures. No occupancy cutoff. Standard plots unchanged.',fontsize=8)
    fig.tight_layout(rect=[0,.035,1,.97]);fig.savefig(HERE/'supplemental_discrepancies.png',dpi=160);fig.savefig(HERE/'supplemental_discrepancies.pdf');plt.close(fig)
    print(json.dumps(dict(common_coverage=data['common_coverage_share'],metrics={k:dict(mean=v['mean_abs'],p95=v['quantiles']['0.95'],p99=v['quantiles']['0.99']) for k,v in data['fields'].items()})))
if __name__=='__main__':main()
