"""Render authenticated saved 2007 stationary comparisons on Torch; zero solves."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path('/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project')
BASE = ROOT/'output/model/fertility_identification_20260928'
NIGHT = BASE/'two_stream_overnight_v1'
OUT = NIGHT/'comparison_v1'
CASES = {
    'Original selected': NIGHT/'run_v1/one_birth/one_birth_024_gn1_0/case',
    'Two-birth selected': NIGHT/'run_v1/two_birth/two_birth_024_gn1_0/case',
}
REF = BASE/'resume_v1/selected_export/primary'
AUDIT = BASE/'measurement_audit_v1'
ACS = OUT/'input/actual2007_age_housing_levels.csv'
ACS_SHA = '9b4433299f43c1c23b8aa28ba224a41ba6ee21d1a7a2265a336f90c831ecb8b9'
NEEDED = ('observers.json','lifecycle_2023.csv','receipt.json','target_fit.csv','parameters.csv')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def js(path):
    return json.loads(path.read_text())

def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))

def write_rows(name, records):
    path = OUT/name
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(records)
    return path

def eq(x,y,tag,tol=1e-10):
    assert math.isfinite(float(x)) and abs(float(x)-float(y))<=tol,(tag,x,y)

def project(accounting,lo,hi):
    ages=np.asarray(accounting['age_cell_start'],float)
    pre=np.asarray(accounting['pre_parity_mass_by_age'],float)
    post=np.asarray(accounting['post_parity_mass_by_age'],float)
    assert pre.shape==post.shape==(len(ages),4)
    assert np.array_equal(ages,18+4*np.arange(len(ages)))
    assert np.max(np.abs(pre.sum(axis=1)-post.sum(axis=1)))<2e-10
    left=np.maximum(ages,lo); right=np.minimum(ages+4,hi+1)
    overlap=np.maximum(right-left,0)/4
    eq(4*overlap.sum(),hi+1-lo,'full window coverage')
    alpha=((left+right)/2-ages)/4
    mass=(overlap[:,None]*((1-alpha[:,None])*pre+alpha[:,None]*post)).sum(axis=0)
    assert mass.sum()>0 and np.all(mass>=0)
    share=mass/mass.sum(); mothers=float(1-share[0]); count=float(share@np.arange(4))
    return {'children_capped3':count,'mother_share':mothers,
            'children_among_mothers_capped3':count/mothers,
            'share0':float(share[0]),'share1':float(share[1]),
            'share2':float(share[2]),'share3plus':float(share[3])}

def main():
    assert os.environ.get('SLURM_JOB_ID','').isdigit() and os.uname().sysname=='Linux'
    OUT.mkdir(parents=True,exist_ok=True)
    inputs={}; observers={}; fits={}; profiles={}; identities={}
    for label,case in {'Reference':REF,**CASES}.items():
        if label=='Reference': hashes=js(case/'artifact_hashes.json')
        else:
            success=js(case.parent/'SUCCESS.json')
            assert success['status']=='passed' and success['candidate_id']==case.parent.name
            hashes={Path(item['path']).relative_to(Path(success['case_path'])).as_posix():item['sha256'] for item in success['artifacts']}
        for name in NEEDED:
            path=case/name
            assert sha(path)==hashes[name],(label,name)
            inputs[str(path)]=sha(path)
        observers[label]=js(case/'observers.json')['fertility']['uniform_birth_time']
        profiles[label]=rows(case/'lifecycle_2023.csv')
        fits[label]={r['moment']:r for r in rows(case/'target_fit.csv')}
        identities[label]=js(case/'receipt.json')
        assert len(fits[label])==14 and len(profiles[label])==17
        assert abs(identities[label]['normalization']['completed_fertility']-2.1)<=5e-4
        eq(sum(float(r['loss_contribution'] or 0) for r in fits[label].values()),identities[label]['loss'],('loss sum',label),1e-8)
    assert identities['Original selected']['target_fingerprint']==identities['Two-birth selected']['target_fingerprint']
    assert identities['Original selected']['config_sha256']==identities['Two-birth selected']['config_sha256']
    assert identities['Original selected']['lane']=='one_birth' and identities['Two-birth selected']['lane']=='two_birth'
    empirical=rows(AUDIT/'fertility_lifecycle_matched_windows.csv')
    first=rows(AUDIT/'first_birth_age_cells.csv')
    assert sha(ACS)==ACS_SHA,'pinned local ACS 2007 input differs after exact transfer'
    acs=rows(ACS)
    for p in (AUDIT/'fertility_lifecycle_matched_windows.csv',AUDIT/'first_birth_age_cells.csv',ACS): inputs[str(p)]=sha(p)
    assert [(int(r['age_lower']),int(r['age_upper'])) for r in empirical]==[(20,24),(25,29),(30,34),(35,39),(40,44)]
    lifecycle=[]
    fmap={'children_capped3':'capped3','mother_share':'mother_share','children_among_mothers_capped3':'given_mother'}
    for row in empirical:
        lo=int(row['age_lower']);hi=int(row['age_upper'])
        for label,observer in observers.items():
            vals=project(observer['accounting'],lo,hi)
            for metric,source in fmap.items():
                if label=='Reference': eq(vals[metric],row['model_'+source],('reference replay',lo,metric))
                lifecycle.append(dict(age_lower=lo,age_upper=hi,metric=metric,series=label,value=vals[metric],source='saved stationary observer'))
        for metric,source in fmap.items():
            lifecycle.append(dict(age_lower=lo,age_upper=hi,metric=metric,series='CPS 2004/2006',value=float(row['data_'+source]),source='pooled CPS fertility supplement'))
    age25=[]
    for label,observer in observers.items():
        v=project(observer['accounting'],25,25)
        eq(v['children_capped3'],fits[label]['early_fertility']['model'],('age25 target',label))
        late=project(observer['accounting'],40,44)
        eq(1-late['mother_share'],fits[label]['cps_childlessness']['model'],('40–44 childlessness',label))
        eq(late['share1']/late['mother_share'],fits[label]['cps_exactly_one']['model'],('40–44 one child among mothers',label))
        age25.append(dict(series=label,**v))
    early=js(AUDIT/'early_fertility_decomposition.json')
    inputs[str(AUDIT/'early_fertility_decomposition.json')]=sha(AUDIT/'early_fertility_decomposition.json')
    age25.append(dict(series='CPS 2004/2006',children_capped3=float(early['data_early_fertility']),mother_share=float(early['data_mother_share']),children_among_mothers_capped3=float(early['data_capped_children_given_mother']),share0='',share1='',share2='',share3plus=''))
    first_rows=[]
    for i,row in enumerate(first):
        lo=int(row['age_lower']);hi=int(row['age_upper'])
        assert (lo,hi)==(18+4*i,21+4*i)
        first_rows.append(dict(age_lower=lo,age_upper=hi,series='NCHS 2003–2006',share=float(row['empirical_share']),tail_mapping=row['empirical_tail_mapping']))
        for label,observer in observers.items():
            a=observer['accounting'];eq(a['age_cell_start'][i],lo,'birth cell')
            share=a['parity_birth_flows_by_age'][i][0]/a['first_birth_flow']
            if label=='Reference':eq(share,row['model_share'],'reference first-birth replay')
            first_rows.append(dict(age_lower=lo,age_upper=hi,series=label,share=share,tail_mapping='model four-year cell'))
    for label in ('NCHS 2003–2006','Reference',*CASES):
        sub=[r for r in first_rows if r['series']==label];eq(sum(r['share'] for r in sub),1,('first shares',label))
        mean=sum((r['age_lower']+2)*r['share'] for r in sub)
        late=sum(r['share'] for r in sub if r['age_lower']>=30)
        field='target' if label=='NCHS 2003–2006' else 'model'
        fit=fits['Reference'] if label=='NCHS 2003–2006' else fits[label]
        eq(mean,fit['nchs_mean_age'][field],('mean first age',label))
        eq(late,fit['nchs_share30'][field],('share first age30+',label))
    costs=[]
    for moment,one in fits['Original selected'].items():
        two=fits['Two-birth selected'][moment]
        assert all(one[k]==two[k] for k in ('role','target','weight'))
        costs.append(dict(moment=moment,role=one['role'],target=one['target'],weight=one['weight'],original=one['model'],two_birth=two['model'],original_gap=one['gap'],two_birth_gap=two['gap'],original_loss=one['loss_contribution'],two_birth_loss=two['loss_contribution'],loss_change=float(two['loss_contribution'] or 0)-float(one['loss_contribution'] or 0)))
    weak=[]
    for lane,label in [('one_birth','Original selected'),('two_birth','Two-birth selected')]:
        path=NIGHT/'run_v1'/lane/'jacobian_1.json';j=js(path);inputs[str(path)]=sha(path)
        matrix=np.asarray(j['matrix'],float);scales=np.asarray(j['coordinate_scales'],float)
        assert matrix.shape==(10,10) and scales.shape==(10,) and np.all(scales>0)
        _,s,vh=np.linalg.svd(matrix,full_matrices=False)
        np.testing.assert_allclose(s,j['diagnostics']['singular_values'],rtol=1e-7,atol=1e-9)
        rank=int(np.count_nonzero(s/s[0]>1e-6));assert rank==j['diagnostics']['numerical_rank']
        v=vh[-1].copy();v*=np.sign(v[np.argmax(np.abs(v))])
        assert np.linalg.norm(matrix@v)-s[-1]<1e-8
        for index,d in enumerate(j['design']):
            weak.append(dict(lane=lane,round_center_candidate=j['center']['candidate_id'],selected_candidate=CASES[label].parent.name,
                parameter=d['parameter'],scaled_coefficient=float(v[index]),physical_displacement_per_unit_vector=float(v[index]*scales[index]),coordinate_scale=float(scales[index]),
                smallest_singular_value=float(s[-1]),rank_at_relative_1e_minus_6=rank,
                caveat='Round-center one-sided finite differences; not a derivative at selected point; scale and step dependent'))
    demo=[]
    for label,r in identities.items():
        demo.append(dict(series=label,loss=r['loss'],completed_fertility=r['normalization']['completed_fertility'],renewal_residual=r['adult_entry_gate']['entry_residual'],psi_child=r['normalization']['psi_child']))
    housing=[]
    acs_by_age={int(r['age_lower']):r for r in acs}
    assert len(acs_by_age)==17
    for label,profile in profiles.items():
        for r in profile:
            lo=int(float(r['age_node']));assert lo in acs_by_age
            owner=float(r['owner_rate']);rooms=float(r['mean_rooms']);wealth=float(r['mean_liquid_wealth'])
            assert 0<=owner<=1 and all(map(math.isfinite,(rooms,wealth)))
            for metric,value in [('homeownership',owner),('rooms_uncapped',rooms),('liquid_wealth_model_units',wealth)]:
                housing.append(dict(age_lower=lo,age_upper=lo+3,metric=metric,series=label,value=value,source='saved stationary lifecycle_2023.csv (inherited filename)'))
    for lo,r in acs_by_age.items():
        housing.append(dict(age_lower=lo,age_upper=lo+3,metric='homeownership',series='ACS 2007',value=float(r['ownership_rate']),source='ACS 2007 household head; 42 metro sample'))
    outputs=[write_rows('plotted_data.csv',lifecycle+[])]
    outputs.extend([write_rows('age25_decomposition_inputs.csv',age25),write_rows('first_birth_age_cells.csv',first_rows),write_rows('other_lifecycle_data.csv',housing),write_rows('target_costs.csv',costs),write_rows('demographic_comparison.csv',demo),write_rows('weak_direction_round_centers.csv',weak)])
    plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
    colors={'CPS 2004/2006':'#bc6043','NCHS 2003–2006':'#bc6043','ACS 2007':'#bc6043','Reference':'#777777','Original selected':'#17678a','Two-birth selected':'#338b65'}
    styles={'CPS 2004/2006':'s--','NCHS 2003–2006':'s--','ACS 2007':'s--','Reference':'o:','Original selected':'o-','Two-birth selected':'^-'}
    fig,axes=plt.subplots(2,2,figsize=(12.8,8.5));axes=axes.flat
    order=['CPS 2004/2006','Original selected','Two-birth selected']
    titles=[('children_capped3','Children ever born per woman (capped at 3)'),('mother_share','Women who are mothers'),('children_among_mothers_capped3','Children ever born among mothers (capped at 3)')]
    for ax,(metric,title) in zip(axes[:3],titles):
        for label in order:
            rr=[r for r in lifecycle if r['metric']==metric and r['series']==label]
            ax.plot(range(5),[r['value'] for r in rr],styles[label],color=colors[label],lw=1.9,ms=5,label='Data (CPS/NCHS)' if label=='CPS 2004/2006' else label)
        ax.set(title=title,xticks=range(5),xticklabels=[f'{r["age_lower"]}–{r["age_upper"]}' for r in rr],xlabel='Age at interview')
        ax.grid(axis='y',alpha=.2)
    ax=axes[3]
    for label in ['NCHS 2003–2006','Original selected','Two-birth selected']:
        rr=[r for r in first_rows if r['series']==label]
        ax.plot(range(7),[100*r['share'] for r in rr],styles[label],color=colors[label],lw=1.9,ms=5,label=label)
    ax.set(title='Share of first births by age cell',xticks=range(7),xticklabels=[f'{r["age_lower"]}–{r["age_upper"]}' for r in rr],xlabel='Age cell',ylabel='Percent of first births');ax.grid(axis='y',alpha=.2)
    fig.suptitle('Fertility lifecycle and first-birth timing',fontsize=16)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=4,frameon=False,bbox_to_anchor=(.5,.055))
    fig.text(.5,.012,'Stationary 2007 model approximation; CPS 2004/2006 cross-sections; NCHS 2003–2006 period first births. Exact age-25 target is separate from the 25–29 bin.\nFirst-birth data tails 12–21 and 42–49 are folded into endpoint model cells. Both selected models normalize completed fertility to 2.1.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=(0,.10,1,.95));p=OUT/'fertility_lifecycle_comparison.png';fig.savefig(p,dpi=165);plt.close(fig);outputs.append(p)
    fig,axes=plt.subplots(1,3,figsize=(14.2,4.7))
    specs=[('homeownership','Homeownership','Share of household heads'),('rooms_uncapped','Rooms, model only','Uncapped rooms'),('liquid_wealth_model_units','Liquid wealth, model only','Model income units')]
    for ax,(metric,title,ylabel) in zip(axes,specs):
        labels=['ACS 2007','Original selected','Two-birth selected'] if metric=='homeownership' else list(CASES)
        for label in labels:
            rr=[r for r in housing if r['metric']==metric and r['series']==label]
            ax.plot([r['age_lower']+1.5 for r in rr],[r['value'] for r in rr],styles[label],color=colors[label],lw=1.9,ms=3.5,label=label)
        ax.set(title=title,xlabel='Age of household head',ylabel=ylabel,xticks=[20,35,50,65,80]);ax.grid(axis='y',alpha=.2)
    axes[0].legend(frameon=False,fontsize=8)
    fig.suptitle('Other lifecycle outcomes of the selected stationary economies',fontsize=15)
    fig.text(.5,.012,'ACS 2007 ownership uses household heads in 42 metros. Available ACS rooms are capped at 9 and cannot be overlaid on saved uncapped model rooms.\nNo matched empirical age profile for liquid wealth is included. Inherited lifecycle_2023.csv is a stationary 2007 export, not a 2023 transition.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=(0,.13,1,.92));p=OUT/'other_lifecycle_comparison.png';fig.savefig(p,dpi=165);plt.close(fig);outputs.append(p)
    original={r['series']:r for r in age25}
    decomposition=[]
    for left,right in [('Original selected','CPS 2004/2006'),('Two-birth selected','CPS 2004/2006'),('Two-birth selected','Original selected')]:
        a=original[left];b=original[right];s1=a['mother_share'];s0=b['mother_share'];k1=a['children_among_mothers_capped3'];k0=b['children_among_mothers_capped3']
        e=(s1-s0)*(k1+k0)/2; i=(k1-k0)*(s1+s0)/2
        eq(e+i,a['children_capped3']-b['children_capped3'],'decomposition')
        decomposition.append(dict(comparison=f'{left} minus {right}',children_difference=a['children_capped3']-b['children_capped3'],motherhood_component=e,conditional_children_component=i))
    outputs.append(write_rows('age25_symmetric_decomposition.csv',decomposition))
    md=['# Saved selected-candidate lifecycle comparison','',
        'These are 2007 stationary approximations. The two-birth rule is experimental; neither selected case is an adopted reference. Both separately normalize completed fertility to 2.1 and pass demographic renewal. The frozen block0506 reference is used only to validate saved measurement logic and is not plotted. There is no 2023 transition result here.','',
        'The two-birth model allows one optional additional birth attempt after a first success within a four-year cell. The original permits at most one birth per cell. The candidates were recalibrated separately, so their difference does not isolate the causal effect of this rule. Earnings, entry distributions, transfers and floors, preference forms, target values and weights are retained.','',
        '## Measurement','',
        'CPS 2004/2006 women provide pooled cross-sectional children-ever-born stocks at five-year interview ages 20–44. Counts are capped at three; motherhood is the share with at least one child ever born; children among mothers is the capped mean conditional on motherhood. The model integrates four-year pre/post birth masses over matching five-year windows, with uniform within-cell timing, then aggregates mass before dividing. The exact age-25 target is a separate single-year projection, not the 25–29 bin. Completed fertility 2.1 is a normalization and is distinct from the capped 40–44 stock.','',
        'NCHS first-birth shares pool 2003–2006 period first births. Data ages 12–21 and 42–49 are folded into the endpoint model cells. The plotted shares sum to one and reproduce each saved mean first-birth age and share of first births at age 30 or later. These data do not track the CPS women.','',
        'The housing plot uses the saved age-node `lifecycle_2023.csv` from each stationary 2007 run; its filename is inherited. ACS 2007 homeownership covers household heads in a 42-metro sample. This age-profile overlay is descriptive, not a new calibration target. Rooms are uncapped in the saved model curves; ACS rooms are capped at nine and are therefore omitted. Liquid wealth is plotted in model income units; no matched empirical age curve is asserted.','',
        '## Quantitative comparison','',
        '| Item | Original selected | Two-birth selected |', '|---|---:|---:|',
        f"| Scored loss | {float(identities['Original selected']['loss']):.3f} | {float(identities['Two-birth selected']['loss']):.3f} |",
        f"| Children ever born at exact age 25 | {original['Original selected']['children_capped3']:.3f} | {original['Two-birth selected']['children_capped3']:.3f} |",
        f"| Mother share at exact age 25 | {original['Original selected']['mother_share']:.3f} | {original['Two-birth selected']['mother_share']:.3f} |",
        f"| Children among mothers at exact age 25 | {original['Original selected']['children_among_mothers_capped3']:.3f} | {original['Two-birth selected']['children_among_mothers_capped3']:.3f} |",'',
        'CPS exact-age-25 comparison: children 0.810, mother share 0.457, children among mothers 1.770. The symmetric decomposition in `age25_symmetric_decomposition.csv` is an arithmetic identity, not causal attribution. `target_costs.csv` holds all 14 target and validation rows, including weights and loss contributions; `demographic_comparison.csv` holds normalization and renewal values. `plotted_data.csv`, `first_birth_age_cells.csv`, and `other_lifecycle_data.csv` preserve all plotted values.','',
        '`weak_direction_round_centers.csv` reports the smallest right singular vector of each saved final-round 10-by-10 Jacobian in scaled and physical parameter units. Its derivatives are at round centers, not at the selected candidates; numerical rank at relative cutoff 1e-6 is scale and step dependent and does not establish statistical identification.','',
        'Source hashes, selected SUCCESS artifact hashes, reference export hashes, and numeric replay checks are in `verification.json`. No checkpoint was opened, no model module imported, and no solve run. The existing 17 standard plots per candidate were not modified.']
    p=OUT/'comparison.md';p.write_text('\n'.join(md)+'\n');outputs.append(p)
    verification={'status':'PASS','slurm_job_id':os.environ['SLURM_JOB_ID'],'model_solves':0,'model_imports':0,'checkpoint_reads':0,'selected_cases':{k:str(v) for k,v in CASES.items()},'input_sha256':inputs,'output_sha256':{p.name:sha(p) for p in outputs},'replays':['selected SUCCESS artifact hashes','reference artifact hashes','reference CPS lifecycle','reference NCHS first-birth shares','age25 scored target','first-birth mean and age30+ share','symmetric age25 decomposition'],'standard_plots_changed':False,'limitations':['Rooms model uncapped; available ACS rooms capped9','Wealth age curve model only','Stationary 2007 approximation versus pooled CPS 2004/2006 and NCHS 2003–2006','No fresh selected-point Jacobian or identification claim']}
    (OUT/'verification.json').write_text(json.dumps(verification,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','output':str(OUT),'files':len(outputs),'model_solves':0}))

if __name__=='__main__':main()
