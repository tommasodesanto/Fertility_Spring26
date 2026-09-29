"""Torch-only readout of saved zero-cost cases; no model or checkpoint loads."""
from __future__ import annotations
import csv, hashlib, importlib.util, json, math, os, sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ORIGINAL=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
REMOTE=Path('/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project')
BASE=HERE.parent
NIGHT=BASE/'two_stream_overnight_v1'
PRIOR=NIGHT/'comparison_v1/build_comparison.py'
PRIOR_SHA256='00dc3f5c31bbb2fe3f25ebad3f017f3cd30246b34231332f8f337e29dbdea957'
RUN=HERE/'run_v1/zero_cost'
OUT=HERE/'readout_v1'
ANCHOR=NIGHT/'run_v1/one_birth/one_birth_024_gn1_0/case'
AUDIT=BASE/'measurement_audit_v1'
ACS=NIGHT/'comparison_v1/input/actual2007_age_housing_levels.csv'
NEEDED=('observers.json','lifecycle_2023.csv','receipt.json','target_fit.csv','parameters.csv')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()

def read(path):return json.loads(Path(path).read_text())

def staged(path):
    path=Path(path)
    return REMOTE/path.relative_to(ORIGINAL) if path.is_relative_to(ORIGINAL) else path

def csv_rows(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))

def write_csv(name,records):
    assert records,name
    path=OUT/name
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');w.writeheader();w.writerows(records)
    return path

def load_prior():
    assert sha(PRIOR)==PRIOR_SHA256,'Validated lifecycle projection source changed'
    spec=importlib.util.spec_from_file_location('saved_lifecycle_projection',PRIOR)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    return module

def case_from_record(record):
    assert record['status']=='success'
    success=record['result'];path=staged(success['case_path'])
    assert path.is_relative_to(RUN) and path.name=='case' and path.parent.name==record['candidate_id']
    assert sha(staged(success['success_path']))==success['success_sha256']
    return path

def main():
    assert sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm required'
    prior=load_prior();records=read(RUN/'records.json')
    centers=[r for r in records if r['role']=='zero_cost_center' and r['status']=='success']
    assert len(centers)==1,'Authenticated zero-cost center unavailable; report first-stage censor separately'
    OUT.mkdir(parents=True,exist_ok=True)
    final=read(RUN/'FINAL.json') if (RUN/'FINAL.json').exists() else None
    certified=bool(final and final['numerical_repeat_screens_passed'])
    selected=final['selected'] if final and final.get('selected') else None
    cases={'Original selected':ANCHOR,'Zero cost, fixed coordinates':case_from_record(centers[0])}
    if selected and selected['candidate_id']!=centers[0]['candidate_id']:
        cases['Restricted selected'+(' (repeat checked)' if certified else ' (provisional)')]=staged(selected['case_path'])
    elif selected:
        # The center itself can be the selected restricted case.
        cases['Restricted selected'+(' (repeat checked)' if certified else ' (provisional)')]=staged(selected['case_path'])
    outputs=[];inputs={str(PRIOR):sha(PRIOR),str(Path(__file__)):sha(__file__)}
    data={}
    for label,case in cases.items():
        success=read(case.parent/'SUCCESS.json')
        assert success['status']=='passed' and staged(success['case_path'])==case.resolve()
        hashes={staged(a['path']).relative_to(case).as_posix():a['sha256'] for a in success['artifacts']}
        for name in NEEDED:
            path=case/name;assert sha(path)==hashes[name],(label,name)
            inputs[str(path)]=hashes[name]
        fit=csv_rows(case/'target_fit.csv');params=csv_rows(case/'parameters.csv')
        assert len(fit)==14 and len(params)==31
        receipt=read(case/'receipt.json')
        assert abs(receipt['normalization']['completed_fertility']-2.1)<=5e-4
        assert abs(receipt['adult_entry_gate']['fertility_gap'])<=5e-4
        assert abs(sum(float(r['loss_contribution'] or 0) for r in fit)-receipt['loss'])<=1e-8
        cost=next(r for r in params if r['parameter']=='first_birth_fixed_cost')
        expected=0.35270914196085973 if label=='Original selected' else 0.0
        assert float(cost['estimate'])==expected and receipt['point']['first_birth_fixed_cost']==expected
        observer=read(case/'observers.json')['fertility']['uniform_birth_time']['accounting']
        profile=csv_rows(case/'lifecycle_2023.csv');assert len(profile)==17
        data[label]=dict(fit=fit,params=params,receipt=receipt,observer=observer,profile=profile)
    labels=list(cases)
    fit_rows=[]
    for i,base in enumerate(data[labels[0]]['fit']):
        contract=('moment','role','target','weight')
        row={k:base[k] for k in contract}
        for label in labels:
            item=data[label]['fit'][i];assert all(item[k]==base[k] for k in contract)
            row[label+' model']=item['model'];row[label+' gap']=item['gap'];row[label+' loss']=item['loss_contribution']
        fit_rows.append(row)
    outputs.append(write_csv('full_target_fit.csv',fit_rows))
    parameter_rows=[]
    for i,base in enumerate(data[labels[0]]['params']):
        row={'parameter':base['parameter'],'original_lower':base['lower'],'original_upper':base['upper']}
        for label in labels:
            item=data[label]['params'][i];assert item['parameter']==base['parameter']
            row[label+' estimate']=item['estimate'];row[label+' lower']=item['lower'];row[label+' upper']=item['upper']
            row[label+' near_bound']=item['near_bound'];row[label+' status']=item['status']
        parameter_rows.append(row)
    outputs.append(write_csv('full_parameters.csv',parameter_rows))
    empirical=csv_rows(AUDIT/'fertility_lifecycle_matched_windows.csv')
    first=csv_rows(AUDIT/'first_birth_age_cells.csv')
    acs=csv_rows(ACS)
    assert sha(ACS)=='9b4433299f43c1c23b8aa28ba224a41ba6ee21d1a7a2265a336f90c831ecb8b9'
    for p in (AUDIT/'fertility_lifecycle_matched_windows.csv',AUDIT/'first_birth_age_cells.csv',ACS):inputs[str(p)]=sha(p)
    ages=[(int(r['age_lower']),int(r['age_upper'])) for r in empirical]
    assert ages==[(20,24),(25,29),(30,34),(35,39),(40,44)]
    lifecycle=[];first_rows=[];other=[]
    fmap={'children_capped3':'capped3','mother_share':'mother_share','children_among_mothers_capped3':'given_mother'}
    for row in empirical:
        lo=int(row['age_lower']);hi=int(row['age_upper'])
        for label in labels:
            vals=prior.project(data[label]['observer'],lo,hi)
            for metric in fmap:lifecycle.append(dict(age_lower=lo,age_upper=hi,metric=metric,series=label,value=vals[metric]))
        for metric,source in fmap.items():
            lifecycle.append(dict(age_lower=lo,age_upper=hi,metric=metric,series='CPS 2004/2006',value=float(row['data_'+source])))
    for label in labels:
        vals=prior.project(data[label]['observer'],25,25)
        early=next(r for r in data[label]['fit'] if r['moment']=='early_fertility')
        assert abs(vals['children_capped3']-float(early['model']))<=1e-10
        lifecycle.append(dict(age_lower=25,age_upper=25,metric='exact_age25_children',series=label,value=vals['children_capped3']))
    for i,row in enumerate(first):
        lo=int(row['age_lower']);hi=int(row['age_upper']);assert (lo,hi)==(18+4*i,21+4*i)
        first_rows.append(dict(age_lower=lo,age_upper=hi,series='NCHS 2003–2006',share=float(row['empirical_share'])))
        for label in labels:
            a=data[label]['observer'];assert a['age_cell_start'][i]==lo
            first_rows.append(dict(age_lower=lo,age_upper=hi,series=label,share=a['parity_birth_flows_by_age'][i][0]/a['first_birth_flow']))
    for label in labels:
        for r in data[label]['profile']:
            lo=int(float(r['age_node']))
            for metric,value in [('homeownership',float(r['owner_rate'])),('rooms_uncapped',float(r['mean_rooms'])),('liquid_wealth_model_units',float(r['mean_liquid_wealth']))]:
                assert math.isfinite(value)
                other.append(dict(age_lower=lo,age_upper=lo+3,metric=metric,series=label,value=value))
    for r in acs:other.append(dict(age_lower=int(r['age_lower']),age_upper=int(r['age_lower'])+3,metric='homeownership',series='ACS 2007',value=float(r['ownership_rate'])))
    outputs.extend([write_csv('fertility_lifecycle.csv',lifecycle),write_csv('first_birth_age_cells.csv',first_rows),write_csv('other_lifecycle.csv',other)])
    colors={'Original selected':'#17678a','Zero cost, fixed coordinates':'#b56332','CPS 2004/2006':'#555555','NCHS 2003–2006':'#555555','ACS 2007':'#555555'}
    for label in labels[2:]:colors[label]='#368963'
    plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12.8,8.5));axes=axes.flat
    specs=[('children_capped3','Children ever born per woman, capped at 3'),('mother_share','Women who are mothers'),('children_among_mothers_capped3','Children among mothers, capped at 3')]
    for ax,(metric,title) in zip(axes[:3],specs):
        for label in ['CPS 2004/2006',*labels]:
            rows=[r for r in lifecycle if r['metric']==metric and r['series']==label]
            ax.plot(range(len(rows)),[r['value'] for r in rows],marker='s' if label.startswith('CPS') else 'o',color=colors[label],label=label)
        ax.set(title=title,xticks=range(5),xticklabels=[f'{lo}–{hi}' for lo,hi in ages],xlabel='Age at interview');ax.grid(axis='y',alpha=.2)
    ax=axes[3]
    for label in ['NCHS 2003–2006',*labels]:
        rows=[r for r in first_rows if r['series']==label]
        ax.plot(range(len(rows)),[100*r['share'] for r in rows],marker='s' if label.startswith('NCHS') else 'o',color=colors[label],label=label)
    ax.set(title='Share of first births by age cell',xticks=range(7),xticklabels=[f'{r["age_lower"]}–{r["age_upper"]}' for r in first_rows if r['series']=='NCHS 2003–2006'],ylabel='Percent');ax.grid(axis='y',alpha=.2)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=4,frameon=False,bbox_to_anchor=(.5,.055))
    fig.text(.5,.012,'2007 stationary approximation; CPS 2004/2006 cross-sections; NCHS 2003–2006 first births. Exact age-25 target is separate from ages 25–29.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.1,1,.97));p=OUT/'fertility_lifecycle.png';fig.savefig(p,dpi=165);plt.close(fig);outputs.append(p)
    fig,axes=plt.subplots(1,3,figsize=(14.2,4.7))
    for ax,(metric,title,ylabel) in zip(axes,[('homeownership','Homeownership','Share of heads'),('rooms_uncapped','Rooms, model only','Uncapped rooms'),('liquid_wealth_model_units','Liquid wealth, model only','Model income units')]):
        series=['ACS 2007',*labels] if metric=='homeownership' else labels
        for label in series:
            rows=[r for r in other if r['metric']==metric and r['series']==label]
            ax.plot([r['age_lower']+1.5 for r in rows],[r['value'] for r in rows],marker='o',ms=3,color=colors[label],label=label)
        ax.set(title=title,xlabel='Age',ylabel=ylabel,xticks=[20,35,50,65,80]);ax.grid(axis='y',alpha=.2)
    axes[0].legend(frameon=False,fontsize=8)
    fig.text(.5,.012,'ACS 2007 ownership: household heads in 42 metros. ACS rooms capped at 9; saved model rooms uncapped. No matched empirical liquid-wealth age profile.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.08,1,.97));p=OUT/'other_lifecycle.png';fig.savefig(p,dpi=165);plt.close(fig);outputs.append(p)
    report=dict(status='repeat_checked' if certified else 'provisional',model_imports=0,model_solves=0,checkpoint_reads=0,
        prior_projection_source=str(PRIOR),prior_projection_sha256=inputs[str(PRIOR)],cases={k:str(v) for k,v in cases.items()},
        selected_repeat_count=final['repeat_count'] if final else 0,limitations=['No global infeasibility conclusion','No selected-point identification claim','Rooms and wealth curves model only'],
        input_sha256=inputs,output_sha256={p.name:sha(p) for p in outputs},standard_17_plots_changed=False)
    (OUT/'verification.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':report['status'],'output':str(OUT),'cases':list(cases),'model_solves':0}))

if __name__=='__main__':main()
