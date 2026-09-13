"""Compare two verified, equally rebated tax paths from identical 2023 households."""
from pathlib import Path
import argparse,csv,hashlib,json,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

POLICIES=('baseline_rebate','tax2_rebate')
METRICS=('asset_price','renter_price','housing_demand','resident_persons','household_heads',
    'birth_children_topcode_adjusted','owner_rate','pension_period_units',
    'equal_transfer_period_units','period_tfr_topcode_adjusted','rooms_per_head')
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(condition,message):
    if not condition:raise ValueError(message)
def numeric(value):return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)
def relative(revenue,outlays):
    require(numeric(revenue) and numeric(outlays),'Nonfinite fiscal account')
    scale=max(abs(revenue),abs(outlays));return (revenue-outlays)/scale if scale else 0.

def check_case(case,state):
    case=Path(case).resolve();state=Path(state).resolve()
    contract=load(case/'contract_receipt.json');conditional=contract.get('conditional_history_count')
    if conditional is not None:
        require(type(conditional) is int and conditional>0 and contract.get('history_refitted') is False,
            'Conditional forecast must identify its historical horizon and absence of refit')
    completed='conditioning_history_complete.json' if conditional is not None else 'finite_history_complete.json'
    files=[case/name for name in ('contract_receipt.json','realized_fit.json',completed)]+[state]
    contract,fit,complete=map(load,files[:3]);proof=load(state)
    require([x['year'] for x in fit]==[2007,2011,2015,2019] and complete['realized']==fit,'Incomplete historical fit')
    for row in fit:
        require(all(numeric(row[k]) for k in ('psi','target','model','gap')) and
            row['gap']==row['model']-row['target'] and abs(row['gap'])<=.005,'Invalid historical fit')
    require(proof.get('status')=='passed' and proof.get('common_initial_g_pre') is True,'Identical native initial populations required')
    require(all(proof.get(k) is True for k in ('common_grid','common_supply_rule','all_other_parameter_fields_exact','worker_income_exact')),
        'Common native grid, supply rule and non-policy parameters required')
    require(proof.get('conditional_history_count')==conditional,'Conditional history label differs from native proof')
    require(proof.get('case_contract_sha256')==sha(files[0]),'State proof belongs to a different case')
    count=contract['count'];require(isinstance(count,int) and count>0,'Invalid forecast count')
    years=list(range(2023,2023+4*count,4));series={}
    for name,tax in zip(POLICIES,(.01,.02)):
        folder=case/'policies'/name;paths=[folder/x for x in ('summary.json','root_receipt.json','rows.json','fertility.json')]
        summary,root,rows,fert=map(load,paths);files+=paths
        require(proof['policy_root_sha256'].get(name)==sha(paths[1]),'Policy root differs from native state proof')
        require(summary.get('finite_converged') is True and summary.get('annual_tax')==tax,'Policy summary not accepted')
        require(root.get('converged') is True and root.get('finite_horizon_market_fiscal_converged') is True and
            root.get('final',{}).get('mapping_valid') is True,'Policy root not accepted')
        require(root.get('horizon_verified') is False and root.get('production_eligible') is False,'Expected provisional finite-horizon receipt')
        require(root.get('case')==contract['case'] and root.get('count')==count and root.get('start_year')==2023 and
            root.get('psi')==fit[-1]['psi'],'Policy case, horizon, start year or preference mismatch')
        replay=root.get('final_reproduction_max_abs');require(numeric(replay) and 0<=replay<=2e-10,'Missing or failed exact replay')
        residual=root['final'].get('residual');require(isinstance(residual,list) and len(residual)==3*(count+1) and
            all(numeric(x) and abs(x)<2e-4 for x in residual),'Final market/fiscal residual fails retained gate')
        require([x['calendar_year'] for x in rows]==years and [x['calendar_year'] for x in fert]==years,'Policy clocks or counts disagree')
        series[name]={}
        for row,f in zip(rows,fert):
            values={k:row[k] for k in METRICS[:-2]}
            values['period_tfr_topcode_adjusted']=f['period_tfr_topcode_adjusted']
            require(all(numeric(x) for x in values.values()) and values['household_heads']>0,'Missing or nonfinite policy quantities')
            require(0<=values['owner_rate']<=1,'Ownership probability outside bounds')
            values['rooms_per_head']=values['housing_demand']/values['household_heads']
            require(abs(row['relative_market_residual'])<2e-4 and
                abs(relative(row['payroll_tax_revenue'],row['pension_outlays']))<1e-6 and
                abs(relative(row['property_tax_revenue'],row['equal_transfer_outlays']))<1e-6,'Dated market or fiscal gate fails')
            series[name][row['calendar_year']]=values
    return files,years,series

def build(case,out,state):
    files,years,series=check_case(case,state);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    conditional=load(Path(case)/'contract_receipt.json').get('conditional_history_count')
    with (out/'comparison.csv').open('w',newline='') as h:
        fields=['decision_year','fertility_window_end','metric','baseline','reform','absolute_change','percent_change','percentage_point_change']
        writer=csv.DictWriter(h,fieldnames=fields);writer.writeheader()
        for year in years:
            for metric in METRICS:
                a,b=[series[name][year][metric] for name in POLICIES];gap=b-a
                writer.writerow(dict(decision_year=year,fertility_window_end=year+4,metric=metric,baseline=a,reform=b,
                    absolute_change=gap,percent_change=100*gap/a if a else '',percentage_point_change=100*gap if metric=='owner_rate' else ''))
    fig,axes=plt.subplots(2,3,figsize=(12,7))
    panels=[('period_tfr_topcode_adjusted','Period fertility',1),('rooms_per_head','Mean physical rooms per head',1),
        ('owner_rate','Homeownership (%)',100),('resident_persons','Resident persons (model units)',1),
        ('asset_price','House price per room',1),('pension_period_units','Pension per model period',1)]
    plotted={}
    for ax,(metric,title,scale) in zip(axes.flat,panels):
        plotted[metric]={name:[scale*series[name][y][metric] for y in years] for name in POLICIES}
        for name,label,color in zip(POLICIES,('1% tax, equal rebate','2% tax, equal rebate'),('#25659a','#bf463a')):
            ax.plot(years,plotted[metric][name],'o-',label=label,color=color)
        ax.set(title=title,xlabel='Decision year');ax.grid(alpha=.2)
    axes[0,0].legend(frameon=False,fontsize=8)
    title='Property tax with equal rebates and balanced pensions'
    note='Provisional finite horizon; horizon adequacy unverified. Fertility decision 2023 corresponds to the 2024–2027 flow.'
    if conditional is not None:
        note=(f'{len(years)}-period forecast conditional on a {conditional}-period fitted history; no historical refit.\n'+note)
    fig.suptitle(title);fig.text(.5,.012,note,ha='center',fontsize=8);fig.tight_layout(rect=(0,.055,1,.95))
    for ext in ('png','pdf'):fig.savefig(out/f'policy_comparison.{ext}',dpi=160)
    plt.close(fig)
    receipt=dict(status='passed',horizon_verified=False,production_eligible=False,years=years,
        source_sha256={str(p):sha(p) for p in files},state_verification=str(Path(state).resolve()),conditional_history_count=conditional)
    receipt['output_sha256']={name:sha(out/name) for name in ('comparison.csv','policy_comparison.pdf','policy_comparison.png')}
    (out/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    (out/'figure_manifest.json').write_text(json.dumps(dict(receipt,title=title,footnote=note,
        supplemental=True,plotted_series=plotted,ownership_plot_units='percent'),indent=2)+'\n')
    return receipt

def main():
    p=argparse.ArgumentParser();p.add_argument('--case-dir',required=True);p.add_argument('--out',required=True)
    p.add_argument('--state-verification',required=True);a=p.parse_args();build(a.case_dir,a.out,a.state_verification)
if __name__=='__main__':main()
