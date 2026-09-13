"""Rebuild the five retained figures from a complete verified native history."""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
ROOT=Path(__file__).resolve().parents[3]
DEFAULT_DATA=ROOT/'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/patch_readout/data'
BLUE='#1f5fa6';RED='#c73e3a'
def read(p):return json.loads(Path(p).read_text())
def csvread(p):
    with Path(p).open() as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def selected_folder(case_dir,fit,contract):
    # Map the explicit recorded window/trial suffix into a locally copied case.
    parts=Path(fit['folder']).parts;marker='window_'+str(fit['year'])
    if parts.count(marker)!=1:raise ValueError('Selected fit lacks its unique dated window')
    base=(case_dir/Path(*parts[parts.index(marker):])).resolve()
    if not base.is_relative_to(case_dir.resolve()):raise ValueError('Fit path escapes case')
    matches=[]
    for p in (base,base/'alternative'):
        if not all((p/name).is_file() for name in ('root_receipt.json','rows.json','fertility.json')):continue
        root=read(p/'root_receipt.json');rows=read(p/'rows.json');fert=read(p/'fertility.json')
        if (root.get('start_year')==fit['year'] and root.get('count')==contract['count']
            and root.get('case')==contract['case'] and root.get('converged') is True
            and root.get('finite_horizon_market_fiscal_converged') is True
            and root.get('final',{}).get('mapping_valid') is True
            and math.isfinite(float(root.get('final_reproduction_max_abs',math.inf)))
            and float(root['final_reproduction_max_abs'])<=2e-10
            and abs(float(root['psi'])-float(fit['psi']))<=1e-12
            and rows and fert and rows[0]['calendar_year']==fit['year']==fert[0]['calendar_year']
            and abs(float(fert[0]['period_tfr_topcode_adjusted'])-float(fit['model']))<=2e-10):matches.append(p)
    if len(matches)!=1:raise ValueError('Expected one matching verified primary/alternative forecast')
    return matches[0]

def verified_inputs(case_dir,readout_dir):
    realized=read(case_dir/'realized_fit.json');complete=read(case_dir/'finite_history_complete.json')
    if [r['year'] for r in realized]!=[2007,2011,2015,2019] or complete.get('realized')!=realized:
        raise ValueError('Exactly four completed ordered historical fits required')
    for r in realized:
        if not all(math.isfinite(float(r[k])) for k in ('model','target','gap','psi')) or abs(r['gap'])>.005:
            raise ValueError('Historical fit outside retained tolerance')
        if abs(r['model']-r['target']-r['gap'])>1e-12:raise ValueError('Fit arithmetic differs')
    contract=read(case_dir/'contract_receipt.json')
    folders=[selected_folder(case_dir,r,contract) for r in realized]
    check=read(readout_dir/'verification.json');model=read(readout_dir/'model_2023.json')
    if (check.get('status')!='PASS' or check.get('finite_converged') is not True
        or check.get('verification_method')!='native_saved_snapshot_aggregate_match'
        or check.get('historical_fit_status',{}).get('realized')!=realized
        or check.get('historical_fit_status',{}).get('complete') is not True
        or model.get('calendar_year')!=2023 or model.get('finite_converged') is not True):
        raise ValueError('Matching complete native2023 readout required')
    digest=sha(folders[-1]/'root_receipt.json')
    if digest!=check.get('root_receipt_sha256') or digest!=model.get('forecast_receipt_sha256'):
        raise ValueError('Final native readout/root linkage differs')
    if complete.get('horizon_verified') is not False or check.get('horizon_verified') is not False:
        raise ValueError('This packet requires explicit unverified finite horizon')
    return realized,folders,check,model

def build(case_dir,readout_dir,initial_readout=None,data_dir=DEFAULT_DATA,out=None,pdf=None):
    case_dir=Path(case_dir).resolve();readout_dir=Path(readout_dir).resolve();data=Path(data_dir)
    realized,folders,check,full=verified_inputs(case_dir,readout_dir)
    initial=2.1 if initial_readout is None else float(read(initial_readout)['fertility']['period_tfr_topcode_adjusted'])
    if not math.isfinite(initial):raise ValueError('Nonfinite initial fertility')
    final_fert=read(folders[-1]/'fertility.json')
    path=[read(p/'rows.json')[0] for p in folders]+read(folders[-1]/'rows.json')[1:]
    if [r['calendar_year'] for r in path]!=list(range(2007,path[-1]['calendar_year']+1,4)):
        raise ValueError('Expected unique consecutive dated quantity rows')
    r=next(x for x in path if x['calendar_year']==2023)
    for x in path:
        for key in ('asset_price','renter_price','housing_demand','household_heads','resident_persons'):
            if not math.isfinite(float(x[key])):raise ValueError('Nonfinite path quantity: '+key)
    out=Path(out or case_dir/'figures');out.mkdir(parents=True,exist_ok=True)
    pdf=Path(pdf or out/'e5f_final_history_figures.pdf');pdf.parent.mkdir(parents=True,exist_ok=True)
    inputs=[case_dir/n for n in ('realized_fit.json','finite_history_complete.json','contract_receipt.json')]
    inputs += [p/n for p in folders for n in ('root_receipt.json','rows.json','fertility.json')]
    inputs += [readout_dir/n for n in ('verification.json','model_2023.json')]
    inputs += [data/n for n in ('actual2023_age_housing_levels.csv','actual2023_female_recent_birth_rates.csv')]
    if initial_readout:inputs.append(Path(initial_readout))
    manifest={'sources':{str(p):sha(p) for p in inputs},'selected_folders':[str(p) for p in folders],
        'horizon_verified':False,'property_tax':'1% annual, equally rebated; PAYGO balanced',
        'initial_fertility':initial,'initial_source':str(initial_readout) if initial_readout else 'retained2.1normalization',
        'forecast_window_clock':'decision2023 is2024–2027 flow','deterministic_fixture':bool(full.get('deterministic_fixture')),'figures':{}}
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':12,'legend.fontsize':10})
    pages=PdfPages(pdf)
    def save(fig,name,values):
        footnote='SYNTHETIC FIXTURE — not model results' if manifest['deterministic_fixture'] else 'Finite horizon unverified; baseline property tax equally rebated; PAYGO balanced.'
        fig.text(.5,-.035,footnote,ha='center',fontsize=7,color='.35')
        fig.savefig(out/f'{name}.pdf',bbox_inches='tight');fig.savefig(out/f'{name}.png',dpi=150,bbox_inches='tight')
        pages.savefig(fig,bbox_inches='tight');plt.close(fig);manifest['figures'][name]=values
    ends=[2007]+[x['year']+4 for x in realized];mv=[initial]+[x['model'] for x in realized]
    future=final_fert[1:];fy=[2023]+[x['calendar_year']+4 for x in future]
    fv=[mv[-1]]+[x['period_tfr_topcode_adjusted'] for x in future]
    fig,ax=plt.subplots(figsize=(9,4.8));ax.plot(ends,mv,'o-',color=BLUE,label='Model: realized history')
    ax.plot(ends[1:],[x['target'] for x in realized],'s--',color=RED,label='Data')
    ax.plot(fy,fv,'o:',color=BLUE,label='Model: continuation');ax.axvline(2023,color='.65',lw=.8,ls=':')
    ax.set(xlabel='End of four-year fertility window',ylabel='Period fertility',title='Fertility: model, data and continuation')
    ax.legend(frameon=False);ax.grid(alpha=.16);fig.tight_layout(rect=(0,.045,1,1))
    save(fig,'historical_fertility',dict(years=ends,model=mv,target=[x['target'] for x in realized],future_years=fy,future=fv))
    years=[int(a['calendar_year']) for a in path]
    fig,axs=plt.subplots(1,3,figsize=(12,3.7))
    for ax,title,key,unit in zip(axs,['House prices','Rents','Housing quantity'],['asset_price','renter_price','housing_demand'],['Price per room','Rent per room / model period','Total physical rooms (model units)']):
        values=[float(a[key]) for a in path];ax.plot(years,values,'o-',color=BLUE);ax.axvline(2023,color='.65',lw=.8,ls=':');ax.set(title=title,xlabel='Year',ylabel=unit,xticks=years[::2]);ax.grid(alpha=.16)
    fig.tight_layout();save(fig,'prices_quantities_path',dict(years=years,series={k:[float(a[k]) for a in path] for k in ['asset_price','renter_price','housing_demand']}))
    if full:
        prof=full['profile'];t=prof['totals'];mass=t['households'];age=prof['rows'];counts=np.asarray(prof['number_children_mass'])
        fig,axs=plt.subplots(1,4,figsize=(14,3.8))
        groups=[(['Price','Rent'],[float(r['asset_price']),float(r['renter_price'])],'Prices','Per room'),(['Households','Persons'],[float(r['household_heads']),float(r['resident_persons'])],'Population','Model population units'),(['Consumption','Rooms'],[t['consumption']/mass,t['rooms']/mass],'Household means','Consumption units / physical rooms'),(['0','1','2','3+'],(100*counts/counts.sum()).tolist(),'Children ever born','Percent of households')]
        for ax,(labels,values,title,ylabel) in zip(axs,groups):
            bars=ax.bar(labels,values,color=BLUE,width=.65);ax.set(title=title,ylabel=ylabel,ylim=(0,max(values)*1.22));ax.tick_params(axis='x',labelsize=9)
            for b,v in zip(bars,values):ax.text(b.get_x()+b.get_width()/2,v+max(values)*.025,f'{v:.2f}',ha='center',fontsize=9)
        fig.tight_layout();save(fig,'equilibrium_2023',dict(groups=groups,calendar_year=2023))
        empirical=csvread(data/'actual2023_age_housing_levels.csv');byage={int(a['age_lower']):a for a in empirical}
        x=np.array([a['age']+1.5 for a in age]);d=[byage[int(a['age'])] for a in age]
        comparisons=[('Homeownership','Percent of households','owners','ownership_rate',100),('Housing size','Physical rooms, capped at 9','capped_rooms','mean_capped_rooms',1),('Children at home','Percent of households','with_children','with_minor_rate',100)]
        fig,axs=plt.subplots(1,3,figsize=(13,4));mseries={}
        for ax,(title,ylabel,mkey,dkey,scale) in zip(axs,comparisons):
            mv=[scale*a[mkey]/a['households'] for a in age];dv=[scale*float(a[dkey]) for a in d];ax.plot(x,mv,'-',lw=2,color=BLUE,label='Model');ax.plot(x,dv,'--',lw=2,color=RED,label='ACS 2023');ax.set(title=title,xlabel='Age of household head',ylabel=ylabel,xticks=[20,35,50,65,80]);ax.grid(alpha=.16);mseries[mkey]=dict(model=mv,data=dv)
        axs[0].legend(frameon=False);fig.tight_layout(rect=(0,.065,1,1));fig.text(.5,.015,'Children at home: model dependents; ACS resident own children under 18.',ha='center',fontsize=9);save(fig,'lifecycle_2023',dict(ages=x.tolist(),series=mseries))
        f=full['fertility'];obs=csvread(data/'actual2023_female_recent_birth_rates.csv')
        fy=np.array(f['age_cell_start']);q=fy<50;mx=fy[q]+2;my=250*np.array(f['age_specific_birth_rate_topcode_adjusted'])[q]
        dx=[int(a['age_lower'])+1.5 for a in obs];dy=[1000*float(a['recent_birth_rate']) for a in obs]
        fig,ax=plt.subplots(figsize=(9,4.6));ax.plot(mx,my,'o-',color=BLUE,label='Model: households');ax.plot(dx,dy,'s--',color=RED,label='ACS 2023: women');ax.set(xlabel='Age',ylabel='Annual births / recent-birth reports per 1,000',title='Fertility by age',xticks=[20,25,30,35,40,45,50]);ax.legend(frameon=False);ax.grid(alpha=.16)
        fig.text(.5,.015,'Model: annualized 2024–2027 birth flow. ACS: women reporting a birth in the previous 12 months.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.055,1,1));save(fig,'fertility_age_2023',dict(model_age=mx.tolist(),model=my.tolist(),data_age=dx,data=dy,measurement='Model household denominator versus ACS women; model four-year flow annualized; descriptive comparison'))
    pages.close()
    (out/'figure_manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    return manifest

def main():
    p=argparse.ArgumentParser();p.add_argument('--case-dir',type=Path,required=True);p.add_argument('--readout-dir',type=Path,required=True)
    p.add_argument('--initial-readout',type=Path);p.add_argument('--data-dir',type=Path,default=DEFAULT_DATA)
    p.add_argument('--out',type=Path);p.add_argument('--pdf',type=Path);a=p.parse_args()
    build(a.case_dir,a.readout_dir,a.initial_readout,a.data_dir,a.out,a.pdf)
if __name__=='__main__':main()
