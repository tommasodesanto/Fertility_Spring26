"""Rebuild presentation figures from the saved2019stationary-start patch and ACS2023."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/patch_readout'
BLUE='#1f5fa6';RED='#c73e3a'
def read(p):return json.loads(Path(p).read_text())
def csvread(p):return list(csv.DictReader(Path(p).open()))
def main():
    source=BASE/'source';data=BASE/'data';out=BASE/'figures';out.mkdir(exist_ok=True)
    path=csvread(source/'expected_transition.csv');fert=read(source/'fertility.json');static=read(source/'stationary_history.json');initial=read(source/'initial.json')
    r=next(r for r in path if int(r['calendar_year'])==2023)
    root=read(source/'root_receipt.json');assert root['finite_horizon_market_fiscal_converged']
    assert abs(float(path[0]['psi_child'])-.09239514522037684)<1e-14
    full=read(source/'model_2023.json') if (source/'model_2023.json').exists() else None
    if full:
        check=read(source/'verification.json');assert check['status']=='PASS' and full['calendar_year']==2023
        assert full['forecast_receipt_sha256']==hashlib.sha256((source/'root_receipt.json').read_bytes()).hexdigest()
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':12,'legend.fontsize':10})
    pdf=ROOT/'output/pdf/e5f_patch_review.pdf';pages=PdfPages(pdf);manifest={}
    def save(fig,name,values):
        fig.savefig(out/f'{name}.pdf',bbox_inches='tight');fig.savefig(out/f'{name}.png',dpi=150,bbox_inches='tight');pages.savefig(fig,bbox_inches='tight');plt.close(fig);manifest[name]=values
    ends=np.array([2007,2011,2015,2019,2023]);model=[initial['fertility']['period_tfr_topcode_adjusted'],*[s['model'] for s in static],fert[0]['period_tfr_topcode_adjusted']];target=[s['data'] for s in static]+[1.64575]
    fig,ax=plt.subplots(figsize=(9,4.8));ax.plot(ends,model,'o-',color=BLUE,label='Model');ax.plot(ends[1:],target,'s--',color=RED,label='Data')
    ax.axvline(2019,color='.65',lw=.8,ls=':');ax.set(xticks=ends,xlabel='End of four-year fertility window',ylabel='Period fertility',ylim=(1.45,2.2),title='Fertility: model and data');ax.legend(frameon=False);ax.grid(alpha=.16)
    fig.text(.5,.01,'Stationary approximation through 2019; transition thereafter. Initial model fertility normalized to 2.1.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.055,1,1));save(fig,'historical_fertility',dict(years=ends.tolist(),model=model,data=target))
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
        fig.text(.5,.015,'Model: annualized 2023-2027 birth flow. ACS: women reporting a birth in the previous 12 months.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.055,1,1));save(fig,'fertility_age_2023',dict(model_age=mx.tolist(),model=my.tolist(),data_age=dx,data=dy,measurement='Model household denominator versus ACS women; model four-year flow annualized; descriptive comparison'))
        valid=[a for a in age if a['age']>=22];den=sum(a['households'] for a in valid);dvalid=[a for a in empirical if int(a['age_lower'])>=22];dden=sum(float(a['hhwt']) for a in dvalid)
        table=[]
        for label,mkey,dkey,scale in [('Ownership (%)','owners','owners',100),('Mean rooms (capped at 9)','capped_rooms','capped_rooms',1),('With children at home (%)','with_children','with_minor',100)]:
            mv=scale*sum(a[mkey] for a in valid)/den;dv=scale*sum(float(a[dkey]) for a in dvalid)/dden;table.append(dict(moment=label,data=dv,model=mv,gap=mv-dv,status='2023 validation; not targeted in historical shock fitting'))
        with (out/'validation_2023.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=table[0]);w.writeheader();w.writerows(table)
        fig,ax=plt.subplots(figsize=(9,3));ax.axis('off');tab=ax.table(cellText=[[a['moment'],f"{a['data']:.2f}",f"{a['model']:.2f}"] for a in table],colLabels=['2023 outcome','Data','Model'],loc='center',cellLoc='left',colWidths=[.64,.18,.18]);tab.auto_set_font_size(False);tab.set_fontsize(12);tab.scale(1,2)
        for (row,col),cell in tab.get_celld().items():
            cell.visible_edges='TB' if row==0 else ('B' if row==len(table) else '')
            cell.set_linewidth(.8)
        ax.set_title('2023 validation',pad=10);save(fig,'validation_2023',table)
        six=csvread(BASE.parent/'age_housing_allocation/comparison_2023/large_owner_data.csv');masses=np.zeros((3,2))
        for a in prof['large_owner_age_cells']:
            for j,(lo,hi) in enumerate([(22,40),(40,60),(60,86)]):
                frac=max(0,min(a['age']+a['age_width'],hi)-max(a['age'],lo))/a['age_width'];masses[j]+=[frac*a['without_children'],frac*a['with_children']]
        mv=100*masses.ravel()/masses.sum();dv=100*np.array([float(a['share']) for a in six]);assert abs(dv.sum()-100)<1e-8
        fig,ax=plt.subplots(figsize=(11,4.5));x=np.arange(6);w=.36
        for offset,vals,color,label in [(-w/2,mv,BLUE,'Model'),(w/2,dv,RED,'ACS 2023')]:
            bars=ax.bar(x+offset,vals,w,color=color,label=label)
            for b,v in zip(bars,vals):ax.text(b.get_x()+w/2,v+.5,f'{v:.1f}',ha='center',fontsize=9)
        ax.set(xticks=x,xticklabels=[f'{a}\n{b}' for a in ['Young (22-39)','Middle (40-59)','Old (60-85)'] for b in ['No children','Children']],ylabel='Share of large owner-occupied homes (%)',ylim=(0,max(mv.max(),dv.max())*1.2));ax.legend(frameon=False);ax.tick_params(axis='x',labelsize=9);fig.tight_layout();save(fig,'intergenerational_allocation_2023',dict(model=mv.tolist(),data=dv.tolist(),rooms_threshold=6))
    pages.close();manifest.update(source='Conditional stationary households through2019, then saved2019forecast',horizon_verified=False,historical_fit_complete=False,initial_checkpoint_sha256=initial['checkpoint_sha256']);(out/'figure_verification.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(dict(pdf=str(pdf),figures=list(manifest)[:6],model2023_available=bool(full))))
if __name__=='__main__':main()
