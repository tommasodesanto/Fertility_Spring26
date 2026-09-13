"""Rebuild patch or carried-history figures and dated empirical comparisons."""
from pathlib import Path
import csv,hashlib,json
import argparse
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
def sequence_fertility(base, sequence_base):
    """Show the carried history separately from the stationary-start patch."""
    source=sequence_base/'source';out=sequence_base/'figures';out.mkdir(parents=True,exist_ok=True)
    prefix=read(source/'realized_fit.json');trial=read(source/'final_window_trial_00.json')
    receipt=read(source/'root_receipt.json');initial=read(base/'source/initial.json')
    assert len(prefix)==3 and [x['year'] for x in prefix]==[2007,2011,2015]
    assert receipt['finite_horizon_market_fiscal_converged'] and receipt['start_year']==trial['year']==2019
    assert receipt['psi']==trial['psi']
    rows=prefix+[trial];ends=[2007]+[x['year']+4 for x in rows]
    model=[initial['fertility']['period_tfr_topcode_adjusted']]+[x['model'] for x in rows]
    data=[x['data'] for x in rows]
    future=read(source/'fertility.json')[1:]
    future_years=[2023]+[x['calendar_year']+4 for x in future]
    future_values=[model[-1]]+[x['period_tfr_topcode_adjusted'] for x in future]
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':12,'legend.fontsize':10})
    fig,ax=plt.subplots(figsize=(9,4.8))
    ax.plot(ends,model,'o-',color=BLUE,label='Model: historical fit')
    ax.plot(future_years,future_values,'o:',color=BLUE,label='Model: constant-preference continuation')
    ax.plot(ends[1:],data,'s--',color=RED,label='Data')
    ax.axvline(2023,color='.6',lw=.8,ls='--');ax.axvspan(2023,future_years[-1],color=BLUE,alpha=.035)
    ax.set(xticks=ends+future_years[1::2],xlabel='End of four-year fertility window',ylabel='Period fertility',ylim=(1.45,2.2),
        title='Fertility: historical fit and continuation')
    ax.legend(frameon=False,loc='lower left');ax.grid(alpha=.16)
    fig.text(.5,.01,'No property-tax rebate. Six-date forecasts; terminal-horizon sensitivity has not been verified.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=(0,.055,1,1))
    for ext in ('png','pdf'):fig.savefig(out/f'historical_fertility.{ext}',dpi=150,bbox_inches='tight')
    plt.close(fig)
    (out/'figure_verification.json').write_text(json.dumps(dict(years=ends,model=model,data=data,future_years=future_years,future_model=future_values,
        root_receipt_sha256=hashlib.sha256((source/'root_receipt.json').read_bytes()).hexdigest(),
        initial_checkpoint_sha256=initial['checkpoint_sha256'],property_tax_rebated=False,
        horizon_verified=False,final_window_fit_accepted=False),indent=2)+'\n')
    print(str(out/'historical_fertility.png'))
def sequence_validation(sequence_base):
    """Read all original moment families at 2023, preserving empirical vintages."""
    source=sequence_base/'source/readout_2023';out=sequence_base/'figures'
    model=read(source/'model_2023.json');check=read(source/'verification.json')
    snapshot_gap=check.get('snapshot_maximum_abs')
    native_verified=(check.get('verification_method')=='native_saved_snapshot_aggregate_match'
        and check.get('replay_performed') is False and check.get('finite_converged') is True
        and isinstance(snapshot_gap,(int,float)) and np.isfinite(snapshot_gap)
        and 0<=snapshot_gap<=2e-10)
    assert check['status']=='PASS' and (check.get('replay_maximum_abs')==0 or native_verified) and model['calendar_year']==2023
    assert not read(source/'measurement_verification.json')['errors']
    assert model['forecast_receipt_sha256']==hashlib.sha256((source/'root_receipt.json').read_bytes()).hexdigest()
    empirical=ROOT/'output/model/e5f_matched_pf_20260909a/design_research'
    paths={
        'cps':ROOT/'code/data/cps_fertility/output/cps_fertility_targets.csv',
        'nchs':ROOT/'code/data/nchs_natality_timing/first_birth_counts_year_age.csv',
        'acs':empirical/'housing/early_housing_target_candidates.csv',
        'rooms':ROOT/'code/data/psid_followup_mar2026/output/sa_rooms_first_birth_household_aligned_v1/target_receipt.csv',
        'wealth':empirical/'wealth/aggregate_wealth_results.csv',
        'old':empirical/'wealth/old_wealth_results.csv',
        'bequest':ROOT/'output/model/e5f_matched_pf_20260909a/initial_calibration_contract/working_weights.csv',
    }
    cps={r['moment_key']:float(r['estimate']) for r in csvread(paths['cps'])}
    counts=[(int(r['age']),float(r['n_first_births'])) for r in csvread(paths['nchs']) if int(r['year'])==2023]
    assert counts and min(a for a,n in counts)==12 and max(a for a,n in counts)==49
    midpoint=lambda a:20 if a<22 else 44 if a>=42 else 20+4*((a-18)//4)
    total=sum(n for a,n in counts)
    nchs_mean=sum(n*midpoint(a) for a,n in counts)/total
    nchs_share=sum(n for a,n in counts if a>=30)/total
    acs={r['moment']:float(r['point']) for r in csvread(paths['acs']) if r['window']=='2023'}
    f=model['fertility_stock_timing'];h=model['housing_wealth']['moments'];flow=model['fertility']
    third=np.asarray(flow['birth_flow_third_bin_entry']);explicit=np.asarray(flow['birth_flow_explicit']);adjusted=np.asarray(flow['birth_flow_topcode_adjusted'])
    top=3+(adjusted-explicit)[third>1e-12]/third[third>1e-12]
    assert np.allclose(top,top[0],rtol=0,atol=1e-12)
    shares=f['parity_shares_40_44'];completed=shares['1']+2*shares['2']+float(top[0])*shares['3plus']
    rows=[]
    def add(key,label,d,m,vintage,source_key,scale=1,decimals=2,note=''):
        assert np.isfinite(d) and np.isfinite(m)
        rows.append(dict(moment_key=key,moment=label,data=scale*d,model=scale*m,gap=scale*(m-d),
            data_vintage=vintage,model_year=2023,decimals=decimals,data_source=str(paths[source_key]),
            empirical_status='External benchmark' if key=='bequest_wealth' else 'Untargeted validation',
            weight='',loss_contribution='',measurement_note=note))
    add('completed_fertility','Completed fertility, ages 40–44',cps['tfr'],completed,'CPS 2024','cps',note='Stock of children ever born at ages 40–44; not the period fertility rate or the initial 2.1 normalization. Model top-bin representative recovered from saved birth accounting.')
    add('childlessness','Childless, ages 40–44 (%)',cps['childless_rate'],f['moments']['childless_rate_40_44'],'CPS 2024','cps',100,1)
    add('exactly_one','Exactly one child among mothers, 40–44 (%)',cps['parity_share_1']/(1-cps['childless_rate']),f['moments']['exactly_one_among_mothers_40_44'],'CPS 2024','cps',100,1)
    add('first_birth_age','Mean age at first birth (years)',nchs_mean,f['moments']['period_mean_age_first_birth'],'NCHS 2023','nchs',note='Count-weighted model-cell midpoints; model 2023–2027 birth flow versus annual 2023 births.')
    add('first_birth_share30','First births at age 30+ (%)',nchs_share,f['moments']['period_share_first_births_age30plus'],'NCHS 2023','nchs',100,1,note='Model 2023–2027 first-birth flow versus annual 2023 births.')
    add('mean_rooms','Mean occupied rooms (capped at 9)',acs['aggregate_mean_occupied_rooms_capped9_18_85'],h['aggregate_mean_occupied_rooms_capped9_18_85'],'ACS 2023','acs')
    add('ownership_30_55','Ownership, heads 30–55 (%)',acs['own_rate_30_55'],h['own_rate_30_55'],'ACS 2023','acs',100,1)
    add('first_birth_rooms','First-birth room response, −1 to +3',float(csvread(paths['rooms'])[0]['estimate']),model['dated_first_birth_rooms']['housing_response'],'PSID pooled','rooms',1,3,note='Sun–Abraham empirical contrast; model matched branch from 2019 into 2023, not a stationary 2023 comparison.')
    add('family_rooms','Rooms: 3+ versus 1–2 resident children',acs['prime30_55_resident_3plus_minus_1to2_rooms_capped9'],h['prime30_55_model_dependent_3plus_minus_1to2_rooms_capped9'],'ACS 2023','acs',1,3,note='Model dependent counts proxy resident own children under age 18.')
    add('recent_parent','Recent-parent ownership gap (pp)',acs['recent_parent_minus_no_resident_child_ownership_30_55'],model['recent_parent']['model_value'],'ACS 2023','acs',100,1,note='Model current birth into empty dependent home versus currently empty home; same flow proxy as initial calibration, not exact ACS oldest-child-age reconstruction.')
    wealth=next(float(r['estimate']) for r in csvread(paths['wealth']) if r['window']=='pooled_2005_2019')
    add('wealth_earnings','Wealth / annual gross earnings',wealth,h['aggregate_wealth_to_annual_gross_labor_earnings'],'PSID 2005–19','wealth')
    bequest=next(float(r['target']) for r in csvread(paths['bequest']) if r['restriction_id']=='bequest_wealth')
    add('bequest_wealth','Annual bequests / wealth (%)',bequest,h['annual_bequest_flow_to_aggregate_wealth'],'External benchmark','bequest',100,2,note='External historical restriction, not a 2023 empirical observation.')
    old=next(float(r['estimate']) for r in csvread(paths['old']) if r['window']=='pooled_1984_2019' and r['moment']=='old_p90_p50')
    add('old_dispersion','Wealth/income p90 / median, ages 76–84',old,h['old_total_wealth_to_annual_income_p90_p50_7684'],'PSID 1984–2019','old',note='Beginning-period wealth among living households; model pension-income proxy versus empirical family income.')
    assert len(rows)==13
    with (out/'validation_2023.csv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    def fmt(r,key):return f"{r[key]:.{r['decimals']}f}"
    # Three-column slide table; row markers retain the data-vintage mapping.
    marks=['a','a','a','b','b','c','c','d','c','c','e','f','g']
    latex=['\\begin{tabularx}{\\textwidth}{@{}Xrr@{}}','\\toprule','Moment & Data & Model \\\\','\\midrule']
    for i,(r,mark) in enumerate(zip(rows,marks)):
        label=r['moment'].replace('–','--').replace('−','$-$').replace('%',r'\%')
        latex.append(label+r'\textsuperscript{'+mark+'} & '+fmt(r,'data')+' & '+fmt(r,'model')+r' \\')
        if i in (4,9):latex.append(r'\addlinespace[4pt]')
    latex += [r'\bottomrule',r'\end{tabularx}']
    (out/'validation_2023.tex').write_text('\n'.join(latex)+'\n')
    fig,ax=plt.subplots(figsize=(11.7,7.3));ax.axis('off')
    cells=[[r['moment'],fmt(r,'data'),fmt(r,'model'),r['data_vintage']] for r in rows]
    tab=ax.table(cellText=cells,colLabels=['Moment','Data','Model 2023','Data vintage'],loc='upper center',cellLoc='left',colWidths=[.56,.10,.12,.22],bbox=[0,.15,1,.80])
    tab.auto_set_font_size(False);tab.set_fontsize(10)
    for (i,j),cell in tab.get_celld().items():
        cell.visible_edges='B' if i in (0,13) else '';cell.set_linewidth(.7)
        if i==0:cell.set_text_props(weight='bold')
        if j in (1,2):cell.set_text_props(ha='right')
    ax.set_title('2023 model: untargeted moment comparison',fontsize=15,pad=12)
    fig.text(.06,.10,'All 13 moment families from the initial calibration. CPS uses the nearest available fertility supplement; PSID moments are pooled.\nCompleted fertility is children ever born at ages 40–44, distinct from period fertility. Birth timing uses model age-cell midpoints.\nModel first-birth rooms follow the 2019–2023 matched branch. Housing samples retain the original 42-metro definition.',fontsize=9,linespacing=1.6)
    fig.text(.06,.025,'Provisional carried-history model; no property-tax rebate. Finite-path checks pass; terminal horizon is not verified.',fontsize=9,color='.3')
    fig.subplots_adjust(top=.91,bottom=.13,left=.055,right=.96)
    for ext in ('pdf','png'):fig.savefig(out/f'validation_2023.{ext}',dpi=150)
    plt.close(fig)
    manifest=dict(rows=13,model_source=str(source/'model_2023.json'),top_bin_representative=float(top[0]),
        model_source_sha256=hashlib.sha256((source/'model_2023.json').read_bytes()).hexdigest(),
        sources={k:dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for k,p in paths.items()},
        weight_or_target_changes=False,horizon_verified=False,replay_check=check)
    (out/'validation_2023_verification.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps([dict(moment=r['moment'],data=r['data'],model=r['model']) for r in rows],indent=2))

def permanent_profile_figures(profile_path, figures_dir):
    """Refresh the two 2023 model-vs-data figures from a permanent-shock profile."""
    profile_path=Path(profile_path)
    model_path=profile_path/'model_2023.json' if profile_path.is_dir() else profile_path
    source=model_path.parent
    verification_path=source/'verification.json'
    model=read(model_path);verification=read(verification_path)
    assert model['calendar_year']==2023
    assert verification['status']=='PASS'
    assert verification.get('row_reproduction',{}).get('passed') is True
    assert verification.get('observed_year')==2023 or verification.get('observed_year2023') is True
    assert verification.get('model_source_sha256')==hashlib.sha256(model_path.read_bytes()).hexdigest()
    assert float(verification['row_reproduction'].get('numeric_max_abs_gap',0))<=2e-10
    prof=model['profile'];age=prof['rows']
    empirical_path=BASE/'data/actual2023_age_housing_levels.csv'
    empirical=csvread(empirical_path);byage={int(r['age_lower']):r for r in empirical}
    assert age and all(int(r['age']) in byage for r in age)
    out=Path(figures_dir);out.mkdir(parents=True,exist_ok=True)
    model_hash=hashlib.sha256(model_path.read_bytes()).hexdigest()
    verification_hash=hashlib.sha256(verification_path.read_bytes()).hexdigest()
    empirical_hash=hashlib.sha256(empirical_path.read_bytes()).hexdigest()
    six_path=BASE.parent/'age_housing_allocation/comparison_2023/large_owner_data.csv'
    six=csvread(six_path);six_hash=hashlib.sha256(six_path.read_bytes()).hexdigest()
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.titlesize':12,'legend.fontsize':10})
    qa={}
    def write_qa(name, plotted, sources):
        payload=dict(figure=name,calendar_year=2023,model_source=str(model_path),
                     model_source_sha256=model_hash,verification_source=str(verification_path),
                     verification_source_sha256=verification_hash,source_hashes=sources,
                     plotted_arrays=plotted,arrays_match_source=True,
                     note='Permanent preference shock; transition not converged.')
        (out/f'{name}_verification.json').write_text(json.dumps(payload,indent=2)+'\n')
        qa[name]=payload
    # Canonical Lifecycle Fit in 2023: homeownership, capped rooms, and with dependents.
    x=np.array([r['age']+1.5 for r in age],dtype=float)
    d=[byage[int(r['age'])] for r in age]
    comparisons=[('Homeownership','Percent of households','owners','ownership_rate',100),
                 ('Housing size','Physical rooms, capped at 9','capped_rooms','mean_capped_rooms',1),
                 ('Children at home','Percent of households','with_children','with_minor_rate',100)]
    fig,axs=plt.subplots(1,3,figsize=(13,4));mseries={}
    for ax,(title,ylabel,mkey,dkey,scale) in zip(axs,comparisons):
        mv=np.asarray([scale*r[mkey]/r['households'] for r in age],dtype=float)
        dv=np.asarray([scale*float(r[dkey]) for r in d],dtype=float)
        line_m,=ax.plot(x,mv,'-',lw=2,color=BLUE,label='Model')
        line_d,=ax.plot(x,dv,'--',lw=2,color=RED,label='ACS 2023')
        assert np.array_equal(np.asarray(line_m.get_ydata()),mv)
        assert np.array_equal(np.asarray(line_d.get_ydata()),dv)
        ax.set(title=title,xlabel='Age of household head',ylabel=ylabel,
               xticks=[20,35,50,65,80]);ax.grid(alpha=.16)
        mseries[mkey]=dict(model=mv.tolist(),data=dv.tolist())
    axs[0].legend(frameon=False);fig.tight_layout(rect=(0,.065,1,1))
    fig.text(.5,.015,'Children at home: model dependents; ACS resident own children under 18. Permanent preference shock; transition not converged.',
             ha='center',fontsize=8.5)
    for ext in ('pdf','png'):fig.savefig(out/f'lifecycle_2023.{ext}',dpi=150,bbox_inches='tight')
    plt.close(fig)
    write_qa('lifecycle_2023',dict(ages=x.tolist(),series=mseries),
             dict(acs_2023_age_housing_levels=dict(path=str(empirical_path),sha256=empirical_hash)))
    # Canonical Intergenerational Allocation Model vs Data: six large-owner groups.
    masses=np.zeros((3,2))
    for row in prof['large_owner_age_cells']:
        for j,(lo,hi) in enumerate([(22,40),(40,60),(60,86)]):
            frac=max(0,min(row['age']+row['age_width'],hi)-max(row['age'],lo))/row['age_width']
            masses[j]+=[frac*row['without_children'],frac*row['with_children']]
    mv=100*masses.ravel()/masses.sum();dv=100*np.asarray([float(r['share']) for r in six])
    assert abs(dv.sum()-100)<1e-8 and len(dv)==6
    fig,ax=plt.subplots(figsize=(11,4.5));xx=np.arange(6);width=.36
    for offset,values,color,label in [(-width/2,mv,BLUE,'Model'),(width/2,dv,RED,'ACS 2023')]:
        bars=ax.bar(xx+offset,values,width,color=color,label=label)
        for bar,value in zip(bars,values):
            assert np.isclose(bar.get_height(),value,rtol=0,atol=1e-12)
            ax.text(bar.get_x()+width/2,value+.5,f'{value:.1f}',ha='center',fontsize=9)
    ax.set(xticks=xx,xticklabels=[f'{a}\n{b}' for a in ['Young (22-39)','Middle (40-59)','Old (60-85)']
                                   for b in ['No children','Children']],
           ylabel='Share of large owner-occupied homes (%)',ylim=(0,max(mv.max(),dv.max())*1.2))
    ax.legend(frameon=False);ax.tick_params(axis='x',labelsize=9);fig.tight_layout(rect=(0,.075,1,1))
    fig.text(.5,.015,'Permanent preference shock; transition not converged.',ha='center',fontsize=8.5)
    for ext in ('pdf','png'):fig.savefig(out/f'intergenerational_allocation_2023.{ext}',dpi=150,bbox_inches='tight')
    plt.close(fig)
    write_qa('intergenerational_allocation_2023',dict(model=mv.tolist(),data=dv.tolist(),rooms_threshold=6),
             dict(large_owner_data=dict(path=str(six_path),sha256=six_hash)))
    (out/'permanent_profile_figures_verification.json').write_text(json.dumps(dict(status='PASS',figures=qa),indent=2)+'\n')
    print(json.dumps(dict(figures_dir=str(out),figures=list(qa),model_source_sha256=model_hash),indent=2))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=BASE,help='Patch packet root (default: frozen patch_readout).')
    ap.add_argument('--pdf',type=Path,default=None,help='Combined review PDF destination.')
    ap.add_argument('--sequence-base',type=Path,help='Opt-in carried-history fertility readout; leaves the patch packet unchanged.')
    ap.add_argument('--permanent-profile',type=Path,help='Opt-in permanent-shock profile JSON or directory containing model_2023.json.')
    ap.add_argument('--figures-dir',type=Path,help='Output directory for --permanent-profile figures and QA sidecars.')
    args=ap.parse_args()
    if args.permanent_profile is not None:
        if args.figures_dir is None: ap.error('--figures-dir is required with --permanent-profile')
        permanent_profile_figures(args.permanent_profile,args.figures_dir)
        return
    if args.sequence_base is not None:
        sequence_fertility(args.base,args.sequence_base)
        if (args.sequence_base/'source/readout_2023/model_2023.json').exists():sequence_validation(args.sequence_base)
        if args.pdf:
            from pypdf import PdfReader,PdfWriter
            writer=PdfWriter()
            for name in ('historical_fertility','validation_2023'):
                packet=PdfReader(args.sequence_base/'figures'/f'{name}.pdf')
                assert len(packet.pages)==1
                writer.add_page(packet.pages[0])
            args.pdf.parent.mkdir(parents=True,exist_ok=True)
            with args.pdf.open('wb') as handle:writer.write(handle)
        return
    base=args.base; source=base/'source';data=base/'data';out=base/'figures';out.mkdir(parents=True,exist_ok=True)
    path=csvread(source/'expected_transition.csv');fert=read(source/'fertility.json');static=read(source/'stationary_history.json');initial=read(source/'initial.json')
    r=next(r for r in path if int(r['calendar_year'])==2023)
    root=read(source/'root_receipt.json');assert root['finite_horizon_market_fiscal_converged']
    # The selected patch's preference is source-controlled by the forecast
    # receipt; the frozen default packet retains its original value, while an
    # opt-in refreshed packet may use a newly fitted final shock.
    assert abs(float(path[0]['psi_child'])-float(root['psi']))<1e-14
    full=read(source/'model_2023.json') if (source/'model_2023.json').exists() else None
    if full:
        check=read(source/'verification.json');assert check['status']=='PASS' and full['calendar_year']==2023
        assert full['forecast_receipt_sha256']==hashlib.sha256((source/'root_receipt.json').read_bytes()).hexdigest()
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':12,'legend.fontsize':10})
    pdf=args.pdf or (ROOT/'output/pdf/e5f_patch_review.pdf');pdf.parent.mkdir(parents=True,exist_ok=True);pages=PdfPages(pdf);manifest={}
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
