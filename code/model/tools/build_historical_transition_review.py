"""Build a saved-output transition slide review; never solve or alter targets."""
from pathlib import Path
import csv
import os
import json
import hashlib
import numpy as np
os.environ.setdefault('MPLCONFIGDIR', '/tmp/fertility_transition_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/'output/model/e5f_matched_pf_20260909a'
WORK = BASE/'current_candidate_transition'
OUT = WORK/'slide_review'
OUT.mkdir(exist_ok=True)
def rows(p):
    return list(csv.DictReader(p.open()))
path = rows(WORK/'return_home_20260911/paths/delta_m005/transition_path.csv')
path = [r for r in path if int(r['calendar_year']) <= 2023]
annual = rows(BASE/'path_pilot_20260910/fertility_data/annual_fertility_2007_2023.csv')
blocks = rows(WORK/'inputs/empirical_blocks.csv')
mappings = json.loads((WORK/'return_home_20260911/paths/delta_m005/dated_household_fertility.json').read_text())['mappings']
receipt = json.loads((WORK/'return_home_20260911/paths/delta_m005/root_receipt.json').read_text())
final = mappings[receipt['final']['payload']['trial']-1]
assert final == mappings[receipt['best']['payload']['trial']-1]
model = {r['calendar_year']:r['diagnostics']['period_tfr_topcode_adjusted'] for r in final}
fit = rows(BASE/'initial_calibration_contract/extended_refinement/collected_17378993/selected_target_fit.csv')
assert len(fit) == 13
BLUE='#234f7d'; RED='#b74634'; GRAY='#777777'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
 'axes.spines.right':False,'axes.labelcolor':'#333333','text.color':'#222222',
 'axes.titleweight':'normal','savefig.facecolor':'white','figure.facecolor':'white'})
def finish(fig,name):
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(OUT/(name+'.png'),dpi=160,bbox_inches='tight')
    plt.close(fig)
fig,(a,b)=plt.subplots(1,2,figsize=(12.2,4.7),gridspec_kw={'width_ratios':[1,2.05]})
years=np.array([int(r['calendar_year']) for r in path]);psi=np.array([float(r['psi_child']) for r in path])
a.plot([2003,2007], [100,100],color=BLUE,lw=2.3)
a.plot(years,100*psi/psi[0],color=BLUE,lw=2.3,marker='o',ms=4)
a.plot([2023,2026], [100*psi[-1]/psi[0]]*2,color=BLUE,lw=2.3)
a.axvline(2007,color=GRAY,ls=':',lw=1)
a.text(2007.5,97,'Path announced\nin 2007',fontsize=9,va='top')
a.set(xlim=(2002,2026),ylim=(65,104),xticks=[2003,2007,2015,2023],ylabel='Child-preference parameter (initial = 100)',title='Announced preference decline')
a.grid(axis='y',alpha=.17)
by=np.array([(int(r['birth_year_start'])+int(r['birth_year_end']))/2 for r in blocks])
data=np.array([float(r['period_tfr_arithmetic_mean']) for r in blocks]);mv=np.array([model[int(r['decision_year'])] for r in blocks])
b.plot([int(r['year']) for r in annual],[float(r['period_tfr_births_per_woman']) for r in annual],color='#b7b7b7',lw=1.3,label='Annual data')
b.plot(by,data,'s-',color='#222222',lw=1.8,ms=5,label='Data: four-year average')
b.errorbar(by,mv,xerr=1.5,fmt='o-',color=BLUE,lw=2,elinewidth=.8,capsize=2,ms=5,label='Model: four-year birth flow')
b.plot([2003,2006.5],[2.1,2.1],color=BLUE,lw=2,ls='--')
b.text(2003,2.13,'Initial steady-state\nnormalization: 2.1',fontsize=9,color=BLUE)
b.axvline(2007,color=GRAY,ls=':',lw=1)
b.annotate(f'{mv[-1]:.3f}',(by[-1],mv[-1]),xytext=(9,4),textcoords='offset points',color=BLUE,fontsize=10)
b.annotate(f'{data[-1]:.3f}',(by[-1],data[-1]),xytext=(9,-13),textcoords='offset points',color='#222222',fontsize=10)
b.set(xlim=(2002.5,2025),ylim=(1.56,2.22),xticks=[2003,2007,2011,2015,2019,2023],ylabel='Births per woman / model fertility index',title='Fertility: model and data')
b.grid(axis='y',alpha=.17)
b.legend(loc='lower left',fontsize=9,frameon=False)
fig.subplots_adjust(wspace=.3,bottom=.12,top=.88)
finish(fig,'fertility_transition')
comparison=[]
for r,m,d in zip(blocks,mv,data):
    comparison.append({'birth_window':f"{r['birth_year_start']}-{r['birth_year_end']}",'decision_year':r['decision_year'],'data':d,'model':m,'gap':m-d})
with (OUT/'fertility_comparison.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(comparison[0]));w.writeheader();w.writerows(comparison)
# National ACS heads18–85; recorded rooms retain the historical source coding.
acs_file = ROOT/'code/data/Spatial_aggregate_withmicrodata/output/national_householder_housing_path/national_householder_housing_path.csv'
acs = rows(acs_file)
assert [int(r['calendar_year']) for r in acs] == list(years)
for r in acs:
    assert np.isclose(float(r['owner_weight'])/float(r['tenure_valid_weight']),float(r['ownership_rate']))
    assert np.isclose(float(r['rooms_weighted_sum'])/float(r['rooms_valid_weight']),float(r['mean_rooms_literal']))
heads=np.array([float(r['adult_population'] or r['household_heads']) for r in path])
rooms=np.array([float(r['housing_demand']) for r in path])/heads
ownership=100*np.array([float(r['owner_rate']) for r in path])
other=[dict(year=int(y),housing_services_per_household=float(h),ownership_percent=float(o),households_initial100=float(n*100)) for y,h,o,n in zip(years,rooms,ownership,heads)]
with (OUT/'other_outcomes.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(other[0]));w.writeheader();w.writerows(other)
fig,axes=plt.subplots(1,2,figsize=(12.2,4.7))
for ax,values,empirical,title,unit in [(axes[0],rooms,[float(r['mean_rooms_literal']) for r in acs],'Housing services per household','Rooms / model room units'),(axes[1],ownership,[100*float(r['ownership_rate']) for r in acs],'Homeownership','Percent of households')]:
    ax.plot(years,values,'o-',color=BLUE,lw=2.3,ms=5,label='Model')
    ax.plot(years,empirical,'s--',color='#333333',lw=1.8,ms=5,label='ACS: recorded rooms' if ax is axes[0] else 'ACS: heads ages18–85')
    ax.legend(frameon=False,fontsize=9,loc='best')
    ax.set(xticks=years,title=title,ylabel=unit,xlim=(2006,2024))
    ax.grid(axis='y',alpha=.17)
    ax.annotate(f'{values[0]:.2f}',(years[0],values[0]),xytext=(7,7),textcoords='offset points',color=BLUE)
    ax.annotate(f'{values[-1]:.2f}',(years[-1],values[-1]),xytext=(-35,9),textcoords='offset points',color=BLUE)
fig.subplots_adjust(wspace=.28,bottom=.12,top=.88)
finish(fig,'housing_transition')
proof={'model_path':'delta_m005','preference_change':float(psi[-1]-psi[0]),'historical_path_fitted':False,'horizon_verified':False,
 'block_data_decline_percent':float(100*(1-data[-1]/data[0])), 'block_model_decline_percent':float(100*(1-mv[-1]/mv[0])),
 'acs_source_sha256':hashlib.sha256(acs_file.read_bytes()).hexdigest(),
 'existing_initial_fit_table':str(BASE/'initial_calibration_contract/extended_refinement/collected_17378993/selected_target_fit.csv'),
 'existing_parameter_table':str(BASE/'initial_calibration_contract/extended_refinement/collected_17378993/selected_parameters.csv')}
(OUT/'verification.json').write_text(json.dumps(proof,indent=2)+'\n')
print(json.dumps(proof,indent=2));print(other)

with (OUT/'acs_comparison.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=['year','acs_ownership_percent','model_ownership_percent','acs_recorded_rooms','model_room_units']);w.writeheader()
    for a,m in zip(acs,other):
        w.writerow(dict(year=a['calendar_year'],acs_ownership_percent=100*float(a['ownership_rate']),model_ownership_percent=m['ownership_percent'],acs_recorded_rooms=a['mean_rooms_literal'],model_room_units=m['housing_services_per_household']))
