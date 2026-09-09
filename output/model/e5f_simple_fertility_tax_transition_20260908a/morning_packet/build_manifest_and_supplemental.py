from pathlib import Path
import json,csv,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
O=R/'output/model/e5f_simple_fertility_tax_transition_20260908a'; P=O/'morning_packet'
def readcsv(p):
 with p.open() as f:return list(csv.DictReader(f))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
rows=readcsv(O/'results/comparison.csv'); assert len(rows)==11
D={k:np.array([float(r[k]) for r in rows]) for k in rows[0]}; years=D['calendar_year']
first,last=rows[0],rows[-1]
def val(row,key):return float(row[key])
def fmt(row,key):return f'{val(row,key):+.3f}'
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2,'savefig.dpi':180})
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
for ax,series,title,unit in [
 (axes[0,0],[('births_per_household_pct_change','Births per household'),('total_adjusted_births_pct_change','Total births')],'Birth flows','% versus 1% tax'),
 (axes[0,1],[('young_mean_rooms_whole_nodes_pct_change','All young households'),('young_dependent_rooms_pct_change','Young parents')],'Occupied housing','% versus 1% tax'),
 (axes[1,0],[('young_ownership_whole_nodes_pp_change','All young households'),('young_dependent_owner_rate_pp_change','Young parents')],'Ownership','Percentage points'),
 (axes[1,1],[('asset_price_pct_change','Housing asset price'),('household_mass_pct_change','Number of households')],'Prices and household mass','% versus 1% tax')]:
 for key,label in series:ax.plot(years,D[key],marker='o',markersize=3,label=label)
 ax.axhline(0,color='grey',lw=.8);ax.set(title=title,ylabel=unit,xlabel='Year');ax.set_xticks(years[::2]);ax.legend(fontsize=9)
fig.suptitle('Supplemental: doubling property tax, with equal rebates in both paths',fontsize=14)
fig.savefig(P/'supplemental/policy_effects.png');plt.close(fig)
cases=['tax1-equal-rebate','tax2-equal-rebate'];paths={case:readcsv(O/'results/full'/case/'path.csv') for case in cases}
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
for ax,key,title in [(axes[0,0],'births_per_household','Births per household (4-year flow)'),(axes[0,1],'household_mass','Household mass (model units)'),(axes[1,0],'housing_demand_per_adult','Occupied rooms per household'),(axes[1,1],'asset_price','Housing asset price (model units)')]:
 for case,label in zip(cases,['Annual 1% tax + rebate','Annual 2% tax + rebate']):ax.plot(years,[float(r[key]) for r in paths[case]],marker='o',markersize=3,label=label)
 ax.set(title=title,xlabel='Year');ax.set_xticks(years[::2]);ax.legend(fontsize=8)
fig.suptitle('Supplemental: levels along the two temporary-equilibrium paths',fontsize=14)
fig.savefig(P/'supplemental/path_levels.png');plt.close(fig)
base=json.loads((P/'calibration_inputs.json').read_text())
m={k:base[k] for k in ['target_fit','parameters','source_files']}
m.update(schema='e5f_rebated_tax_morning_report_v1',title='Rebated property tax: overnight results',date='Morning review • 9 September 2026',subtitle='Fixed simultaneous fertility-nest calibration. Annual property tax rises from 1% to 2%; each path rebates its own revenue equally to households.',footer='Fertility and housing • diagnostic tax comparison • September 2026')
def table(cols,rs,notes=()):return {'columns':[{'key':k,'label':label,'width':w} for k,label,w in cols],'rows':rs,'notes':list(notes)}
summarykeys=[('Births per household','births_per_household_pct_change','%'),('Total births','total_adjusted_births_pct_change','%'),('Young ownership','young_ownership_whole_nodes_pp_change','pp'),('Young occupied rooms','young_mean_rooms_whole_nodes_pct_change','%'),('Young-parent ownership','young_dependent_owner_rate_pp_change','pp'),('Young-parent occupied rooms','young_dependent_rooms_pct_change','%'),('Housing asset price','asset_price_pct_change','%')]
short=[dict(outcome=name,impact=fmt(first,k)+' '+unit,end=fmt(last,k)+' '+unit) for name,k,unit in summarykeys]
m['sections']=[{'title':'Start here', 'paragraphs':[
 'Both tax paths completed all eleven dates, 2023–2063. The selected calibration remains fixed at loss 23.791955; the full fit and every estimated parameter follow. Numerical completion does not resolve the economic validation issues below.',
 f"The initial fertility gain is {fmt(first,'births_per_household_pct_change')}%; by 2063 the birth flow per household changes {fmt(last,'births_per_household_pct_change')}%. Young households occupy less space: {fmt(first,'young_mean_rooms_whole_nodes_pct_change')}% initially and {fmt(last,'young_mean_rooms_whole_nodes_pct_change')}% in 2063.",
 'The mechanism is clear at impact: lower purchase prices ease the down payment, but the ongoing cost per room rises. The rebate supports births; it does not undo the reduction in housing space.',
 'Immediate discussion: (1) the first-birth housing response remains too small; (2) the policy raises ongoing housing costs despite cheaper assets; (3) the future household-entry rule and almost-universal late-life ownership remain unresolved. These deserve attention before a final presentation claim.'
 ],'tables':[table([('outcome','Outcome',3),('impact','2023 impact',1.5),('end','2063',1.5)],short, ['Young = unchanged model age nodes 26, 30 and 34. Young parents have dependent children. These are group averages, with composition effects. Percentage changes and percentage-point changes are distinct.'])]},
 {'title':'The complete dated policy comparison','paragraphs':['Every row compares the 2% tax path with the 1% tax path at the same date. Both return all property-tax revenue through an equal rebate per household. Birth flows include the maintained adjustment for top-coded parity.'], 'tables':[table([('year','Year',.65),('birth','Births/HH %',1.1),('total','Total births %',1.1),('hh','HH mass %',1.1),('own','Young-parent own. pp',1.3),('rooms','Young-parent rooms %',1.3)], [dict(year=int(float(r['calendar_year'])),birth=fmt(r,'births_per_household_pct_change'),total=fmt(r,'total_adjusted_births_pct_change'),hh=fmt(r,'household_mass_pct_change'),own=fmt(r,'young_dependent_owner_rate_pp_change'),rooms=fmt(r,'young_dependent_rooms_pct_change')) for r in rows], ['Household mass is identical across paths before 2043. The first policy-induced births enter the household population in 2043 under the maintained twenty-year entry lag.'])]},
 {'title':'Why ownership can rise while housing space falls','paragraphs':[
 'At impact the asset price falls 4.802%, but the four-year rent/user cost per room rises from 0.12106 to 0.14306: +18.178%. The recurring tax more than offsets capitalization in the flow cost. A cheaper purchase therefore does not imply cheaper housing services.',
 'An exact eight-cell decomposition averages each channel over all orderings of changes in the tax rate, asset price and equal rebate. Mixed cells are conditional household responses, not separate equilibria. The components add to the verified impact effect; they are not long-run decompositions.'
 ],'tables':[]},
 {'title':'Policy effects through 2063','paragraphs':['These supplemental summaries complement the unchanged standard diagnostic set. They use the verified comparison table; no additional model is fitted.'], 'figures':[{'path':str(P/'supplemental/policy_effects.png'),'caption':'Effects relative to the contemporaneous 1% tax path, with equal rebates in both paths. Lines join four-year observations.'}]},
 {'title':'Levels along each path','paragraphs':['The reference path also changes as the inherited age distribution advances and new household cohorts enter. A percentage difference between policies should not be read as a population forecast.'], 'figures':[{'path':str(P/'supplemental/path_levels.png'),'caption':'Model household units and birth flows. The same selected 2023 preference parameters are held fixed throughout.'}]},
 {'title':'What is verified, and what remains open','paragraphs':[
 'Verified: unchanged scientific source and target fingerprints; common inherited population and entry queues; reproduction of both impact equilibria and both two-date smoke paths; market clearing, balanced rebate ledgers, mass and queue accounting, feasibility, choice probabilities and occupied-state value monotonicity at each date. Full saved checkpoints remain on Torch.',
 'Reporting repairs: ownership-policy graphs now sum all five owned sizes and weight conception outcomes; first-birth hazards use the pre-choice population at risk. The corrections reproduce the native transition operator. Fifteen income states have distinct colors; housing and ownership rise with permanent income in independently reconstructed group totals.',
 'Still open: the first-birth housing response is 0.4439 against 0.7202 rooms. Late-life ownership is close to 100%. Young ownership is 33.1957% in the selected calibration versus an ACS 34.1166%, but exact annual-age alignment remains unresolved. A saved-state check traces ownership-probability dips to renting overtaking four-room ownership before six-room ownership becomes attractive. All inspected origin states have zero population mass. Sensitivity to the numerical grid remains untested; this is neither a demonstrated solver failure nor a proof of grid accuracy.',
 'The forward closure is a closed-national household-unit diagnostic: no outside entry, full retention, births divided by 2.1 for household formation, and a twenty-year lag. Household headship and the handoff from historical empirical age weights to endogenous entrants remain unresolved. It is a sequence of temporary equilibria, with each date’s prices perceived as permanent; it is not a perfect-foresight path, resident-population forecast or welfare calculation.',
 'No calibration, weights, target definitions, policy closure, numerical tolerance or production benchmark was changed during this monitoring run. The standard appendix shows both 2063 endpoints; the companion folder contains all 17 standard graphs for every date and both policies (374 images), plus audit receipts.'
 ]}]
channel=readcsv(R/'output/model/e5f_simple_fertility_tax_channels_20260908a/results/shapley_components.csv')
ch=[]
for metric,label,pp in [('births_per_household','Births/HH (%)',False),('rooms_per_household','Rooms/HH (%)',False),('young_ownership_whole_nodes','Young ownership (pp)',True),('young_mean_rooms_whole_nodes','Young rooms (%)',False)]:
 subset=[r for r in channel if r['metric']==metric]; d={r['component']:float(r['contribution_percentage_points' if pp else 'contribution_percent_of_baseline']) for r in subset}
 ch.append(dict(outcome=label,tax=f"{d['tax_rate']:+.3f}",price=f"{d['asset_price']:+.3f}",rebate=f"{d['equal_rebate']:+.3f}",total=f"{sum(d.values()):+.3f}"))
m['sections'][2]['tables']=[table([('outcome','Outcome',2.5),('tax','Tax',1),('price','Price',1),('rebate','Rebate',1),('total','Net',1)],ch)]
names=json.loads((O/'graphs/v2/smoke/tax1-equal-rebate/date_2023/graph_manifest.json').read_text())['standard_filenames']
m['diagnostic_groups']=[]
for name in names:
 caption='Both policy paths at 2063. All other dates are in the companion graph set.'
 if name.startswith('policy_'):caption+=' Consumption and housing are conditional on renting with no birth; fertility is the attempt probability. Ownership integrates over all housing products and conception outcomes.'
 if name=='fertility_by_age.png':caption+=' First births divided by the pre-choice mass at risk, per four-year period.'
 if name=='fertility_policy_by_age_income_state.png':caption+=' Attempt probabilities differ from realized birth probabilities because conception is uncertain.'
 if name=='income_state_outcomes.png':caption+=' Income states retain their permanent-group ordering; the horizontal values are not globally sorted. Parity is top-coded at 3+.'
 m['diagnostic_groups'].append({'title':name[:-4].replace('_',' ').capitalize(),'caption':caption,'layout':'one_per_page' if name.startswith('policy_') else 'two_per_page','figures':[{'path':str(O/'graphs/v2/full'/case/'date_2063/standard_diagnostics'/name),'label':label+' • 2063'} for case,label in zip(cases,['Annual 1% tax + equal rebate','Annual 2% tax + equal rebate'])]})
for p in [O/'results/comparison.csv',O/'results/comparison_receipt.json',P/'verification.json',R/'output/model/e5f_simple_fertility_tax_channels_20260908a/results/shapley_components.csv',O/'income_audit/income_audit.json',O/'owner_shape_audit/extraction_status.json']:
 m['source_files'].append({'path':str(p),'sha256':digest(p)})
(P/'report_manifest.regenerated_draft.json').write_text(json.dumps(m,indent=2)+'\n')
print('Manifest prepared; source verification and graph checks required before build.')
