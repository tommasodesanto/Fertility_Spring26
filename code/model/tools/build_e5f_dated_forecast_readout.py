"""Dated diagnostic figures from saved forecast receipts; never fit or solve."""
from pathlib import Path
import argparse, csv, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
SRC=ROOT/'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/forecast_28_completed'
SUB='Single permanent-shock diagnostic; historical path not fitted'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pdf',action='store_true');a=ap.parse_args()
    rows=list(csv.DictReader((SRC/'expected_transition.csv').open()))
    fert=json.loads((SRC/'fertility.json').read_text())
    receipt=json.loads((SRC/'root_receipt.json').read_text())
    summary=json.loads((SRC/'summary.json').read_text())
    assert receipt['finite_horizon_market_fiscal_converged'] and not receipt['terminal_distance_passed']
    assert not summary['historical_fit_complete'] and not summary['horizon_verified']
    assert len(rows)==len(fert)==28
    assert [int(r['calendar_year']) for r in rows]==[int(r['calendar_year']) for r in fert]
    selected=[r for r in rows if int(r['calendar_year'])==2023];assert len(selected)==1
    r=selected[0];values={k:float(r[k]) for k in ['asset_price','renter_price','housing_demand','housing_supply','owner_rate','pension_period_units']}
    assert abs(values['housing_demand']/values['housing_supply']-1)<2e-4
    assert 0<=values['owner_rate']<=1
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    def save(fig,name):
        fig.savefig(SRC/f'{name}.png',dpi=160)
        if a.pdf:fig.savefig(SRC/f'{name}.pdf')
        plt.close(fig)
    fig,ax=plt.subplots(1,3,figsize=(12,4.6))
    ax[0].axis('off');table=ax[0].table(cellText=[['House asset price',f"{values['asset_price']:.4f}"],['Rent per room',f"{values['renter_price']:.4f}"],['Period pension',f"{values['pension_period_units']:.4f}"]],colLabels=['Price / benefit','Model units'],loc='center',cellLoc='left',colWidths=[.68,.32]);table.auto_set_font_size(False);table.set_fontsize(10);table.scale(1,1.8)
    v=[values['housing_demand'],values['housing_supply']];bars=ax[1].bar(['Demand','Supply'],v,color=['#b04a4a','#4776a8']);ax[1].set(title='Housing market',ylabel='Physical rooms (model population units)',ylim=(0,max(v)*1.2))
    for bar,n in zip(bars,v):ax[1].text(bar.get_x()+bar.get_width()/2,n+.12,f'{n:.4f}',ha='center',fontsize=10)
    v=[100*values['owner_rate'],100*(1-values['owner_rate'])];bars=ax[2].bar(['Owners','Renters'],v,color=['#b04a4a','#4776a8']);ax[2].set(title='Tenure allocation',ylabel='Percent of household heads',ylim=(0,100))
    for bar,n in zip(bars,v):ax[2].text(bar.get_x()+bar.get_width()/2,n+1,f'{n:.1f}%',ha='center')
    fig.suptitle('2023 equilibrium: diagnostic forecast',fontsize=16,y=.99);fig.text(.5,.905,SUB,ha='center',fontsize=10);fig.text(.5,.025,'Housing and PAYGO budgets balance on this finite path; terminal horizon is not certified.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.07,1,.89));save(fig,'equilibrium_2023')
    target={2011:1.974875,2015:1.861,2019:1.755375,2023:1.64575}
    observed={int(f['calendar_year'])+4:float(f['period_tfr_topcode_adjusted']) for f in fert}
    ends=np.array(list(target));model=np.array([observed[int(y)] for y in ends]);data=np.array(list(target.values()))
    assert abs(model[0]-float(summary['model']))<1e-12
    fig,ax=plt.subplots(figsize=(9,5));ax.plot(ends,model,'o-',color='#b04a4a',label='Model forecast');ax.plot(ends,data,'s--',color='#4776a8',label='Published female TFR: four-year average')
    ax.scatter([2007],[2.1],marker='D',color='gray',label='Initial model normalization: 2.1')
    ax.set(xticks=[2007,*ends],xlabel='End of four-year birth window',ylabel='Period fertility',ylim=(1.55,2.17));ax.grid(alpha=.18);ax.legend(frameon=False,loc='lower left',fontsize=9)
    fig.suptitle('Historical fertility: data and diagnostic forecast',fontsize=15,y=.99);fig.text(.5,.925,SUB,ha='center',fontsize=10)
    fig.text(.5,.035,'Model: household-rate analogue of female TFR. Initial 2.1 is a normalization, not an observed data point.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.08,1,.91));save(fig,'historical_fertility')
    result=dict(calendar_year=2023,equilibrium=values,window_ends=ends.tolist(),model_fertility=model.tolist(),data_fertility=data.tolist(),finite_converged=True,horizon_verified=False,historical_fit_complete=False,source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [SRC/'expected_transition.csv',SRC/'fertility.json',SRC/'root_receipt.json',SRC/'summary.json']})
    (SRC/'figure_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
