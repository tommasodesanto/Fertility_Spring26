"""Supplemental native forecast, never a fitted realized history."""
from pathlib import Path
import csv,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
p=Path(__file__).resolve().parent
rows=list(csv.DictReader((p/'expected_transition.csv').open()))
fertility=json.loads((p/'fertility.json').read_text())
receipt=json.loads((p/'root_receipt.json').read_text());assert receipt['finite_horizon_market_fiscal_converged']
years=np.array([int(r['calendar_year']) for r in rows]);series=lambda k:np.array([float(r[k]) for r in rows])
fig,axs=plt.subplots(2,2,figsize=(11,8))
axs[0,0].plot(years+4,[r['period_tfr_topcode_adjusted'] for r in fertility],'o--',label='Expected path: current preference stays fixed')
axs[0,0].plot([2011,2015,2019,2023],[1.974875,1.861,1.755375,1.64575],'ks-',label='Data: four-year average')
axs[0,0].set(title='Fertility: diagnostic forecast versus data',xlabel='End of birth window',ylabel='Period fertility');axs[0,0].legend(fontsize=8)
axs[0,1].plot(years,series('asset_price'),'o-');axs[0,1].set(title='Expected house price',xlabel='Year',ylabel='Model price units')
axs[1,0].plot(years,100*series('owner_rate'),'o-');axs[1,0].set(title='Expected ownership',xlabel='Year',ylabel='Percent of household heads')
axs[1,1].semilogy(years,np.maximum(np.abs(receipt['final']['fiscal_residual']),1e-16),'o-',label='PAYGO relative budget gap')
axs[1,1].axhline(1e-6,color='k',linestyle=':',label='Unchanged tolerance');axs[1,1].set(title='Pension budget verification',xlabel='Year',ylabel='Absolute relative residual');axs[1,1].legend(fontsize=8)
for ax in axs.flat:
    ax.grid(alpha=.2);ax.set_xticks(years+4 if ax is axs[0,0] else years);ax.tick_params(axis='x',labelsize=9)
fig.suptitle('Six-date native smoke: markets and pensions pass; horizon does not',fontsize=14)
fig.text(.5,.02,'Single permanent preference decline of 0.045. No later surprises are anticipated. Not a fitted historical path.\nFertility retains the documented household-rate approximation to female TFR.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.065,1,.95));fig.savefig(p/'native_forecast.png',dpi=160);fig.savefig(p/'native_forecast.pdf');plt.close(fig)
