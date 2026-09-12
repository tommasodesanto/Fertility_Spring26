#!/usr/bin/env python3
"""Render the controlled timing comparison from completed aggregate receipts."""
import csv
import json
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'code/data/psid_followup_mar2026/output/first_birth_correction_review'
ARMS = ['original_native','original_common','aligned_common']
RESULTS = OUT/'timing_local'


def read(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def curve(arm):
    rows = read(RESULTS/arm/'coefficients.csv')
    values = {r['coefficient']:float(r['estimate']) for r in rows}
    covariance = {(r['coefficient_i'],r['coefficient_j']):float(r['covariance'])
                  for r in read(RESULTS/arm/'covariance.csv')}
    values['F2event'] = 0.0
    def event(name):
        return int(name[1:-5]) * (-1 if name.startswith('F') else 1)
    names = sorted(values,key=event)
    x = np.array([event(n) for n in names])
    y = np.array([values[n] for n in names])
    se = np.sqrt([max(0,covariance.get((n,n),0)) for n in names])
    ref = values['F1event']
    rebased_se = np.sqrt([max(0,covariance.get((n,n),0)+covariance['F1event','F1event']
                              -2*covariance.get((n,'F1event'),0)) for n in names])
    return x,y,se,y-ref,rebased_se


def main(baseline_only=False):
    arms = ['original_native'] if baseline_only else ARMS
    receipts = {a:json.loads((RESULTS/a/'run_receipt.json').read_text()) for a in arms}
    assert all(r['status']=='pass' for r in receipts.values())
    if not baseline_only:
        assert receipts['original_common']['sample_keys_sha256']==receipts['aligned_common']['sample_keys_sha256']
    assert len({r['estimator_do_sha256'] for r in receipts.values()})==1
    curves = {a:curve(a) for a in arms}
    fit = {a:read(RESULTS/a/'fit_receipt.csv')[0] for a in arms}
    may = read(OUT/'may_plot_estimates.csv')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig, axes = plt.subplots(1,1 if baseline_only else 3,figsize=(8,5) if baseline_only else (15,4.6),sharex=True)
    axes = np.atleast_1d(axes)
    navy='#285a7d'; orange='#b76539'; gray='#777777'
    axes[0].plot([float(r['relative_time']) for r in may],[float(r['b']) for r in may],
                 'o--',color=gray,lw=1.4,ms=3,label='Saved May figure')
    x,y,se,_,_=curves['original_native']
    axes[0].plot(x,y,'o-',color=navy,lw=1.6,ms=3,label='Recognized code, current input')
    axes[0].fill_between(x,y-1.96*se,y+1.96*se,color=navy,alpha=.10)
    axes[0].set_title('Reproduction of the original specification')
    comparisons = [] if baseline_only else [('original_common',gray,'Original rooms assignment'),('aligned_common',orange,'Rooms moved to next interview')]
    for arm,color,label in comparisons:
        x,y,se,rebased,rse=curves[arm]
        axes[1].plot(x,y,'o-',color=color,lw=1.6,ms=3,label=label)
        axes[1].fill_between(x,y-1.96*se,y+1.96*se,color=color,alpha=.08)
        axes[2].plot(x,rebased,'o-',color=color,lw=1.6,ms=3,label=label)
        axes[2].fill_between(x,rebased-1.96*rse,rebased+1.96*rse,color=color,alpha=.08)
    if not baseline_only:
        axes[1].set_title('Timing only, identical observations')
        axes[2].set_title('Same estimates, each measured from year −1')
    for ax in axes:
        ax.axhline(0,color='#aaaaaa',lw=.7)
        ax.axvline(-.5,color='#aaaaaa',lw=.7,ls=':')
        ax.set_xlabel('Years relative to first birth; tails pooled')
        ax.set_ylabel('Rooms coefficient')
        ax.set_xticks([-7,-5,-3,-1,1,3,5,7,9,11])
        ax.legend(frameon=False,fontsize=8,loc='best')
    if not baseline_only:
        fig.suptitle('Controlled rooms-timing diagnostic — original sample definition, controls, codes and estimator',fontsize=13)
    fig.text(.5,.02,'Diagnostic only. Original −2/−6 reference restrictions and last-cohort comparison remain. Shading: 95% intervals.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.06,1,.92])
    stem = 'original_reproduction' if baseline_only else 'timing_only_comparison'
    fig.savefig(OUT/(stem+'.png'),dpi=180)
    fig.savefig(OUT/(stem+'.pdf'))
    summaries=[]
    for a in arms:
        x,y,*_=curves[a]
        item={'arm':a,'observations':int(float(fit[a]['observations'])),'clusters':int(float(fit[a]['clusters'])),
              'coefficient_m1':float(y[x==-1][0]),'coefficient_p3':float(y[x==3][0]),
              'contrast_p3_m1':float(fit[a]['contrast_p3_m1']),'contrast_se':float(fit[a]['contrast_se']),
              'runtime_seconds':float(fit[a]['runtime_seconds'])}
        assert abs(item['coefficient_p3']-item['coefficient_m1']-item['contrast_p3_m1'])<1e-10
        summaries.append(item)
    with (OUT/(stem+'_summary.csv')).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(summaries[0]),lineterminator='\n');writer.writeheader();writer.writerows(summaries)
    (OUT/(stem+'_verification.json')).write_text(json.dumps({'identical_common_sample':None if baseline_only else True,'summaries':summaries},indent=2)+'\n')
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-label',default='')
    parser.add_argument('--baseline-only',action='store_true')
    args=parser.parse_args()
    RESULTS=RESULTS/args.run_label
    main(args.baseline_only)
