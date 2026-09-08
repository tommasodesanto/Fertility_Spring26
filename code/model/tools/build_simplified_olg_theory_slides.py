#!/usr/bin/env python3
"""Build the current seven-slide theory extract and its illustrative diagrams.

No equilibrium or calibration is solved. Utility curves illustrate the proved
local reallocation. Fertility paths are assumed; population follows exactly
from the cohort law. They are not computed policy equilibrium paths.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/model/simplified_olg_amendments'
TMP=ROOT/'tmp/pdfs/utilitarian_slides'
PDF=ROOT/'output/pdf'
SOURCE=ROOT/'latex/september_14_presentation.tex'
TMP.mkdir(parents=True,exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(TMP/'mplconfig'))
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
BLUE,RED,GRAY='#26527a','#b0493d','#737373'


def save(fig,name):
    for extension in ('pdf','png'):
        fig.savefig(OUT/f'{name}.{extension}',dpi=160,bbox_inches='tight')
    plt.close(fig)


def style():
    plt.rcParams.update({'font.family':'serif','font.serif':['DejaVu Serif'],
        'mathtext.fontset':'cm','font.size':15,'axes.spines.top':False,
        'axes.spines.right':False,'axes.linewidth':1.,'pdf.fonttype':42,
        'legend.frameon':False})


def allocation_figure():
    # Schematic bundles, not a numerical equilibrium or proof by example.
    alpha,gamma,child_space=1.,.70,.30
    hy,ho,move=.90,1.40,.12
    young=lambda h:alpha/(h-child_space)
    old=lambda h:gamma/h
    gain=alpha*np.log1p(move/(hy-child_space))+gamma*np.log1p(-move/ho)
    assert hy>child_space and ho>move and gain>0 and young(hy)>old(ho)
    fig,ax=plt.subplots(figsize=(6.6,4.2),layout='constrained')
    grid=np.linspace(.68,1.72,350)
    ax.plot(grid,young(grid),color=BLUE,lw=2.7,label='Young owner')
    ax.plot(grid,old(grid),color=RED,lw=2.7,ls=(0,(6,3)),label='Old owner')
    for h,val,color in ((hy,young(hy),BLUE),(ho,old(ho),RED)):
        ax.plot(h,val,'o',color=color,ms=7)
        ax.plot([h,h],[.30,val],color=color,lw=.9,ls=':')
    for h0,h1,fun,color in ((hy,hy+move,young,BLUE),(ho,ho-move,old,RED)):
        ax.annotate('',xy=(h1,fun(h1)),xytext=(h0,fun(h0)),
            arrowprops=dict(arrowstyle='->',color=color,lw=1.8,
                            mutation_scale=16,shrinkA=5,shrinkB=0))
    middle=1.14
    ax.plot([hy,middle],[young(hy),young(hy)],color=GRAY,lw=.9,ls=':')
    ax.plot([middle,ho],[old(ho),old(ho)],color=GRAY,lw=.9,ls=':')
    ax.annotate('',xy=(middle,old(ho)),xytext=(middle,young(hy)),
        arrowprops=dict(arrowstyle='<->',color=GRAY,lw=1.1))
    ax.text(middle-.035,1.07,'Utility\ngap',ha='right',va='center',color=GRAY,fontsize=13)
    ax.set(xlim=(.65,1.75),ylim=(.30,2.85),xlabel='Housing',ylabel='Marginal utility of housing')
    ax.set_xticks([hy,ho],[r'$h$',r'$h^2$'])
    ax.set_yticks([])
    ax.legend(loc='upper right',fontsize=13)
    save(fig,'theory_slides_utilitarian_allocation')
    return {'scope':'Schematic raw-utility curves at fixed consumption, fertility and estates; not an equilibrium or a compensation calculation.',
            'young_initial_marginal_utility':young(hy),'old_initial_marginal_utility':old(ho),
            'housing_transfer':move,'utility_change':float(gain),
            'illustrative_values':dict(alpha=alpha,gamma=gamma,child_space=child_space,hy=hy,ho=ho)}


def transition_figure():
    # Assumed fertility paths; only the population accounting is a model identity.
    dates=np.arange(-4,181)
    policy_date=5
    base=np.ones(len(dates))
    after_shock=dates>=0
    base[after_shock]=1-.14*np.exp(-dates[after_shock]/6)
    policy=base.copy()
    after_policy=dates>=policy_date
    policy[after_policy]+=.035*np.exp(-(dates[after_policy]-policy_date)/6)
    y=np.ones((2,len(dates)))
    for k,growth in enumerate((base,policy)):
        for j in range(len(dates)-1):
            y[k,j+1]=growth[j]*y[k,j]
    old=np.concatenate((np.ones((2,1)),y[:,:-1]),axis=1)
    population=(y+old)/2
    assert np.all(policy>=base) and np.all(policy<=1)
    assert np.array_equal(y[0,dates<=policy_date],y[1,dates<=policy_date])
    assert np.all(population[1]>=population[0])
    ratios=np.concatenate(([1.],np.cumprod(policy[:-1]/base[:-1])))
    error=float(np.max(abs(y[1]/y[0]-ratios)))
    assert error<2e-14
    fig,axes=plt.subplots(1,2,figsize=(12.,4.25))
    fig.subplots_adjust(left=.075,right=.985,bottom=.19,top=.81,wspace=.30)
    shown=dates<=34
    for ax in axes:
        ax.axvline(0,color=GRAY,lw=.8,ls=':')
        ax.axvline(policy_date,color=GRAY,lw=.8,ls=':')
        ax.set_xlim(-4,34)
        ax.set_xticks([0,policy_date],[r'$t_0$',r'$t_p$'])
        ax.set_xlabel('Date')
    # Duplicate dates locate the two discrete changes at their stated dates.
    def fertility_line(values, with_policy=False):
        plot_dates, plot_values = [], []
        for date, value in zip(dates[shown], values[shown]):
            if date == 0:
                plot_dates.append(date)
                plot_values.append(1.)
            if with_policy and date == policy_date:
                plot_dates.append(date)
                plot_values.append(base[dates == date][0])
            plot_dates.append(date)
            plot_values.append(value)
        return plot_dates, plot_values
    axes[0].plot(*fertility_line(base),color=BLUE,lw=2.5,label='Without intervention')
    axes[0].plot(*fertility_line(policy,True),color=RED,lw=2.5,ls=(0,(6,3)),label='With intervention')
    axes[0].axhline(1,color=GRAY,lw=1.,ls=':')
    axes[0].set_ylim(.845,1.025)
    axes[0].set_yticks([1],[r'$1/\nu$'])
    axes[0].set_ylabel('Fertility')
    axes[0].set_title('Fertility comparison',fontsize=15,pad=10)
    for k,color,ls in ((0,BLUE,'-'),(1,RED,(0,(6,3)))):
        axes[1].plot(dates[shown],population[k,shown],color=color,lw=2.5,ls=ls)
        axes[1].axhline(population[k,-1],color=color,lw=.8,ls=':')
    axes[1].set_ylim(.35,1.07)
    axes[1].set_yticks([population[0,-1],population[1,-1],1],
                      [r'$N_B^\infty$',r'$N_P^\infty$',r'$N_0$'])
    axes[1].set_ylabel('Adult households')
    axes[1].set_title('Population implied by fertility',fontsize=15,pad=10)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.52,1.01),
               ncol=2,fontsize=14,columnspacing=2.8,handlelength=2.6)
    save(fig,'theory_slides_utilitarian_transition')
    return {'scope':'Assumed fertility decline and later positive policy fertility effect; exact cohort accounting. Not a solved equilibrium or a demonstrated policy response.',
        'shock_date':0,'policy_date':policy_date,'same_inherited_cohorts':True,
        'cohort_ratio_identity_max_error':error,
        'baseline_limit_relative_to_initial':float(population[0,-1]),
        'policy_limit_relative_to_initial':float(population[1,-1]),
        'baseline_fertility_tail_bound':float(.14*np.exp(-180/6)/(1-np.exp(-1/6))),
        'limit_scope':'Geometrically summable fertility deficits yield positive cohort-product limits; these are not established stationary equilibria of the household model.'}


def compile_decks():
    source=SOURCE.read_text()
    preamble=source.split(r"\begin{document}",1)[0]
    block=source.split(r"\section{Simple Theory}",1)[1].split(r"\section{Quantitative Model}",1)[0]
    block=block.split(r"\end{frame}",1)[1]  # Omit the section divider.
    assert block.count(r"\begin{frame}")==7
    block=re.sub(r"\\hyperlink\{[^}]+\}\{\\beamerbutton\{[^}]*\}\}", "", block)
    extract=TMP/"simplified_olg_theory_slides.tex"
    extract.write_text(preamble+r"\hypersetup{pdftitle={Housing Allocation across Generations}}"+"\n"+
                       r"\begin{document}"+"\n"+block+r"\end{document}"+"\n")
    engine=shutil.which("pdflatex") or "/Library/TeX/texbin/pdflatex"
    env=dict(os.environ)
    env["PATH"]=str(Path(engine).parent)+os.pathsep+env.get("PATH","")
    for src in (SOURCE,extract):
        for run in (1,2):
            with (TMP/f"{src.stem}.pass{run}.txt").open("w") as log:
                subprocess.run([engine,"-interaction=nonstopmode","-halt-on-error","-file-line-error",
                                f"-output-directory={TMP}",str(src)],cwd=ROOT/"latex",env=env,
                               stdout=log,stderr=subprocess.STDOUT,check=True)
        shutil.copy2(TMP/f"{src.stem}.pdf",PDF/f"{src.stem}.pdf")
    shutil.copy2(PDF/"september_14_presentation.pdf",SOURCE.with_suffix(".pdf"))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figures-only',action='store_true')
    parser.add_argument('--compile-only',action='store_true')
    args=parser.parse_args()
    if args.figures_only and args.compile_only:
        parser.error('Choose one build scope')
    for folder in (OUT,TMP,PDF):
        folder.mkdir(parents=True,exist_ok=True)
    receipt_path=OUT/'theory_slides_utilitarian_checks.json'
    if not args.compile_only:
        style()
        receipt={'allocation':allocation_figure(),'transition':transition_figure()}
    else:
        receipt=json.loads(receipt_path.read_text())
    if not args.figures_only:
        compile_decks()
    receipt['generated_utc']=datetime.now(timezone.utc).isoformat()
    receipt['sources']={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (Path(__file__),SOURCE,ROOT/'latex/JMP_DS_suggestions/simplified_olg_utilitarian.tex')}
    receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
    print('Theory illustrations checked'+('; slide PDFs compiled.' if not args.figures_only else '.'))


if __name__=='__main__':
    main()
