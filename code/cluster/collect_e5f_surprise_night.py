"""Collect the bounded campaign even if some Slurm tasks fail."""
from pathlib import Path
import argparse,csv,json,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def read(p):return json.loads(p.read_text())
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/night_surprises_20260912'));root=parser.parse_args().root
    out=root/'report';out.mkdir(exist_ok=True)
    plan=read(root/'plan.json');score=Path(plan['initial_score_path'])
    shutil.copy2(score,out/'historical_initial_score_full.json')
    fits=[];status=[]
    with PdfPages(out/'overnight_diagnostics.pdf') as pdf:
        fig,ax=plt.subplots(figsize=(10,6));targets=list(csv.DictReader(Path(plan['empirical_blocks']).open()))
        ax.plot([int(t['birth_year_end']) for t in targets],[float(t['period_tfr_arithmetic_mean']) for t in targets],'ko-',label='Data: four-year average')
        for arm in sorted((root/'results').glob('arm_*')):
            current=read(arm/'realized_fit.json') if (arm/'realized_fit.json').exists() else []
            failure=read(arm/'failure.json') if (arm/'failure.json').exists() else None
            status.append(dict(arm=arm.name,accepted_windows=len(current),failure=failure))
            if current:
                ax.plot([r['year']+4 for r in current],[r['model'] for r in current],'o-',label=arm.name+' accepted windows')
            for p in sorted(arm.glob('trial_*/fit.json')):fits.append(dict(arm=arm.name,**read(p)))
        ax.set(xlabel='End of birth window',ylabel='Period fertility',title='Successive surprises: accepted realized windows only')
        ax.legend();ax.grid(alpha=.2);fig.tight_layout();pdf.savefig(fig);fig.savefig(out/'fertility_model_vs_data.png',dpi=160);plt.close(fig)
        for arm in sorted((root/'results').glob('arm_*')):
            casefiles={}
            for case in ('baseline','equal-rebate-1pct','equal-rebate-2pct'):
                files=list((arm/'policies'/case).glob('policy_path_*.csv'))
                if files:casefiles[case]=max(files,key=lambda p:int(p.stem.split('_')[-1]))
            if not casefiles:continue
            rows={case:list(csv.DictReader(p.open())) for case,p in casefiles.items()}
            keys=[k for k in ('asset_price','renter_price','housing_demand','housing_supply','total_persons','births','pension') if all(k in v[0] for v in rows.values())]
            for key in keys:
                fig,ax=plt.subplots(figsize=(10,6))
                for case,values in rows.items():ax.plot([2023+4*i for i in range(len(values))],[float(r[key]) for r in values],label=case)
                ax.set(xlabel='Year',ylabel=key,title=arm.name+': finite-horizon policy diagnostics; horizon not certified')
                ax.legend();ax.grid(alpha=.2);fig.tight_layout();pdf.savefig(fig);plt.close(fig)
        if not fits:
            fig,ax=plt.subplots(figsize=(10,6));ax.axis('off');ax.text(.03,.93,'No historical window passed all acceptance gates.\nSee failure receipts; the data line is not a model fit.',va='top',fontsize=15);pdf.savefig(fig);plt.close(fig)
    (out/'status.json').write_text(json.dumps(status,indent=2)+'\n')
    (out/'all_completed_fertility_fits.json').write_text(json.dumps(fits,indent=2)+'\n')
    cal=Path('/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/night_20260912/search')
    for p in cal.glob('round_*/results/selected_*'):
        if p.is_file():shutil.copy2(p,out/(p.parents[1].name+'_'+p.name))
    (out/'README.md').write_text('# Overnight evidence\n\nPDF and PNG show only accepted realized windows. All completed fitting trials and failures are saved separately. Full historical seed score is preserved; newer calibration tables are separate and were not substituted into a running history. Standard 17-graph packets remain beside each accepted solution. No path is certified for production without horizon verification.\n')
    print(str(out))
if __name__=='__main__':main()
