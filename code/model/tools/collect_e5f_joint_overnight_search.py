#!/usr/bin/env python3
"""Finite overnight collector: preserve verified results and draft the morning PDF.

This never submits or restarts a model job. The cluster controller owns the search.
Visual review of the automatically rendered PDF remains explicitly pending.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[3]
REL=Path('output/model/e5f_full_joint_overnight_20260906a')
OUT=ROOT/REL
REMOTE=Path('/scratch/td2248/projects/Fertility_Spring26_full_joint_overnight_20260906a')
PDF_NAME='e5f_full_joint_overnight_review.pdf'
LABELS=['Completed fertility, age 42','Childless share, age 42','Age at first birth',
        'First births at age 30+','First-birth rooms response','Parent rooms: 3+ vs. 1-2',
        'Parent ownership gap','Ownership, ages 30-55','Mean rooms, ages 18-85',
        'Wealth / annual earnings','Bequests / wealth, annual','Old wealth/income: p90/p50']
PARAMETER_LABELS={'beta_annual':'Annual patience','kappa_fert':'First-birth dispersion',
 'kappa_fert_continuation':'Later-birth dispersion','chi':'Owner premium','H0':'Supply intercept',
 'theta0':'Bequest strength','theta1':'Bequest shifter','hbar_child_rooms':'Per-child housing floor',
 'first_birth_fixed_cost':'First-birth fixed cost','hbar_first_child_jump':'First-child housing term',
 'psi_child_change_2023':'Dated child-value change','psi_child_2007':'Old child-value intercept',
 'psi_child_2023':'2023 child-value intercept','tenure_choice_kappa':'Tenure dispersion',
 'housing_supply_elasticity':'Dated supply elasticity'}


def read(path):return json.loads(Path(path).read_text())
def rows(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)
def command(args,timeout=180):
    r=subprocess.run(args,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=timeout)
    if r.returncode:raise RuntimeError(f'{args[0]} exited {r.returncode}: {r.stderr[-1800:]}')
    return r.stdout


def queue_state(job_id):
    # Listing this user's queue also works after Slurm purges a completed job ID.
    queue=command(['ssh','-o','BatchMode=yes','torch','squeue -h -u td2248 -o "%i %T"'],timeout=30)
    return next((line.split()[1] for line in queue.splitlines()
                 if line.split() and line.split()[0]==str(int(job_id))), '')


def pull():
    command(['rsync','-az','--exclude=*.tmp','--exclude=*worker*','--include=*/',
        '--include=*.json','--include=*.csv','--include=MORNING_SUMMARY.md','--exclude=*',
        f'torch:{REMOTE/REL}/',str(OUT)+'/'],timeout=300)


def local(remote_path):
    return ROOT/Path(remote_path).relative_to(REMOTE)


def verify_collected():
    count=0
    for receipt in OUT.rglob('case_receipt.json'):
        data=read(receipt)
        if data['status']!='complete':raise RuntimeError(f'Incomplete receipt {receipt}')
        for name,sha in data['artifact_sha256'].items():
            path=receipt.parent/name
            if path.exists():
                if digest(path)!=sha:raise RuntimeError(f'Collected artifact mismatch {path}')
                count+=1
    return count


def build_report():
    state=read(OUT/'search/search_state.json')
    best=read(OUT/'search/best_so_far.json')['best']
    case=local(best['summary']).parent
    command(['rsync','-az',f"torch:{Path(best['summary']).parent/'standard_diagnostics'}/",str(case/'standard_diagnostics')+'/'])
    verified=verify_collected()
    s=read(case/'summary.json');r=read(case/'case_receipt.json')
    fit=rows(case/'target_fit_long.csv');params=rows(case/'parameter_table.csv')
    if len(fit)!=12 or sum(x['is_free_parameter']=='True' for x in params)!=11:
        raise RuntimeError('Wrong target or estimated-parameter count')
    if abs(sum(float(x['loss_contribution']) for x in fit)-best['loss'])>1e-8:
        raise RuntimeError('Loss does not match the complete table')
    proof_path=OUT/'search/final_verification.json'
    verified_repeat=state['status']=='complete_verified' and proof_path.exists() and read(proof_path)['status']=='pass'
    moments={x['moment']:float(x['model']) for x in fit}
    contract=read(OUT/'contract.json')
    smoke_count=contract.get('smoke_histories',2)
    text=['---','title: "Overnight calibration: results and remaining gaps"',
          'subtitle: "All eleven parameters searched under the unchanged twelve targets"',
          'author: "Prepared for Tommaso De Santo"','date: "6 September 2026"','---','',
          '# The short read','',f"**Best evaluated loss: {best['loss']:.4f}.** "+
          ('Two original-generator repetitions passed exactly.' if verified_repeat else 'Final exact verification is incomplete; this is a provisional result.'),'',
          f"The first-birth housing response is **{moments['housing_increment_0to1']:.4f} rooms against 0.7202**. "+
          f"Mean rooms are **{moments['aggregate_mean_occupied_rooms_18_85']:.4f} against 5.7800**, and ownership is **{100*moments['own_rate']:.2f}% against 57.55%**.",'',
          '| Point | Loss | First-birth rooms | Mean rooms |','|---|---:|---:|---:|',
          '| Retained production | 30.4830 | 0.4364 | 6.3173 |',
          '| Evening, old bounds | 21.7980 | 0.4766 | 6.3325 |',
          '| Housing-only starting diagnostic | 20.6953 | 0.5158 | 6.3491 |',
          f"| Overnight joint search | {best['loss']:.4f} | {moments['housing_increment_0to1']:.4f} | {moments['aggregate_mean_occupied_rooms_18_85']:.4f} |",'',
          'All eleven estimated parameters could move during the overnight search. The first-child housing term could range from 0 to 2 instead of 0 to 0.5; every other bound, all twelve empirical targets and all weights stayed fixed. The model equations and numerical gates were retained.',
          '',f"The controller stopped with status `{state['status']}` after {state['completed_histories']} completed histories, counting the {smoke_count} initial smoke evaluations. Production remains the September 4 calibration. These results have not been used to update policy or population forecasts.",
          '',('The first broad attempt stopped at candidate 17 when the housing-market residual exceeded the unchanged 0.0002 tolerance; none of its fourteen completed trials improved the seed. This separately recorded recovery uses smaller steps: each round probes all eleven parameters, then evaluates combinations of improving changes under the complete objective. It retains the original overnight cutoff.' if contract.get('schema')=='e5f_joint_local_recovery_v1' else 'The complete failure and status receipts accompany this attempt.'),
          '', '**Next discussion:** assess the complete fit and parameter restrictions below before selecting the presentation calibration. This finite search does not establish a global optimum. The empirical childbirth housing target has not been lowered.',
          '', '\\clearpage', '', '# Complete target fit', '',
          'Every gap is model minus target; the loss contribution is weight times squared gap. Provisional weight scales mean the total is an optimization objective, not a chi-squared statistic.','',
          '| Moment | Target | Model | Gap | Weight | Loss |','|---|---:|---:|---:|---:|---:|']
    for label,x in zip(LABELS,fit):
        text.append('| '+label+' | '+' | '.join(f"{float(x[k]):.6g}" if k=='weight' else f"{float(x[k]):.4f}" for k in ['target','model','gap','weight','loss_contribution'])+' |')
    text+=['', 'Completed fertility and childlessness refer to the cohort centered at age 42 in 2023. The first-birth rooms response compares matched 2019-to2023 branches, allowing later births in the treated branch. Parent room gaps use ages 30--55. The final wealth statistic is p90/p50 of wealth divided by annual income at ages 76--84.',
           '', '\\clearpage', '', '# Estimated parameters and restrictions', '',
           '| Parameter | Estimate | Lower | Upper | Near bound? |','|---|---:|---:|---:|---|']
    for x in params:
        if x['is_free_parameter']=='True':
            text.append('| '+PARAMETER_LABELS.get(x['parameter'],x['parameter'])+' | '+
                ' | '.join(f"{float(x[k]):.6f}" for k in ['value','lower_bound','upper_bound'])+' | '+('Yes' if x['near_bound']=='True' else 'No')+' |')
    text+=['','Near-bound flags use the saved rule: within 2% of the raw parameter range. They need not imply proximity in the logarithmic search coordinate.',
           '', '| Other recorded object | Value |','|---|---:|']
    for x in params:
        if x['is_free_parameter']!='True':text.append('| '+PARAMETER_LABELS.get(x['parameter'],x['parameter'])+f" | {float(x['value']):.6f} |")
    text+=['','Dated supply elasticity 0.63 and tenure dispersion 0.005 are externally fixed. Old initialization retains elasticity 1.75; old fertility is normalized to 2.1. The dated child-value change is estimated. The inherited household-entry/population closure is retained and remains a limitation for interpreting long forecasts.',
           '', '\\clearpage', '', '# Verification and evidence', '',
           f"The source, target, history, market, accounting, feasibility and population checks passed for the selected point. Its occupied adjacent-wealth value screen flags **{r['policy_array_diagnostic']['occupied_negative_steps']} decreases**; the exposed pre-choice mass share is {r['policy_array_diagnostic']['share_pre_choice_mass_at_negative_steps']:.4g}. Its reported-budget excess share is {r['budget_diagnostic']['budget_excess_share']:.4g}.",'',
           ('Both final repetitions exactly match the twelve fit rows, physical parameters, numeric historical entries and all seventeen standard PNGs.' if verified_repeat else 'The best point is not certified by two final exact repetitions. Inspect the failure/status receipt before using it.'),'',
           f"The local collector verified {verified} available original artifact hashes. Large dated checkpoints remain on Torch. All completed histories retain their standard seventeen-graph packet; the selected packet is collected locally. Visual review of this automatically generated PDF and the selected diagnostic plots remains pending.",'',
           f"[Experiment design and reproducibility]({OUT/'README.md'}). [Complete case ledger]({OUT/'search/all_cases.csv'}). [All target fits]({OUT/'search/all_target_fits.csv'}). [All parameter tables]({OUT/'search/all_parameters.csv'}). [Selected diagnostics]({case/'standard_diagnostics/summary.json'}).",'',
           'The earlier local Jacobian had a weak bequest direction and a nonsmooth old wealth-to-income quantile. A lower loss alone does not resolve identification, age-window measurement, old-age ownership saturation or population-closure questions. These remain separate from the target-and-equation-preserving search.']
    source=OUT/'MORNING_REVIEW.md';source.write_text('\n'.join(text)+'\n')
    pdf=ROOT/'output/pdf'/PDF_NAME
    command([sys.executable,str(ROOT/'code/model/tools/build_e5f_independent_audit_pdf.py'),
             '--source',str(source),'--output',str(pdf),'--heading','FERTILITY / OVERNIGHT CALIBRATION','--date','6 September 2026','--no-source-index'],timeout=180)
    render=ROOT/'tmp/pdfs'/pdf.stem;render.mkdir(parents=True,exist_ok=True)
    command(['pdftoppm','-scale-to','1400','-png',str(pdf),str(render/'page')],timeout=180)
    write(OUT/'automatic_report_receipt.json',{'status':'rendered_visual_review_pending','source':str(source),
          'pdf':str(pdf),'pdf_sha256':digest(pdf),'available_artifact_hashes_verified':verified,
          'exact_final_repeats':verified_repeat,'selected_case':str(case),'production_promoted':False})


def main():
    global REL, OUT, REMOTE, PDF_NAME
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--watch',action='store_true');parser.add_argument('--job-id')
    parser.add_argument('--recovery',action='store_true',help='Collect the separately preserved local-step recovery.')
    args=parser.parse_args()
    if args.recovery:
        REL=REL/'recovery_01';OUT=ROOT/REL
        REMOTE=Path('/scratch/td2248/projects/Fertility_Spring26_full_joint_overnight_20260906b')
        PDF_NAME='e5f_full_joint_recovery_review.pdf'
    deadline=read(OUT/'contract.json')['absolute_finish_epoch']+1800
    failures=0;last_count=-1
    while True:
        try:
            if args.job_id:
                queue=queue_state(args.job_id)
                if queue == 'PENDING':
                    write(OUT/'watch_status.json',{'status':'waiting_for_cluster_priority','epoch':time.time(),'job_id':args.job_id})
                    if not args.watch or time.time()>deadline:return
                    time.sleep(60)
                    continue
            raw=command(['ssh','-o','BatchMode=yes','torch',f'cat {REMOTE/REL}/search/search_state.json'],timeout=30)
            state=json.loads(raw);count=state['completed_histories']
            write(OUT/'watch_status.json',{'status':'watching','remote_state':state,'epoch':time.time()})
            terminal=state['status'] in ('complete_verified','failed','time_budget_exhausted_without_final_repeats')
            if count!=last_count or terminal:
                pull();last_count=count
            if terminal:
                build_report();return
            if args.job_id:
                queue=queue_state(args.job_id)
                if not queue:
                    pull()
                    write(OUT/'watch_status.json',{'status':'job_ended_without_terminal_report','epoch':time.time(),'remote_state':state})
                    return
            failures=0
        except Exception as error:
            failures+=1;write(OUT/'watch_status.json',{'status':'collection_error','error':str(error),'attempts':failures,'epoch':time.time()})
            if not args.watch:raise
        if not args.watch or time.time()>deadline:return
        # Finite collection process, not a recurring automation or model retry.
        time.sleep(60 if failures else 300)

if __name__=='__main__':
    import sys
    main()
