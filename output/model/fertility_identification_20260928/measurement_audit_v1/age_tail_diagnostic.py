"""Saved-only age 25/26 and teenage first-birth diagnostic; run on Torch."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path

ROOT = Path('/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project')
BASE = ROOT / 'output/model/fertility_identification_20260928'
NIGHT = BASE / 'two_stream_overnight_v1'
OUT = BASE / 'measurement_audit_v1/age_tail_v1'
EMP = OUT / 'input/early_fertility_target.json'
NCHS = ROOT / 'code/data/nchs_natality_timing/first_birth_counts_year_age.csv'
CASES = {
    'Original selected': NIGHT / 'run_v1/one_birth/one_birth_024_gn1_0/case',
    'Two-birth selected': NIGHT / 'run_v1/two_birth/two_birth_024_gn1_0/case',
}
REF = BASE / 'resume_v1/selected_export/primary'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def model_at_age(accounting, age):
    starts = [int(x) for x in accounting['age_cell_start']]
    i = starts.index(18 + 4 * ((age - 18) // 4))
    start = starts[i]
    # Integer interview age A is the interval [A,A+1); interpolate at midpoint.
    alpha = (age + .5 - start) / 4
    pre = accounting['pre_parity_mass_by_age'][i]
    post = accounting['post_parity_mass_by_age'][i]
    assert len(pre) == len(post) == 4 and 0 <= alpha <= 1
    mass = [(1-alpha)*float(a)+alpha*float(b) for a,b in zip(pre,post)]
    total = math.fsum(mass)
    assert total > 0 and all(x >= -1e-12 for x in mass)
    share = [x/total for x in mass]
    mother = 1-share[0]
    count = math.fsum(n*s for n,s in enumerate(share))
    return dict(cell_start=start, post_weight=alpha, children_capped3=count,
                mother_share=mother, children_among_mothers_capped3=count/mother,
                share0=share[0], share1=share[1], share2=share[2], share3plus=share[3])


def main():
    assert os.environ.get('SLURM_JOB_ID','').isdigit() and os.uname().sysname == 'Linux'
    OUT.mkdir(parents=True, exist_ok=True)
    empirical = read_json(EMP)
    assert empirical['authoritative_candidate'] == '25'
    assert empirical['estimates']['25']['n'] == 1774
    inputs = {str(EMP): digest(EMP), str(NCHS): digest(NCHS)}
    model = {}
    for label, case in [('Reference',REF), *CASES.items()]:
        path = case / 'observers.json'
        if label == 'Reference':
            hashes = read_json(case/'artifact_hashes.json')
        else:
            success = read_json(case.parent/'SUCCESS.json')
            assert success['status'] == 'passed' and success['candidate_id'] == case.parent.name
            hashes = {Path(x['path']).relative_to(Path(success['case_path'])).as_posix(): x['sha256'] for x in success['artifacts']}
        assert digest(path) == hashes['observers.json'], label
        inputs[str(path)] = digest(path)
        observer = read_json(path)['fertility']['uniform_birth_time']
        account = observer['accounting']
        model[label] = {str(age): model_at_age(account, age) for age in (25,26)}
    with NCHS.open(newline='') as f:
        rows = list(csv.DictReader(f))
    counts = {age:math.fsum(int(r['n_first_births']) for r in rows if 2003 <= int(r['year']) <= 2006 and int(r['age']) == age) for age in range(12,50)}
    total = math.fsum(counts.values())
    assert abs(total-6611269) < 1, total
    teen = {str(cut):dict(first_births_before_age=math.fsum(n for age,n in counts.items() if age < cut),
                          share_of_period_first_births=math.fsum(n for age,n in counts.items() if age < cut)/total)
            for cut in (18,20)}
    # The earlier published pre-18 figure is an independent replication gate.
    assert abs(teen['18']['share_of_period_first_births'] - .07731359894749404) < 1e-12
    data = {str(age):dict(n=empirical['estimates'][str(age)]['n'],
                           children_capped3=empirical['estimates'][str(age)]['mean_children_ever_born_capped3'],
                           mother_share=empirical['estimates'][str(age)]['share_with_any_birth'],
                           children_among_mothers_capped3=empirical['estimates'][str(age)]['mean_children_ever_born_capped3']/empirical['estimates'][str(age)]['share_with_any_birth'],
                           children_capped3_bootstrap_se=empirical['estimates'][str(age)]['capped3_bootstrap_se'])
            for age in (25,26)}
    table=[]
    for age in (25,26):
        for label in ('Reference',*CASES):
            m=model[label][str(age)]; d=data[str(age)]
            table.append(dict(age=age,series=label,model_cell_start=m['cell_start'],model_post_weight=m['post_weight'],
                              data_n=d['n'],data_children_capped3=d['children_capped3'],model_children_capped3=m['children_capped3'],children_gap=m['children_capped3']-d['children_capped3'],
                              data_mother_share=d['mother_share'],model_mother_share=m['mother_share'],mother_gap=m['mother_share']-d['mother_share'],
                              data_children_among_mothers=d['children_among_mothers_capped3'],model_children_among_mothers=m['children_among_mothers_capped3'],
                              conditional_gap=m['children_among_mothers_capped3']-d['children_among_mothers_capped3']))
    with (OUT/'age25_26_fit.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(table[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(table)
    result=dict(status='PASS',slurm_job_id=os.environ['SLURM_JOB_ID'],model_solves=0,model_imports=0,checkpoint_reads=0,
                model=model,data=data,nchs_period_years=[2003,2004,2005,2006],nchs_first_births_ages12_49=total,
                nchs_teen_first_births=teen,inputs_sha256=inputs,
                definitions={'CPS':'Women age 25 or 26 at June 2004/2006 interview; FREVER 0–20, positive FRSUPPWT; pooled weighted mean min(FREVER,3). Bootstrap SE is year-stratified person bootstrap, not survey design.',
                             'model':'Uniform-birth-time linear pre/post state interpolation at midpoint of completed-integer-age interview interval. Age 25 uses 0.875 post weight in 22–26 cell; age 26 uses 0.125 post weight in 26–30 cell.',
                             'NCHS':'Pooled 2003–2006 period first births among mothers ages 12–49. Teen shares divide first births, not women or children per woman. Different population/cohort from CPS.'},
                limitations=['No linked CPS first-birth histories for the age-25/26 women, so missing pre-entry births cannot be identified from the NCHS period share.',
                             'Age 26 diagnostic changes both model and data interview age; no target change or adoption.',
                             'NCHS first-birth shares do not measure second or third births by 25.'])
    (OUT/'age_tail_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',job_id=os.environ['SLURM_JOB_ID'],table=table,teen=teen)))


if __name__ == '__main__':
    main()
