"""Torch-only explicit readout watchers for the two additional 100-period starts."""
import hashlib
import json
import math
from pathlib import Path
import shlex
import subprocess
import time

B=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
template=json.loads((B/'readout_collector_corrected_seeded_v1_submission.json').read_text())
for group,science,helper in (
    ('corrected_cpu32_100','corrected_auto_v2','auto_v1'),
    ('corrected_extended_100','corrected_extended_seed_v1','extended_seed_v1')):
    receipt_path=B/f'readout_collector_{group}_submission.json'
    if receipt_path.exists():raise RuntimeError('Refusing duplicate submission: '+str(receipt_path))
    argv=template['reader_command'];argv=argv[:argv.index('--case-dir')]
    sm=B/f'history_manifest_{science}.json';output=B/f'readout_collection_{group}'
    for flag,value in [('--output',str(output)),('--scientific-manifest',str(sm)),
            ('--scientific-manifest-sha256',sha(sm)),('--helper-root',str(B/f'history_source_{helper}'))]:
        argv[argv.index(flag)+1]=value
    for case in ('A0','Aplus'):argv.extend(['--case-dir',str(B/f'histories_{group}'/f'{case}_100')])
    for name,digest in template['files'].items():
        if Path(name).parent==B/'readout_source':assert sha(name)==digest
    subprocess.run(argv+['--once'],check=True,capture_output=True,text=True,timeout=120)
    status=json.loads((output/'latest.json').read_text())
    assert all(row['status'] in ('PENDING','READY') for row in status['cases'])
    script=B/f'readout_collector_{group}.sbatch'
    env='\n'.join('export '+k+'=1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'))
    script.write_text('#!/bin/bash\nset -euo pipefail\n'+env+'\n'+shlex.join(argv)+'\n')
    subprocess.run(['bash','-n',str(script)],check=True)
    minutes=math.ceil((1789322400-time.time())/60);assert 0<minutes<720
    command=['sbatch','--parsable','--account=torch_pr_570_general','--cpus-per-task=1','--mem=12G',
        '--time='+str(minutes),'--job-name=e5f_readout_collect',
        f'--output={B}/readout_collect_{group}_%j.out',f'--error={B}/readout_collect_{group}_%j.err',str(script)]
    result=subprocess.run(command,check=True,capture_output=True,text=True)
    receipt=dict(job=result.stdout.strip().split(';')[0],command=command,reader_command=argv,
        files={str(p):sha(p) for p in [script,sm]+[Path(x) for x in template['files'] if Path(x).parent==B/'readout_source']},
        preflight_status=[row['status'] for row in status['cases']])
    receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(group=group,job=receipt['job'],preflight_status=receipt['preflight_status'])))
