#!/usr/bin/env python3
"""Execute on Torch: validate a pinned manifest and submit one bounded array."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess


def prepare(manifest, dependency=None):
    workers=int(manifest['max_workers']); hours=float(manifest['horizon_hours'])
    if not 1 <= workers <= 18 or not 0 < hours <= 12:
        raise ValueError('Explicit positive limits of at most18 workers/12hours required')
    pins=manifest['source_pins']
    if not pins:raise ValueError('Source pins required')
    for name,digest in pins.items():
        with Path(name).open('rb') as stream:
            if hashlib.file_digest(stream,'sha256').hexdigest()!=digest:
                raise ValueError('Source pin failed: '+name)
    stages=manifest['stages'];root=Path(manifest['runroot'])
    if not stages or not root.is_absolute():raise ValueError('Stages and absolute runroot required')
    names=[s['name'] for s in stages]
    if len(set(names))!=len(names) or any(not re.fullmatch(r'[A-Za-z0-9_-]+',n) for n in names):
        raise ValueError('Unique simple stage names required')
    seconds=max(int(s['seconds']) for s in stages)
    if any(not 0<int(s['seconds'])<=hours*3600 or s.get('cpus',1)!=1 for s in stages):
        raise ValueError('Every stage requires one CPU and a bounded runtime')
    env=dict(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
        NUMBA_NUM_THREADS='1',NUMBA_DISABLE_JIT='0',PYTHONUNBUFFERED='1',MPLBACKEND='Agg')
    env.update(manifest.get('environment',{}))
    if any(not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*',k) for k in env):raise ValueError('Invalid environment key')
    if any(env[k]!='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS')):
        raise ValueError('Nested numerical threading forbidden')
    lines=['#!/bin/bash','set -euo pipefail','module load anaconda3/2025.06']
    lines += [f'export {k}={shlex.quote(str(v))}' for k,v in env.items()]
    lines += ['case "$SLURM_ARRAY_TASK_ID" in']
    for i,stage in enumerate(stages):
        argv=stage['command']
        if not isinstance(argv,list) or not argv or any(not isinstance(v,str) for v in argv):raise ValueError('Command argv required')
        destination=root/stage['name']
        lines += [f'{i})',f'mkdir -p {shlex.quote(str(destination))}',
            shlex.join(argv)+f' >{shlex.quote(str(destination/"driver.log"))} 2>&1 ;;']
    lines += ['*) exit 64 ;;','esac']
    root.mkdir(parents=True,exist_ok=True);script=root/'final_night_array.sh'
    script.write_text('\n'.join(lines)+'\n')
    subprocess.run(['bash','-n',str(script)],check=True)
    wall=f'{seconds//3600:02d}:{seconds%3600//60:02d}:{seconds%60:02d}'
    command=['sbatch','--parsable','--account=torch_pr_570_general',
        f'--array=0-{len(stages)-1}%{workers}','--cpus-per-task=1',
        f'--mem={manifest.get("memory","24G")}',f'--time={wall}',
        '--job-name=e5f_final_history',f'--output={root}/slurm_%A_%a.out',f'--error={root}/slurm_%A_%a.err']
    if dependency:
        if not re.fullmatch(r'[0-9]+',dependency):raise ValueError('Numeric prerequisite job id required')
        command.append('--dependency=afterany:'+dependency)
    return command+[str(script)]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('manifest',type=Path)
    parser.add_argument('--submit',action='store_true');parser.add_argument('--dependency')
    args=parser.parse_args();manifest=json.loads(args.manifest.read_text())
    command=prepare(manifest,args.dependency)
    record=dict(command=command,manifest_sha256=hashlib.sha256(args.manifest.read_bytes()).hexdigest())
    if args.submit:
        record['job']=subprocess.check_output(command,text=True).strip().split(';')[0]
        (Path(manifest['runroot'])/'submission.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))

if __name__=='__main__':main()
