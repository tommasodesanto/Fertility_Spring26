"""Check the analytic initial-pension marginal on a pinned saved full state.

Read-only distribution audit: balancing the pension on the old distribution
does not re-solve households or establish a new equilibrium.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import pickle
import sys
import numpy as np
from e5f_stationary_paygo import bind_initial_balanced_pension, certify_initial_pension


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract',type=Path,required=True)
    parser.add_argument('--contract-sha256',required=True)
    parser.add_argument('--source-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if digest(args.contract)!=args.contract_sha256: raise ValueError('Changed inherited audit contract')
    c=json.loads(args.contract.read_text())
    source=args.source_root.resolve()
    for name,pin in c['source_sha256'].items():
        p=(source/name).resolve()
        if not p.is_relative_to(source) or digest(p)!=pin: raise ValueError(f'Changed inherited source: {name}')
    if digest(c['normalized_checkpoint'])!=c['normalized_checkpoint_sha256']:
        raise ValueError('Changed inherited checkpoint')
    sys.path[:0]=[str(source/'code/model/tools'),str(source/'code/model')]
    with gzip.open(c['normalized_checkpoint'],'rb') as stream: old=pickle.load(stream)['old']
    if float(old.parameters.tau_pay)!=.179: raise ValueError('Inherited payroll tax changed')
    P, prediction=bind_initial_balanced_pension(old.parameters,payroll_tax=.179)
    certificate=certify_initial_pension(old.stationary_g_pre,P,
        marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
    packet=dict(status='passed_saved_distribution_marginal_audit',model_solves=0,
        equilibrium_certified=False, inherited_source_files_verified=len(c['source_sha256']),
        inherited_checkpoint_sha256=c['normalized_checkpoint_sha256'],
        inherited_pension_period=float(old.parameters.pension),
        analytic_balanced_pension_period=float(P.pension),
        prediction=prediction,certificate=certificate,
        interpretation='Analytic age/earnings law verified on old actual stationary heads; revised household equilibrium still required')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream: json.dump(packet,stream,indent=2,allow_nan=False)
    print(json.dumps(packet,allow_nan=False),flush=True)


if __name__=='__main__': main()
