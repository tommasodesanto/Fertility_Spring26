"""Diagnostic-only replay: save original operator inputs/outputs on failure.

No numerical result, tolerance or exception is changed. Use the original frozen
scientific snapshot and a separate, explicitly pinned one-case replay plan.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ.setdefault(key,'1')
import argparse,sys,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'code/model'),str(ROOT/'code/model/tools')]
import numpy as np
import run_e5f_bounded_calibration_refinement as adapter
from intergen_eqscale_seq_optimized import two_shock_choice as operator

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--plan-sha256',required=True)
    ap.add_argument('--case-id',type=int,required=True);ap.add_argument('--capture-dir',type=Path,required=True)
    args=ap.parse_args();args.capture_dir.mkdir(parents=True,exist_ok=True)
    plan=adapter.load_plan(args.plan,args.plan_sha256)
    adapter.verify(__file__,plan["diagnostic_wrapper_sha256"])
    original=operator.choose
    def observed(*pos,**kw):
        try:return original(*pos,**kw)
        except Exception as error:
            tb=error.__traceback__;found=None
            while tb:
                if tb.tb_frame.f_code is original.__code__:found=tb.tb_frame.f_locals
                tb=tb.tb_next
            if found is not None:
                arrays={k:np.array(found[k],copy=True) for k in ('q','sf','v','p','occupied') if k in found}
                np.savez_compressed(args.capture_dir/'original_failure.npz',**arrays)
                report={'error':repr(error),'housing_scale':float(found['housing_scale']),
                    'tolerance':float(found['tolerance']),'wrapper_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
                if 'p' in arrays:
                    p=arrays['p'];mass=p.sum(axis=(1,2));occ=arrays['occupied']
                    bad=(p<0).any(axis=(1,2))|(p>1+10*report['tolerance']).any(axis=(1,2))|(occ&(abs(mass-1)>10*report['tolerance']))
                    report.update(min_probability=float(p.min()),max_probability=float(p.max()),
                        max_occupied_sum_error=float(abs(mass[occ]-1).max()),bad_indices=np.flatnonzero(bad).tolist())
                adapter.write_json(args.capture_dir/'failure_diagnostic.json',report)
            raise
    operator.choose=observed
    try:adapter.run_case(args)
    finally:operator.choose=original

if __name__=='__main__':main()
