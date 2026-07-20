#!/usr/bin/env python3
"""Paired E1 saved-price parity and one-package strict-GE benchmark."""
from __future__ import annotations
import argparse, json, os, time
import numpy as np
from .e1_profile import E1_LOSS, E1_PRICE, e1_overrides, e1_target_system
os.environ.setdefault("NUMBA_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

def api(package):
    stem = "intergen_eqscale_seq" if package == "old" else "intergen_eqscale_seq_optimized"
    module = __import__(stem + ".solver", fromlist=["run_model_cp_dt"])
    cal = __import__(stem + ".calibration", fromlist=["extract_moments", "diagnostic_loss"])
    return module.run_model_cp_dt, cal.extract_moments, cal.diagnostic_loss
def run(package, mode):
    solve, extract, lossfn = api(package); o=e1_overrides(tight=True, optimized=package=="optimized")
    if mode == "fixed-price": o.update(solve_mode="pe", p_fixed=np.array([E1_PRICE]), w_fixed=np.array([1.]), entry_shares_fixed=np.array([1.]))
    t=time.perf_counter(); sol,P,p=solve(o,verbose=False); elapsed=time.perf_counter()-t; system=e1_target_system(); moments=extract(sol,P)
    return {"package":package,"mode":mode,"elapsed_seconds":elapsed,"price":np.asarray(p).tolist(),"loss":float(lossfn(moments,targets=system.targets_dict(),weights=system.weights_dict())),"canonical_loss":E1_LOSS,"market_residual":float(sol.best_max_abs_rel_excess),"mass":float(np.sum(sol.g)),"timings":sol.timings,"moments":{n:float(moments[n]) for n in system.moment_names}},sol
def main():
    pa=argparse.ArgumentParser();pa.add_argument("--package",choices=("old","optimized","both"),default="both");pa.add_argument("--mode",choices=("fixed-price","tight-ge"),default="fixed-price"); args=pa.parse_args()
    if args.package=="both" and args.mode!="fixed-price": raise ValueError("both only supports fixed-price")
    names=("old","optimized") if args.package=="both" else (args.package,); records={}; sols={}
    for n in names: records[n],sols[n]=run(n,args.mode)
    out={"runs":records}
    if len(sols)==2:
      arrays=("V","c_pol","hR_pol","bp_pol","tenure_choice","tenure_probs","loc_probs","fert_probs","fert_value","g")
      out["array_parity"]={n:{"exactly_equal":bool(np.array_equal(getattr(sols["old"],n),getattr(sols["optimized"],n))), "max_abs_difference":float(np.max(np.abs(getattr(sols["old"],n)-getattr(sols["optimized"],n))))} for n in arrays}
      out["speedup"]=records["old"]["elapsed_seconds"]/records["optimized"]["elapsed_seconds"]
      out["max_abs_target_moment_difference"]=max(abs(records["old"]["moments"][n]-records["optimized"]["moments"][n]) for n in records["old"]["moments"])
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__ == "__main__": main()
