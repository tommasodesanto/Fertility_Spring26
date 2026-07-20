#!/usr/bin/env python3
"""Checkpointed E1 promotion tasks for local smoke and one-CPU Torch runs."""
from __future__ import annotations
import argparse, json, os, time
from pathlib import Path
from typing import Any
import numpy as np
for key in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ.setdefault(key, "1")
from .benchmark_e1 import api
from .e1_profile import E1_PRICE, E1_THETA, e1_overrides, e1_target_system

def jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {str(k): jsonable(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)): return [jsonable(v) for v in value]
    return value

def evaluate(package: str, *, ge: bool, tight_ge: bool = False, price: float = E1_PRICE, changes: dict[str, Any] | None = None, p_init: float | None = None, smoke: bool = False) -> dict[str, Any]:
    solve, extract, loss_fn = api(package); overrides=e1_overrides(tight=ge and tight_ge, optimized=package == "optimized")
    if changes: overrides.update(changes)
    if smoke: overrides.update(Nb=30, max_iter_eq=4 if ge else overrides["max_iter_eq"], tol_eq=1e-3 if ge else overrides["tol_eq"])
    if p_init is not None: overrides["p_init_override"] = np.array([p_init])
    if not ge: overrides.update(solve_mode="pe", p_fixed=np.array([price]), w_fixed=np.array([1.0]), entry_shares_fixed=np.array([1.0]))
    started=time.perf_counter(); sol,P,p=solve(overrides,verbose=False); elapsed=time.perf_counter()-started
    moments=extract(sol,P); system=e1_target_system()
    timing=dict(sol.timings); refine=dict(timing.get("scalar_market_refine", {}))
    return {"package":package,"elapsed_seconds":elapsed,"price":np.asarray(p),"market_residual":float(sol.best_max_abs_rel_excess),"loss":float(loss_fn(moments,targets=system.targets_dict(),weights=system.weights_dict())),"moments":{n:float(moments[n]) for n in system.moment_names},"housing_demand":float(sol.housing_demand[0]),"housing_supply":float(sol.housing_supply[0]),"timings":timing,"directional_fallback_used":bool(refine.get("directional_fallback_used",False))}

def paired_fixed_price(price: float, smoke: bool) -> dict[str, Any]:
    old, optimized = evaluate("old",ge=False,price=price,smoke=smoke), evaluate("optimized",ge=False,price=price,smoke=smoke)
    old["excess_demand"] = old["housing_demand"] - old["housing_supply"]; optimized["excess_demand"] = optimized["housing_demand"] - optimized["housing_supply"]
    return {"old":old,"optimized":optimized,"excess_abs_difference":abs(old["excess_demand"]-optimized["excess_demand"]),"sign_equal":bool(np.signbit(old["excess_demand"])==np.signbit(optimized["excess_demand"]))}

def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("task",choices=("demand-map","root-case","strict-repeat","throughput")); parser.add_argument("--output",type=Path,required=True); parser.add_argument("--smoke",action="store_true"); args=parser.parse_args()
    if args.task == "demand-map":
        prices=np.linspace(E1_PRICE*.8,E1_PRICE*1.2,3 if args.smoke else 17); payload={"task":args.task,"smoke":args.smoke,"points":[{"price":float(p),**paired_fixed_price(float(p),args.smoke)} for p in prices]}
    elif args.task == "root-case":
        forced=evaluate("optimized",ge=True,tight_ge=True,p_init=E1_PRICE*1.30,changes={"scalar_market_refine_max_expand":0},smoke=args.smoke)
        if not forced["directional_fallback_used"]:
            raise AssertionError("forced wrong-direction root did not use bilateral fallback")
        payload={"task":args.task,"smoke":args.smoke,"ordinary":evaluate("optimized",ge=True,tight_ge=True,smoke=args.smoke),"forced_wrong_direction":forced,"forced_wrong_direction_exercised":True}
    elif args.task == "strict-repeat":
        first,second=evaluate("optimized",ge=True,tight_ge=True,smoke=args.smoke),evaluate("optimized",ge=True,tight_ge=True,smoke=args.smoke)
        payload={"task":args.task,"smoke":args.smoke,"first":first,"second":second,"price_bitwise_equal":np.array_equal(first["price"],second["price"]),"loss_bitwise_equal":first["loss"]==second["loss"],"moments_bitwise_equal":first["moments"]==second["moments"]}
    else:
        factors=np.linspace(.99,1.01,3 if args.smoke else 10); rows=[]
        for factor in factors:
            changes={"beta":E1_THETA["beta"]*float(factor)}
            rows.append({"factor":float(factor),"fixed_price":{"old":evaluate("old",ge=False,changes=changes,smoke=args.smoke),"optimized":evaluate("optimized",ge=False,changes=changes,smoke=args.smoke)},"ge_loose":{"old":evaluate("old",ge=True,changes=changes,smoke=args.smoke),"optimized":evaluate("optimized",ge=True,changes=changes,smoke=args.smoke)}})
        payload={"task":args.task,"smoke":args.smoke,"candidates":rows,"fixed_price_total_seconds":{package:sum(r["fixed_price"][package]["elapsed_seconds"] for r in rows) for package in ("old","optimized")},"ge_loose_total_seconds":{package:sum(r["ge_loose"][package]["elapsed_seconds"] for r in rows) for package in ("old","optimized")}}
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(jsonable(payload),indent=2,sort_keys=True)+"\n"); print(args.output)
if __name__ == "__main__": main()
