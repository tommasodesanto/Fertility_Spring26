#!/usr/bin/env python3
"""Six fixed-price E1 perturbation parity cases."""
from __future__ import annotations
import json, time
import numpy as np
from .benchmark_e1 import api
from .e1_profile import E1_PRICE, e1_overrides, e1_target_system
CASES=(("beta_plus_1pct",{"beta":0.9522043265603773 * 1.01}), ("alpha_cons",{"alpha_cons":.7006}), ("gamma_e",{"gamma_e":.3869}), ("delta_alpha",{"delta_alpha":.0427}), ("kappa_fert",{"kappa_fert":10.877}), ("H0",{"H0":7.4808}))
def solve(which, changes):
 s,e,_=api(which); o=e1_overrides(optimized=which=="optimized"); o.update(changes,solve_mode="pe",p_fixed=np.array([E1_PRICE]),w_fixed=np.array([1.]),entry_shares_fixed=np.array([1.])); t=time.perf_counter();x,P,_=s(o,verbose=False);return x,e(x,P),time.perf_counter()-t
def main():
 rows=[]; system=e1_target_system()
 for name,change in CASES:
  a,am,ta=solve("old",change);b,bm,tb=solve("optimized",change)
  rows.append({"case":name,"changes":change,"old_seconds":ta,"optimized_seconds":tb,"speedup":ta/tb,"max_abs_moment_difference":max(abs(am[k]-bm[k]) for k in system.moment_names),"max_abs_distribution_difference":float(np.max(np.abs(a.g-b.g)))})
 print(json.dumps({"cases":rows},indent=2))
if __name__=="__main__": main()
