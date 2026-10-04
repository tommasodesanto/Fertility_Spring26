"""Re-normalize psi_child so total births equal chain-13 baseline, then read the gradient (fixed price)."""
import os, sys, json
sys.dont_write_bytecode=True
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import gradient_fixes as gf
from production.inputs import load_inputs, DEFAULT_PRICE, DEFAULT_PARAMETERS
from production.equilibrium import solve_at_price
import numpy as np
TARGET=0.11528810162653443; PSI0=float(DEFAULT_PARAMETERS["psi_child"])
def run_psi(psi,G0=0.0,Gn=0.0,kscale=1.0):
    orig=gf.load_inputs
    def patched(parameters=None,**kw):
        p=dict(parameters or {}); p["psi_child"]=psi; return orig(parameters=p,**kw)
    gf.load_inputs=patched
    try: return gf.run(G0=G0,Gn=Gn,kscale=kscale)
    finally: gf.load_inputs=orig
R=[]
for c in (dict(G0=0.75,Gn=0.3),dict(G0=0.75,Gn=0.6),dict(kscale=2.0)):
    lo,hi=0.001,PSI0
    for _ in range(25):
        mid=0.5*(lo+hi); r=run_psi(mid,**c)
        if r["births"]>TARGET: hi=mid
        else: lo=mid
        if abs(r["births"]/TARGET-1)<2e-3: break
    r["psi"]=mid; R.append(r); json.dump(R,open(os.path.join(gf.OUT,"gradient_fixes_renorm.json"),"w"),indent=1)
    f=lambda d:"/".join(f"{d[k]['ceb']:.2f}" for k in gf.GROUPS); g=lambda d:"/".join(f"{100*d[k]['childless']:.0f}%" for k in gf.GROUPS)
    print(f"{c} psi={mid:.4f} | births {r['births']:.4f} age25 {r['age25']:.3f} | CEB@26 {f(r['at26'])} | CEB@42-45 {f(r['at42_45'])} | childless@42-45 {g(r['at42_45'])} | transfers {r['transfer_pct_mean_period_earnings']:.2f}%",flush=True)
