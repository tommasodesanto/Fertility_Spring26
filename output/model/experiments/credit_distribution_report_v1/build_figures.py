"""Rebuild saved-distribution tables and paired figures; never solves the model."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "code/model"))
sys.path.insert(0, str(ROOT / "code/model/tools"))
from production.storage import load_case  # noqa: E402
from model_policy_tools import aggregate_solution  # noqa: E402

REF = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/credit_mechanism/credit_relaxation"
LOW = ROOT / "output/model/experiments/low_productivity_credit_v1"
CASES = {
    "reference_phi08": (REF / "phi_080", "reference", "phi=0.8", False),
    "reference_phi10": (REF / "phi_100", "reference", "phi=1.0", False),
    "low_phi08": (LOW / "phi_08", "low", "phi=0.8", True),
    "low_phi10": (LOW / "phi_10", "low", "phi=1.0", True),
}
COLORS = {"phi=0.8": "#1f5a99", "phi=1.0": "#d47b24"}
COUNTS = ["0", "1", "2", "3+"]
ROOMS = [2, 4, 6, 8, 10]


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def weighted_quantile(x, w, q):
    x, w = np.asarray(x), np.asarray(w)
    keep = w > 0
    if not np.any(keep):
        return None
    ix = np.argsort(x[keep], kind="stable")
    sx, sw = x[keep][ix], w[keep][ix]
    return float(sx[np.searchsorted(np.cumsum(sw), q * sw.sum(), side="left")])


def load(case_path, native):
    if native:
        result, _ = load_case(case_path)
        return result.solution, result.P, case_path / "native_result.npz", case_path / "metadata.json"
    with np.load(case_path / "solution_arrays.npz", allow_pickle=False) as z:
        needed = ("g", "g_beginning_distribution", "g_stay_distribution", "b_grid", "p_eq", "hR_pol",
                  "c_pol", "c_pol_stay", "bp_pol", "bp_pol_stay", "V")
        sol = SimpleNamespace(**{k: z[k].copy() for k in needed})
    par = json.loads((case_path / "executed_P.json").read_text())
    return sol, SimpleNamespace(**par), case_path / "solution_arrays.npz", case_path / "executed_P.json"


def positive_ratio(n, d):
    return None if d <= 0 else float(n / d)


def policy_distribution(sol, g, stay, field, field_stay, *, bins=160):
    buyer = np.maximum(g - stay, 0)
    v1, v2 = np.asarray(getattr(sol, field)), np.asarray(getattr(sol, field_stay))
    a, b = buyer > 0, stay > 0
    values = np.concatenate((v1[a], v2[b]))
    weights = np.concatenate((buyer[a], stay[b]))
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(weights)):
        raise ValueError("nonfinite occupied policy distribution")
    lo, hi = np.quantile(values, [0, 1])
    if hi <= lo:
        edges = np.array([lo - .5, lo + .5])
    else:
        edges = np.linspace(lo, hi, bins + 1)
    hist, edges = np.histogram(values, bins=edges, weights=weights)
    rows = [{"left": float(edges[i]), "right": float(edges[i+1]), "mass": float(hist[i]),
             "probability": float(hist[i]/weights.sum()), "cdf": float(hist[:i+1].sum()/weights.sum())}
            for i in range(len(hist))]
    order=np.argsort(values,kind="stable")
    sorted_values=values[order];cum=np.cumsum(weights[order]);cum/=cum[-1]
    probs=np.linspace(0,1,1001)
    sampled=sorted_values[np.searchsorted(cum,probs,side="left").clip(max=len(cum)-1)]
    ecdf=[{"value":float(v),"cdf":float(q)} for v,q in zip(sampled,probs)]
    quant = {f"q{int(q*100):02d}": float(sorted_values[np.searchsorted(cum,q,side="left")])
             for q in [.01, .10, .25, .50, .75, .90, .99]}
    quant["mean"] = float(np.dot(values, weights)/weights.sum())
    return rows, ecdf, quant


def extract(case_id, path, population, label, native):
    sol, P, source, parameters = load(path, native)
    g = np.asarray(sol.g); beginning = np.asarray(sol.g_beginning_distribution)
    stay = np.asarray(sol.g_stay_distribution); b = np.asarray(sol.b_grid)
    if g.ndim != 7 or g.shape != beginning.shape or g.shape != stay.shape or len(b) != g.shape[0]:
        raise ValueError(f"distribution shape error: {case_id}")
    if any(not np.all(np.isfinite(x)) or np.min(x) < -1e-12 for x in (g, beginning, stay)):
        raise ValueError(f"distribution mass invalid: {case_id}")
    if np.max(stay-g) > 1e-10 or stay[:, 0].sum() > 1e-10:
        raise ValueError(f"stayer mass invalid: {case_id}")
    if abs(g.sum()-1) > 1e-8 or abs(beginning.sum()-g.sum()) > 1e-10:
        raise ValueError(f"total mass invalid: {case_id}")
    age_mass = g.sum(axis=(0,1,2,4,5,6)); begin_age = beginning.sum(axis=(0,1,2,4,5,6))
    if np.max(np.abs(age_mass-begin_age)) > 1e-10:
        raise ValueError(f"age mass disagreement: {case_id}")
    ages = np.arange(g.shape[3])*float(P.period_years)+float(P.age_start)
    own = g.copy(); own[:, 0] = 0
    n_age = g.sum(axis=(0,1,2,4,6)); m_age = g.sum(axis=(0,1,2,4,5))
    own_n_age = own.sum(axis=(0,1,2,4,6)); own_m_age = own.sum(axis=(0,1,2,4,5))
    count_rows=[]; age_rows=[]
    for j, age in enumerate(ages):
        rent_age=g[:,0,:,j]
        rent_policy=np.asarray(sol.hR_pol)[:,0,:,j]
        rent_mask=rent_age>0
        room_total=float(np.dot(np.asarray(P.H_own,dtype=float),g[:,1:,:,j].sum(axis=(0,2,3,4,5))))
        room_total+=float(np.dot(rent_age[rent_mask],rent_policy[rent_mask]))
        age_rows.append({"case":case_id,"age_start":int(age),"age_end":int(age+P.period_years-1),
                         "mass":float(age_mass[j]),"ownership":positive_ratio(own[:,:,:,j].sum(),age_mass[j]),
                         "mean_rooms":positive_ratio(room_total,age_mass[j]),
                         "mean_children_ever_born_capped":positive_ratio(np.dot(np.arange(g.shape[5]),n_age[j]),age_mass[j]),
                         "mean_children_at_home":positive_ratio(np.dot(np.arange(g.shape[6]),m_age[j]),age_mass[j])})
        for k in range(g.shape[5]):
            count_rows.append({"case":case_id,"age_start":int(age),"kind":"ever_born","count":COUNTS[k],
                               "cell_mass":float(n_age[j,k]),"age_mass":float(age_mass[j]),
                               "share":positive_ratio(n_age[j,k],age_mass[j]),"cdf":positive_ratio(n_age[j,:k+1].sum(),age_mass[j]),
                               "ownership":positive_ratio(own_n_age[j,k],n_age[j,k])})
        for k in range(g.shape[6]):
            count_rows.append({"case":case_id,"age_start":int(age),"kind":"at_home","count":str(k),
                               "cell_mass":float(m_age[j,k]),"age_mass":float(age_mass[j]),
                               "share":positive_ratio(m_age[j,k],age_mass[j]),"cdf":positive_ratio(m_age[j,:k+1].sum(),age_mass[j]),
                               "ownership":positive_ratio(own_m_age[j,k],m_age[j,k])})
    count_path=HERE/"tables"/f"{case_id}_age_children.csv"
    age_path=HERE/"tables"/f"{case_id}_age_summary.csv"
    write_csv(count_path,list(count_rows[0]),count_rows);write_csv(age_path,list(age_rows[0]),age_rows)
    # Income state is a numerical productivity state. Suppress empty denominators.
    z_mass=g.sum(axis=(0,1,2,3,5,6)); z_own=own.sum(axis=(0,1,2,3,5,6))
    z_rows=[{"case":case_id,"state":int(z),"mass":float(z_mass[z]),"ownership":positive_ratio(z_own[z],z_mass[z])}
            for z in range(len(z_mass))]
    z_path=HERE/"tables"/f"{case_id}_productivity.csv";write_csv(z_path,list(z_rows[0]),z_rows)
    # Realized owned rooms are fixed by the occupied tenure branch.
    ten_mass=g.sum(axis=(0,2,3,4,5,6)); owner_mass=float(ten_mass[1:].sum())
    room_rows=[{"case":case_id,"rooms":ROOMS[t-1],"owner_mass":float(ten_mass[t]),
                "share_among_owners":positive_ratio(ten_mass[t],owner_mass)} for t in range(1,len(ten_mass))]
    room_path=HERE/"tables"/f"{case_id}_owned_rooms.csv";write_csv(room_path,list(room_rows[0]),room_rows)
    # Renter services are conditional on realized renter mass; check feasible occupied policies separately.
    renter=g[:,0]; rent_values=np.asarray(sol.hR_pol)[:,0]; mask=renter>0
    parent_floor=float(getattr(P,"hbar_first_child_jump",0.0))
    parent_mask=renter[...,1:]>0
    renter_ok=bool(np.all(np.isfinite(rent_values[mask])) and np.all(rent_values[mask]>0) and
                   np.all(rent_values[mask]<=6.0+1e-8) and
                   np.all(rent_values[...,1:][parent_mask]>=parent_floor-1e-8) and
                   np.all(np.isfinite(np.asarray(sol.V)[:,0][mask])) and
                   np.all(np.asarray(sol.V)[:,0][mask]>-1e9))
    renter_path=None; renter_quant=None
    if renter_ok:
        vals=rent_values[mask]; weights=renter[mask]
        edges=np.linspace(0,6,121); hist,_=np.histogram(vals,bins=edges,weights=weights)
        renter_rows=[{"case":case_id,"left":float(edges[i]),"right":float(edges[i+1]),
                      "mass":float(hist[i]),"share_among_renters":float(hist[i]/weights.sum())}
                     for i in range(len(hist))]
        renter_path=HERE/"tables"/f"{case_id}_rental_rooms.csv";write_csv(renter_path,list(renter_rows[0]),renter_rows)
        renter_quant={"mean":float(np.dot(vals,weights)/weights.sum()),"median":weighted_quantile(vals,weights,.5)}
    wealth=beginning.sum(axis=(1,2,3,4,5,6)); wealth_age=beginning.sum(axis=(1,2,4,5,6))
    wealth_rows=[{"case":case_id,"asset_node":float(b[i]),"mass":float(wealth[i]),
                  "probability":float(wealth[i]/wealth.sum()),"cdf":float(wealth[:i+1].sum()/wealth.sum())}
                 for i in range(len(b))]
    wealth_path=HERE/"tables"/f"{case_id}_beginning_assets.csv";write_csv(wealth_path,list(wealth_rows[0]),wealth_rows)
    wealth_age_rows=[]
    for j,age in enumerate(ages):
        w=wealth_age[:,j]
        wealth_age_rows.append({"case":case_id,"age_start":int(age),"mass":float(w.sum()),
             "mean":float(np.dot(b,w)/w.sum()),**{f"q{int(q*100):02d}":weighted_quantile(b,w,q) for q in [.1,.25,.5,.75,.9]}})
    wealth_age_path=HERE/"tables"/f"{case_id}_assets_by_age.csv";write_csv(wealth_age_path,list(wealth_age_rows[0]),wealth_age_rows)
    # Beginning tenure is inherited; owner housing has gross value price*rooms.
    # This is a pre-transaction accounting distribution, without selling costs.
    price=float(np.asarray(sol.p_eq).reshape(-1)[0]) if hasattr(sol,"p_eq") else float(json.loads((path/"metadata.json").read_text())["price"])
    housing_values=np.r_[0.0,np.asarray(P.H_own,dtype=float)*price]
    net_values=(b[:,None]+housing_values[None,:]).reshape(-1)
    net_weights_age=beginning.sum(axis=(2,4,5,6)).reshape(-1,len(ages))
    uniq,inv=np.unique(net_values,return_inverse=True)
    net_mass=np.bincount(inv,weights=net_weights_age.sum(axis=1),minlength=len(uniq))
    net_rows=[{"case":case_id,"gross_net_worth":float(v),"mass":float(net_mass[i]),
               "probability":float(net_mass[i]/net_mass.sum()),
               "cdf":float(net_mass[:i+1].sum()/net_mass.sum())} for i,v in enumerate(uniq)]
    net_path=HERE/"tables"/f"{case_id}_beginning_gross_net_worth.csv"
    write_csv(net_path,list(net_rows[0]),net_rows)
    by_ten=beginning.sum(axis=(2,3,4,5,6))
    net_ten_rows=[]
    for group,group_weight in (("inherited_renter",np.column_stack([by_ten[:,0],np.zeros_like(by_ten[:,1:])]).reshape(-1)),
                               ("inherited_owner",np.column_stack([np.zeros_like(by_ten[:,0]),by_ten[:,1:]]).reshape(-1))):
        gm=np.bincount(inv,weights=group_weight,minlength=len(uniq))
        if gm.sum()<=0:continue
        gcdf=np.cumsum(gm)/gm.sum()
        net_ten_rows.extend({"case":case_id,"group":group,"gross_net_worth":float(v),"mass":float(gm[i]),
                             "conditional_cdf":float(gcdf[i])} for i,v in enumerate(uniq))
    net_ten_path=HERE/"tables"/f"{case_id}_gross_net_worth_by_inherited_tenure.csv"
    write_csv(net_ten_path,list(net_ten_rows[0]),net_ten_rows)
    net_age_rows=[]
    for j,age in enumerate(ages):
        w=net_weights_age[:,j]
        net_age_rows.append({"case":case_id,"age_start":int(age),"mass":float(w.sum()),
            "mean":float(np.dot(net_values,w)/w.sum()),
            **{f"q{int(q*100):02d}":weighted_quantile(net_values,w,q) for q in [.1,.25,.5,.75,.9]}})
    net_age_path=HERE/"tables"/f"{case_id}_gross_net_worth_by_age.csv"
    write_csv(net_age_path,list(net_age_rows[0]),net_age_rows)
    policy_status="validated";policy_error=None;policy_tables=[];policy_summary={}
    try:
        agg=aggregate_solution(sol,houses=P.H_own,age_start=int(P.age_start),period_years=int(P.period_years))
        for name,fields in (("consumption",("c_pol","c_pol_stay")),("next_assets",("bp_pol","bp_pol_stay"))):
            rows,ecdf,stats=policy_distribution(sol,g,stay,*fields)
            p=HERE/"tables"/f"{case_id}_{name}_histogram.csv";write_csv(p,list(rows[0]),rows)
            ep=HERE/"tables"/f"{case_id}_{name}_ecdf.csv";write_csv(ep,list(ecdf[0]),ecdf)
            policy_tables.extend([str(p),str(ep)]);policy_summary[name]=stats
        policy_summary["aggregate_overall"]=agg["overall"]
        policy_summary["aggregate_by_age"]=agg["by_age"]
    except ValueError as e:
        policy_status="unavailable_native_reporter_gate";policy_error=str(e)
        if not case_id.endswith("phi10"):
            raise
    checks={"realized_mass":float(g.sum()),"beginning_mass":float(beginning.sum()),
            "max_age_mass_gap":float(np.max(np.abs(age_mass-begin_age))),
            "min_g":float(g.min()),"min_beginning":float(beginning.min()),
            "max_stayer_excess":float(np.max(stay-g)),
            "max_count_share_sum_gap":float(np.max(np.abs(n_age.sum(axis=1)/age_mass-1))),
            "min_age_mass":float(age_mass.min()),"rent_policy_occupied_feasible":renter_ok}
    buyer_mass=np.maximum(g-stay,0.0)
    bad_buyer=(buyer_mass>0)&((~np.isfinite(np.asarray(sol.V)))|(np.asarray(sol.V)<=-1e9))
    bad_index=np.argwhere(bad_buyer)
    checks["infeasible_occupied_buyer_cells"]=int(len(bad_index))
    checks["infeasible_occupied_buyer_mass"]=float(buyer_mass[bad_buyer].sum())
    checks["first_infeasible_buyer_index"]=bad_index[0].tolist() if len(bad_index) else None
    checks["first_infeasible_buyer_value"]=float(np.asarray(sol.V)[tuple(bad_index[0])]) if len(bad_index) else None
    d={"case":case_id,"population":population,"label":label,"source_path":str(source),
       "source_sha256":digest(source),"parameter_path":str(parameters),"parameter_sha256":digest(parameters),
       "price":price,
       "checks":checks,"policy_status":policy_status,"policy_error":policy_error,
       "policy_summary":policy_summary,"renter_quantiles":renter_quant,
       "tables":[str(p) for p in (count_path,age_path,z_path,room_path,renter_path,wealth_path,wealth_age_path,net_path,net_age_path,net_ten_path) if p is not None]+policy_tables,
       "standard_plot_paths":[str(p) for p in sorted((path/"standard_diagnostics").glob("*.png"))]}
    assert len(d["standard_plot_paths"])==17
    return d


def read_csv(path):
    with open(path,newline="") as f:return list(csv.DictReader(f))


def rows(case, suffix):
    return read_csv(HERE/"tables"/f"{case}_{suffix}.csv")


def style(ax, title, xlabel="", ylabel=""):
    ax.set_title(title,loc="left",weight="bold",fontsize=10)
    ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
    ax.grid(alpha=.18);ax.spines[["top","right"]].set_visible(False)


def savefig(name):
    plt.tight_layout(pad=1.8)
    p=HERE/"figures"/f"{name}.png";p.parent.mkdir(parents=True,exist_ok=True)
    plt.savefig(p,dpi=185,bbox_inches="tight",facecolor="white");plt.close()
    return str(p)


def make_figures(pair, title, valid_policy):
    cases=[f"{pair}_phi08",f"{pair}_phi10"]
    labels=["phi=0.8","phi=1.0"]
    figs=[]
    # Both children dimensions are stocks, with age-specific shares.
    fig,axs=plt.subplots(2,2,figsize=(10.4,7.4))
    for case,label in zip(cases,labels):
        r=rows(case,"age_children")
        for kind, row in (("ever_born",0),("at_home",1)):
            for k in range(4):
                sub=[x for x in r if x["kind"]==kind and x["count"]==(COUNTS[k] if row==0 else str(k))]
                x=[int(x["age_start"]) for x in sub]
                axs[row,0].plot(x,[float(v["share"]) for v in sub],color=COLORS[label],ls=["-","--",":","-."][k],lw=1.6,label=f"{label}, {COUNTS[k] if row==0 else k}")
                if k<3: axs[row,1].plot(x,[float(v["cdf"]) for v in sub],color=COLORS[label],ls=["-","--",":"][k],lw=1.6,label=f"{label}, <= {k}")
    for ax,t in zip(axs.flat,["Children ever born: share","Children ever born: cumulative share","Children at home: share","Children at home: cumulative share"]):style(ax,t,"Age-cell start","Share")
    for ax in axs.flat:ax.set_ylim(-.02,1.03)
    axs[0,0].legend(fontsize=7,ncol=2);axs[0,1].legend(fontsize=7,ncol=2)
    fig.suptitle(title+" | Children",fontsize=14,weight="bold")
    figs.append((savefig(pair+"_01_children"),"Children ever born and currently at home by four-year age cell. Each line divides by that case's own age-cell mass; 3+ is the top ever-born state."))
    fig,axs=plt.subplots(1,3,figsize=(12.2,3.7))
    for case,label in zip(cases,labels):
        ar=rows(case,"age_summary")
        axs[0].plot([int(x["age_start"]) for x in ar],[float(x["ownership"]) for x in ar],color=COLORS[label],lw=2,label=label)
        cr=rows(case,"age_children")
        for k in range(4):
            v=[x for x in cr if x["kind"]=="ever_born" and x["count"]==COUNTS[k]]
            mass=sum(float(x["cell_mass"]) for x in v)
            numerator=sum(float(x["ownership"])*float(x["cell_mass"]) for x in v if x["ownership"])
            axs[1].plot(k,numerator/mass if mass>0 else np.nan,marker="o",color=COLORS[label])
        z=rows(case,"productivity")
        zz=[x for x in z if float(x["mass"])>1e-12]
        axs[2].plot([int(x["state"]) for x in zz],[float(x["ownership"]) for x in zz],marker="o",color=COLORS[label],label=label)
    style(axs[0],"Ownership by age","Age-cell start","Owner share")
    style(axs[1],"Ownership by children ever born","Children ever born","Owner share")
    style(axs[2],"Ownership by occupied productivity state","Numerical state index","Owner share")
    axs[1].set_xticks(range(4),COUNTS)
    for ax in axs:ax.set_ylim(-.02,1.02)
    axs[0].legend(fontsize=8);axs[2].legend(fontsize=8)
    fig.suptitle(title+" | Ownership",fontsize=14,weight="bold")
    figs.append((savefig(pair+"_02_ownership"),"Ownership among all households by age, within child-count groups, and within occupied productivity states. Empty groups have no plotted rate."))
    fig,axs=plt.subplots(1,3,figsize=(12.2,3.7))
    width=.34
    for i,(case,label) in enumerate(zip(cases,labels)):
        rr=rows(case,"owned_rooms");axs[0].bar(np.arange(5)+(i-.5)*width,[float(x["share_among_owners"]) for x in rr],width,color=COLORS[label],label=label)
        rp=HERE/"tables"/f"{case}_rental_rooms.csv"
        if rp.exists():
            rt=read_csv(rp);axs[1].plot([(float(x["left"])+float(x["right"]))/2 for x in rt],np.cumsum([float(x["share_among_renters"]) for x in rt]),color=COLORS[label],lw=2,label=label)
        a=rows(case,"age_summary");axs[2].plot([int(x["age_start"]) for x in a],[float(x["mean_rooms"]) for x in a],color=COLORS[label],lw=2,label=label)
    axs[0].set_xticks(range(5),ROOMS);style(axs[0],"Owned room sizes","Rooms","Share of owners")
    style(axs[1],"Rental rooms: cumulative share","Rooms","Share of renters")
    style(axs[2],"Mean occupied rooms over the lifecycle","Age-cell start","Rooms")
    axs[0].legend(fontsize=8);axs[1].legend(fontsize=8)
    fig.suptitle(title+" | Housing",fontsize=14,weight="bold")
    figs.append((savefig(pair+"_03_housing"),"Owned rooms use realized owner tenure. Rental rooms use the saved conditional renter housing policy only on occupied renter mass after a finite and physical-range check."))
    fig,axs=plt.subplots(1,3,figsize=(12.2,3.7))
    allw=[rows(c,"beginning_assets") for c in cases]
    central=[]
    for r in allw:
        vals=[float(x["asset_node"]) for x in r];cdf=[float(x["cdf"]) for x in r]
        central += [vals[np.searchsorted(cdf,q)] for q in [.01,.99]]
    lo,hi=min(central[::2]),max(central[1::2])
    for case,label,r in zip(cases,labels,allw):
        x=[float(v["asset_node"]) for v in r];y=[float(v["cdf"]) for v in r]
        axs[0].step(x,y,where="post",color=COLORS[label],lw=2,label=label)
        axs[1].step(x,y,where="post",color=COLORS[label],lw=2,label=label)
        ar=rows(case,"assets_by_age")
        ages=[int(v["age_start"]) for v in ar]
        axs[2].plot(ages,[float(v["q50"]) for v in ar],color=COLORS[label],lw=2,label=label)
        axs[2].fill_between(ages,[float(v["q25"]) for v in ar],[float(v["q75"]) for v in ar],color=COLORS[label],alpha=.14)
    axs[0].set_xlim(lo,hi);style(axs[0],"Beginning assets: central 1–99%","Financial wealth","Cumulative share")
    style(axs[1],"Beginning assets: full grid","Financial wealth","Cumulative share")
    style(axs[2],"Median assets and middle half by age","Age-cell start","Financial wealth")
    axs[0].legend(fontsize=8)
    fig.suptitle(title+" | Beginning financial wealth",fontsize=14,weight="bold")
    figs.append((savefig(pair+"_04_wealth"),"Financial wealth is beginning-of-period net financial assets using the post-fertility, pre-tenure distribution. Units retain the original annual gross-earnings normalization in both populations."))
    fig,axs=plt.subplots(1,3,figsize=(12.2,3.8))
    net_central=[];group_central=[]
    for case,label in zip(cases,labels):
        nw=rows(case,"beginning_gross_net_worth")
        na=rows(case,"gross_net_worth_by_age")
        nt=rows(case,"gross_net_worth_by_inherited_tenure")
        vx=[float(v["gross_net_worth"]) for v in nw];vy=[float(v["cdf"]) for v in nw]
        net_central.append((vx[np.searchsorted(vy,.01)],vx[np.searchsorted(vy,.99)]))
        axs[0].step(vx,vy,where="post",color=COLORS[label],lw=1.8,label=label)
        for group,linestyle in (("inherited_renter","--"),("inherited_owner","-")):
            r=[v for v in nt if v["group"]==group]
            tx=[float(v["gross_net_worth"]) for v in r];ty=[float(v["conditional_cdf"]) for v in r]
            group_central.append((tx[np.searchsorted(ty,.01)],tx[np.searchsorted(ty,.99)]))
            axs[1].step(tx,ty,where="post",color=COLORS[label],lw=1.5,ls=linestyle,
                        label=f"{label}, {'renter' if group=='inherited_renter' else 'owner'}")
        axs[2].plot([int(v["age_start"]) for v in na],[float(v["q50"]) for v in na],color=COLORS[label],lw=2,label=label)
    axs[0].set_xlim(min(x[0] for x in net_central),max(x[1] for x in net_central))
    axs[1].set_xlim(min(x[0] for x in group_central),max(x[1] for x in group_central))
    style(axs[0],"Gross net worth: central 1–99%","Gross net worth","Cumulative share")
    style(axs[1],"By inherited tenure: central 1–99%","Gross net worth","Conditional cumulative share")
    style(axs[2],"Median gross net worth by age","Age-cell start","Gross net worth")
    axs[0].legend(fontsize=8);axs[1].legend(fontsize=7)
    fig.suptitle(title+" | Gross net worth",fontsize=14,weight="bold")
    figs.append((savefig(pair+"_05_net_worth"),"Pre-transaction gross net worth equals beginning net financial assets plus held house price times inherited owned rooms; renters have zero owned housing. CDF panels display each group's central 1–99% range without renormalizing probabilities; full tails are in the companion CSV. It excludes selling costs and is not an empirical wealth-target observer."))
    if valid_policy:
        fig,axs=plt.subplots(2,2,figsize=(10.4,7.2))
        for case,label in zip(cases,labels):
            if manifest_cases[case]["policy_status"]!="validated":
                continue
            for col,name,age_field in ((0,"consumption","mean_consumption"),(1,"next_assets","mean_next_assets")):
                r=rows(case,name+"_ecdf")
                x=[float(v["value"]) for v in r];y=[float(v["cdf"]) for v in r]
                axs[0,col].plot(x,y,color=COLORS[label],lw=2,label=label)
                d=manifest_cases[case]["policy_summary"]["aggregate_by_age"]
                axs[1,col].plot([int(v["age"]) for v in d],[float(v[age_field]) for v in d],color=COLORS[label],lw=2,label=label)
        style(axs[0,0],"Consumption distribution","Four-year flow","Cumulative share")
        style(axs[0,1],"Next financial assets distribution","Financial assets","Cumulative share")
        style(axs[1,0],"Mean consumption by age","Age-cell start","Four-year flow")
        style(axs[1,1],"Mean next assets by age","Age-cell start","Financial assets")
        for col,name in ((0,"consumption"),(1,"next_assets")):
            q=manifest_cases[cases[0]]["policy_summary"][name]
            axs[0,col].set_xlim(q["q01"],q["q99"])
        axs[0,0].legend(fontsize=8)
        fig.suptitle(title+" | Validated policies",fontsize=14,weight="bold")
        figs.append((savefig(pair+"_06_policies"),"Only the phi=0.8 aggregate policy reporter passes its native value-function gate. Phi=1 is unavailable because positive buyer mass occupies initialized-infeasible value entries; no mass is filtered. Consumption is a four-year flow and next financial assets are a stock, not national-account saving."))
    return figs


def compact_tables(cases):
    specs=[]
    overview=[];age=[];room=[];quant=[]
    selected_ages={18,30,42,62,82}
    for cid,d in cases.items():
        a=rows(cid,"age_summary")
        co=rows(cid,"age_children")
        ow=rows(cid,"owned_rooms")
        wr=rows(cid,"beginning_assets")
        owner_mass=sum(float(v["owner_mass"]) for v in ow)
        owner_mean=sum(float(v["rooms"])*float(v["owner_mass"]) for v in ow)/owner_mass if owner_mass else None
        at42=[v for v in co if int(v["age_start"])==42 and v["kind"]=="ever_born"]
        overview.append({"case":cid,"population":d["population"],"phi":d["label"],
                         "price_held":d["price"],"ownership_all_ages":owner_mass,
                         "childless_share_age42":next(float(v["share"]) for v in at42 if v["count"]=="0"),
                         "two_plus_share_age42":sum(float(v["share"]) for v in at42 if v["count"] in ("2","3+")),
                         "owner_mean_rooms":owner_mean,
                         "renter_mean_rooms":d["renter_quantiles"]["mean"] if d["renter_quantiles"] else None,
                         "policy_reporter":d["policy_status"]})
        for v in a:
            if int(v["age_start"]) not in selected_ages:continue
            ar=[x for x in co if x["kind"]=="ever_born" and x["age_start"]==v["age_start"]]
            age.append({"case":cid,"age_cell":f"{v['age_start']}-{v['age_end']}",
                        "ownership":v["ownership"],"mean_rooms":v["mean_rooms"],
                        "children_0":next(x["share"] for x in ar if x["count"]=="0"),
                        "children_1":next(x["share"] for x in ar if x["count"]=="1"),
                        "children_2":next(x["share"] for x in ar if x["count"]=="2"),
                        "children_3plus":next(x["share"] for x in ar if x["count"]=="3+"),
                        "children_at_home_mean":v["mean_children_at_home"]})
        for v in ow:
            room.append({"case":cid,"rooms":v["rooms"],"share_among_owners":v["share_among_owners"]})
        wa=rows(cid,"assets_by_age")
        nw=rows(cid,"beginning_gross_net_worth")
        pooled={"q10":weighted_quantile([float(v["asset_node"]) for v in wr],[float(v["mass"]) for v in wr],.1),
                "q50":weighted_quantile([float(v["asset_node"]) for v in wr],[float(v["mass"]) for v in wr],.5),
                "q90":weighted_quantile([float(v["asset_node"]) for v in wr],[float(v["mass"]) for v in wr],.9)}
        quant.append({"case":cid,"beginning_assets_q10":pooled["q10"],"beginning_assets_q50":pooled["q50"],
                      "beginning_assets_q90":pooled["q90"],
                      "gross_net_worth_q50":weighted_quantile([float(v["gross_net_worth"]) for v in nw],
                                                           [float(v["mass"]) for v in nw],.5),
                      "consumption_q50":d["policy_summary"].get("consumption",{}).get("q50"),
                      "next_assets_q50":d["policy_summary"].get("next_assets",{}).get("q50"),
                      "policy_reporter":d["policy_status"]})
    for stem,title,data in (("overview","Four-case outcome overview",overview),
                            ("selected_ages","Age-cell distribution summary: selected ages",age),
                            ("owner_rooms","Owned room-size shares",room),
                            ("asset_quantiles","Financial-asset and validated policy quantiles",quant)):
        p=HERE/"tables"/f"comparison_{stem}.csv"
        write_csv(p,list(data[0]),data)
        specs.append({"path":str(p),"title":title,"section":"comparison"})
    return specs


if __name__ == "__main__":
    HERE.mkdir(parents=True,exist_ok=True)
    manifest_cases={}
    for case_id,(path,pop,label,native) in CASES.items():
        print("extract",case_id,flush=True)
        manifest_cases[case_id]=extract(case_id,path,pop,label,native)
    figures=[]
    for pair,title,valid in (("reference","Heterogeneous productivity",True),("low","Permanently lowest productivity",True)):
        title_map={"children":"Children ever born and currently at home",
                   "ownership":"Ownership across age and family size",
                   "housing":"Housing services and room sizes",
                   "wealth":"Beginning financial wealth",
                   "net_worth":"Pre-transaction gross net worth",
                   "policies":"Validated consumption and next-asset policies"}
        for i,(path,caption) in enumerate(make_figures(pair,title,valid),1):
            stem=Path(path).stem.split("_",2)[-1]
            figures.append({"id":f"{pair}_{i:02d}","population":pair,"path":path,
                            "title":title_map.get(stem,stem.replace("_"," ").title()),"caption":caption})
    out={"title":"Household distributions under two fixed-price credit comparisons",
         "cases":manifest_cases,"figures":figures,"tables":compact_tables(manifest_cases),
         "notes":["All four cases hold the house price at 0.7760569760205563 and are fixed-price comparisons, not general equilibria or recalibrations.",
                  "Within each population, only the uniform financed share phi differs: 0.8 versus 1.0. Phi=1 retains the collateral floor and other purchase screens.",
                  "The low-productivity experiment also changes every household to the original lowest productivity node, retains the original marginal entrant wealth distribution, and recomputes the balanced pension under the unchanged payroll tax.",
                  "Both phi=1 aggregate reporters fail because positive buyer mass occupies initialized-infeasible value entries. The low-productivity case has a flagged buyer cell of mass 2.681e-15; the reference phi=1 case has total flagged buyer mass 5.258e-38. Consumption and next-asset distributions are unavailable for both phi=1 arms; no gate was altered or mass filtered.",
                  "The 17 inherited standard diagnostic plots per case are supplementary. Their housing-market residual panels evaluate a held price, not market clearing."],
         "units":{"assets":"stock in original mean annual gross-earnings units","consumption":"four-year flow in original mean annual gross-earnings units","rooms":"physical rooms","age":"four-year model cells starting at 18"}}
    (HERE/"manifest.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print("wrote",HERE/"manifest.json",len(figures),"figures",flush=True)
