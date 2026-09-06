#!/usr/bin/env python3
"""Build the seven-slide theory extract and its two economic diagrams.

Run from the repository root with Python/NumPy/SciPy/Matplotlib and TeX installed.
The main deck is the only slide source. Its quantitative sections are not edited.
The diagrams reuse the earlier marginal-value and two-panel transition layouts,
with the current model's compensation argument and analytical local transition.
This is a presentation build, not a calibration or a finite-horizon simulation.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import schur

from verify_simplified_olg_mixed_transition import anchor, linearization
from verify_simplified_olg_local_transition import (
    young_choices, old_choices, complex_jacobian,
)

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/simplified_olg_amendments"
TMP = ROOT / "tmp/pdfs/september_14_theory_compact"
PDF = ROOT / "output/pdf"
SOURCE = ROOT / "latex/september_14_presentation.tex"
BLUE, RED, GRAY = "#26527a", "#b0493d", "#737373"


def save(fig, name):
    for extension in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{extension}", dpi=190, bbox_inches="tight")
    plt.close(fig)


def style():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "cm", "font.size": 16,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 1., "pdf.fonttype": 42,
        "legend.frameon": False,
    })


def allocation_figure(p, hh):
    # The same verified mixed-tenure stationary economy as the transition.
    x, h, n, _, c2, h2, _ = hh["owner"]["z"]
    alpha, gamma, kappa = p["alpha"], p["gamma"], p["kappa"]
    s = h-kappa*n
    u = hh["u0"]

    def compensation(epsilon):
        return c2*((h2/(h2-epsilon))**gamma-1)

    def marginal_values(epsilon):
        D = compensation(epsilon)
        return alpha*(x-D)/(s+epsilon), gamma*(c2+D)/(h2-epsilon)

    epsilon = np.linspace(-.32, .28, 400)
    young, old = marginal_values(epsilon)
    young0, old0 = marginal_values(0.)
    move = .025
    young1, old1 = marginal_values(move)
    gain = np.log1p(-compensation(move)/x)+alpha*np.log1p(move/s)
    old_gain = np.log1p(compensation(move)/c2)+gamma*np.log1p(-move/h2)
    assert abs(old0-u)<1e-12 and young0>old0
    assert gain>0 and abs(old_gain)<1e-12
    assert x-compensation(move)>0 and s+move>0 and h2-move>0

    fig, ax = plt.subplots(figsize=(6.6, 4.3), layout="constrained")
    ax.plot(h+epsilon, young, color=BLUE, lw=2.7, label="Young buyer")
    ax.plot(h2-epsilon, old, color=RED, lw=2.7, ls=(0,(6,3)), label="Old owner")
    ax.axhline(u, color=GRAY, lw=1., ls=":")
    ax.text(1.29, u+.015, r"$q r_t$", color=GRAY, ha="right")
    for size, value, color in ((h, young0, BLUE), (h2, old0, RED)):
        ax.plot(size, value, "o", color=color, ms=7)
        ax.plot([size,size],[.30,value], color=color, lw=.9, ls=":")

    # Arrows follow one exactly compensated move of the same housing quantity.
    ax.annotate("", xy=(h+move,young1), xytext=(h,young0),
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.7,
                                mutation_scale=14, shrinkA=5, shrinkB=0))
    ax.annotate("", xy=(h2-move,old1), xytext=(h2,old0),
                arrowprops=dict(arrowstyle="->", color=RED, lw=1.7,
                                mutation_scale=14, shrinkA=5, shrinkB=0))
    middle = .73
    ax.plot([h2,middle],[old0,old0], color=GRAY, lw=.9, ls=":")
    ax.plot([middle,h],[young0,young0], color=GRAY, lw=.9, ls=":")
    ax.annotate("", xy=(middle,old0), xytext=(middle,young0),
                arrowprops=dict(arrowstyle="<->", color=GRAY, lw=1.2))
    ax.text(middle-.03,(old0+young0)/2,"Value\ngap",ha="right",va="center",
            color=GRAY,fontsize=13)
    ax.set(xlim=(.25,1.32),ylim=(.30,.95),xlabel="Housing",ylabel="Marginal value in goods")
    ax.set_xticks([h2,h],[r"$h_j^2$",r"$h_i$"])
    ax.set_yticks([])
    ax.legend(loc="upper right",fontsize=14)
    save(fig,"theory_slides_misallocation")
    return dict(young_initial=float(young0),old_initial=float(old0),rent_pv=float(u),
                housing_move=move,old_compensation=float(compensation(move)),
                young_utility_gain=float(gain),old_utility_change=float(old_gain),
                curve_definition="Exact old-utility compensation, individual fertility and future real allocations fixed. Curves are plotted against each household's own housing, not against the transfer.")


def transition_figure(p, Z, hh):
    result=linearization(p,Z,hh)
    J,stat=result["J"],result["derivative"]
    tri,E,k=schur(J,output="real",sort=lambda re,im:np.hypot(re,im)<1)
    assert k==4
    Es,stable=E[:,:k],tri[:k,:k]
    a=np.linalg.solve(result["B"]@Es,-result["B"]@stat)
    states=np.array([stat+Es@np.linalg.matrix_power(stable,t)@a for t in range(41)])
    fertility=np.diff(states[:,2])/(p["nu"]*Z[2])
    population=states[:-1,2]+states[:-1,3]
    price=states[:-1,0]
    p1,n1=stat[0],stat[2]+stat[3]
    map_error=np.max(abs(states[1:]-states[:-1]@J.T-result["Dphi"]))
    boundary_error=np.max(abs(result["B"]@states[0]))
    assert map_error<1e-11 and boundary_error<1e-11
    assert fertility[0]>0 and n1>0 and abs(population[0])<1e-12
    assert np.max(abs(states[-1]-stat))<1e-10

    curves = equilibrium_curves(p, Z, hh, states, result)
    # These are first-order conditional housing-market schedules, not a
    # one-state law of motion. Initial and impact schedules keep initial old
    # assets and the corresponding equilibrium future prices/rebate fixed.
    # The final schedule uses constant prices and repeated lifetime choices. Each marked
    # point satisfies the full equilibrium. A monotone vertical transformation
    # makes the small fertility overshoot legible on explicitly schematic axes;
    # it preserves curve ordering, slopes' signs, intersections, and the points.
    fy = lambda n: np.arcsinh(np.asarray(n)/.055)
    fig,axes=plt.subplots(1,2,figsize=(12.5,4.45))
    fig.subplots_adjust(left=.075,right=.985,bottom=.19,top=.77,wspace=.32)
    Ngrid=np.linspace(-.62,3.05,250)
    Pgrid=np.linspace(-.90,4.98,350)
    schedules=[("Initial",BLUE,"-",curves["impact"],0.),
               ("Impact",GRAY,":",curves["impact"],1.),
               ("Long run",RED,(0,(6,3)),curves["stationary"],1.)]
    handles=[]
    for label,color,ls,curve,shock in schedules:
        line,=axes[0].plot(Ngrid,(Ngrid-curve["N_shift"]*shock)/curve["N_P"],
                          color=color,ls=ls,lw=2.6,label=label)
        axes[1].plot(Pgrid,fy(curve["n_P"]*Pgrid+curve["n_shift"]*shock),
                     color=color,ls=ls,lw=2.6)
        handles.append(line)
    fig.legend(handles=handles,labels=[r"Initial: $\phi_0$",r"Impact: $\phi_1$",r"Long run: $\phi_1$"],
               loc="upper center",bbox_to_anchor=(.52,1.01),ncol=3,fontsize=15,
               handlelength=2.6,columnspacing=2.4)
    panels=[(np.array([population[0],price[0]]),np.array([n1,p1])),
            (np.array([price[0],fy(fertility[0])]),np.array([p1,0.]))]
    for j,(ax,(impact,endpoint)) in enumerate(zip(axes,panels)):
        ax.annotate("",xy=impact,xytext=(0,0),
                    arrowprops=dict(arrowstyle="->",color="#999999",lw=1.7,
                                    shrinkA=7,shrinkB=7,mutation_scale=14),zorder=1)
        ax.annotate("",xy=endpoint,xytext=impact,
                    arrowprops=dict(arrowstyle="->",color="#555555",lw=1.7,
                                    shrinkA=7,shrinkB=7,mutation_scale=14),zorder=4)
        ax.scatter([0],[0],color=BLUE,s=70,zorder=5)
        ax.scatter([impact[0],endpoint[0]],[impact[1],endpoint[1]],color=RED,s=70,zorder=5)
        ax.annotate(r"$A$",(0,0),xytext=(8,10),textcoords="offset points",color=BLUE,fontsize=18)
        ax.annotate(r"$I$",impact,xytext=(-16,9) if j==0 else (-18,6),textcoords="offset points",color=RED,fontsize=18)
        ax.annotate(r"$A'$",endpoint,xytext=(9,-20),textcoords="offset points",color=RED,fontsize=18)

    ax=axes[0]
    ax.set(title="(a) Housing-market clearing",xlabel=r"Adult households, $Y+O$",ylabel=r"House price, $P$")
    ax.set_xticks([0,n1],[r"$Y_0+O_0$",r"$Y_1^*+O_1^*$"])
    ax.set_yticks([0,p1],[r"$P_0^*$",r"$P_1^*$"])
    ax.set_xlim(-.62,3.05);ax.set_ylim(-.5,5.5)
    for xv,yv in [(0,0),(n1,p1)]:
        ax.plot([xv,xv],[-.5,yv],color=GRAY,ls=":",lw=.85)
        ax.plot([-.62,xv],[yv,yv],color=GRAY,ls=":",lw=.85)
    ax.annotate(r"$\widetilde P$",(0,price[0]),xytext=(9,-21),textcoords="offset points",fontsize=15)
    ax=axes[1]
    ax.axhline(0,color=GRAY,lw=1,ls=":")
    ax.set(title="(b) House prices and fertility",xlabel=r"House price, $P$",ylabel=r"Mean fertility, $\bar n$")
    ax.set_xticks([0,p1],[r"$P_0^*$",r"$P_1^*$"])
    ax.set_yticks([0,fy(fertility[0])],[r"$1/\nu$",r"$\widetilde n$"])
    ax.set_xlim(-.90,4.98);ax.set_ylim(-1.1,4.25)
    ax.plot([-.90,price[0]],[fy(fertility[0])]*2,color=GRAY,lw=.8,ls=":")
    ax.plot([p1,p1],[-1.1,0],color=GRAY,lw=.8,ls=":")
    for ax in axes:
        ax.title.set_fontsize(16)
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        ax.tick_params(length=4,pad=6,labelsize=16)
    save(fig,"theory_slides_transition")
    return dict(taste_scale=p["sigma"],owner_share=float(hh["pi"]),
                first_fertility_derivative=float(fertility[0]),
                stationary_population_derivative=float(n1),
                stationary_price_derivative=float(p1),
                map_residual=float(map_error),initial_boundary_residual=float(boundary_error),
                states=states.tolist(),fertility=fertility.tolist(),population=population.tolist(),
                curves=curves,
                interpretation="Conditional local equilibrium curves with initial, impact and final equilibria. Initial/impact curves hold inherited old states and their respective equilibrium future prices/rebate fixed; the long-run curve uses stationary lifetime choices. Arrows compare the states, not every transition date. Symbolic axes are schematic; fertility uses a monotone arcsinh transformation. No one-state dynamics or monotonic adjustment claim.")


def equilibrium_curves(p, Z, hh, states, result):
    """Derive both panels from the original household and market equations.

    F is excess housing demand, with equal young/old masses N/2. Eliminating
    N from F=0 yields fertility along the SAME conditional market schedule.
    The old-age asset distribution is generated by repeated lifetime choices
    for stationary curves, and inherited from the pre-reform economy at impact.
    Only the marked A/I/A' points also satisfy all demographic and expectation
    conditions. Off-point schedules deliberately relax those conditions.
    """
    q,tau,H=p["q"],p["tau"],p["Hbar"]
    T=q*tau*Z[0]*H/(Z[2]+Z[3])

    def stationary(v):
        P,N,phi=v
        rebate=q*tau*P*H/N
        choices=young_choices([P]*3,[rebate]*2,dict(p,phi=phi))
        old_h=sum(prob*choices[tenure]["z"][5] for tenure,prob in
                  (("owner",choices["pi"]),("renter",1-choices["pi"])))
        return np.array([N/2*(choices["housing"]+old_h)-H,choices["fertility"]])

    def impact(v):
        P,N,phi,Pn,Pnn,Tn=v
        rebate=q*tau*P*H/N
        pp=dict(p,phi=phi)
        choices=young_choices([P,Pn,Pnn],[rebate,Tn],pp)
        old_h=sum(prob*old_choices(hh[tenure]["assets"],hh[tenure]["z"][1],
                                  P,Pn,rebate,pp,tenure=="owner")[1]
                  for tenure,prob in (("owner",hh["pi"]),("renter",1-hh["pi"])))
        return np.array([N/2*(choices["housing"]+old_h)-H,choices["fertility"]])

    Js=complex_jacobian(stationary,[Z[0],Z[2]+Z[3],p["phi"]])
    Ji=complex_jacobian(impact,[Z[0],Z[2]+Z[3],p["phi"],Z[0],Z[0],T])
    Tn=q*tau*H*(states[1,0]/2-Z[0]*(states[1,2]+states[1,3])/4)
    direction=np.array([0.,0.,1.,states[1,0],states[2,0],Tn])

    def eliminate(J,shift):
        N_P=-J[0,0]/J[0,1]
        N_shift=-shift[0]/J[0,1]
        return dict(N_P=float(N_P),N_shift=float(N_shift),
                    n_P=float(J[1,0]+J[1,1]*N_P),
                    n_shift=float(shift[1]+J[1,1]*N_shift))

    cs=eliminate(Js,Js[:,2]); ci=eliminate(Ji,Ji@direction)
    Pss=-cs["n_shift"]/cs["n_P"]
    Nss=cs["N_P"]*Pss+cs["N_shift"]
    Pi=-ci["N_shift"]/ci["N_P"]
    ni=ci["n_P"]*Pi+ci["n_shift"]
    errors=[abs(Pss-result["derivative"][0]),
            abs(Nss-sum(result["derivative"][2:4])),
            abs(Pi-states[0,0]),abs(ni-(states[1,2]-states[0,2])/p["nu"])]
    assert max(errors)<1e-11 and cs["N_P"]>0 and ci["N_P"]>0
    assert cs["n_P"]<0 and ci["n_P"]<0
    # Independent central differences check the derivatives used to draw curves.
    finite_errors=[]
    for fun,point,J in ((stationary,np.array([1.,2.,p["phi"]]),Js),
                        (impact,np.array([1.,2.,p["phi"],1.,1.,T]),Ji)):
        step=1e-5
        for j in range(len(point)):
            dx=np.zeros_like(point);dx[j]=step
            finite_errors.append(float(np.max(abs((fun(point+dx)-fun(point-dx))/(2*step)-J[:,j]))))
    assert max(finite_errors)<1e-7
    return dict(stationary=cs,impact=ci,
                equilibrium_point_max_error=float(max(errors)),
                finite_difference_max_error=max(finite_errors),
                stationary_jacobian=Js.tolist(),impact_jacobian=Ji.tolist(),
                impact_exogenous_direction=direction.tolist(),
                scope="Local tangents at the certified mixed-tenure economy, not global schedules. Equal age shares on both plotted schedules; inherited old assets and transition expectations at impact, repeated lifetime choices and constant prices for stationary curves. Only marked points impose every equilibrium condition.")


def compile_decks():
    source=SOURCE.read_text()
    preamble=source.split(r"\begin{document}",1)[0]
    block=source.split(r"\section{Simple Theory}",1)[1].split(r"\section{Quantitative Model}",1)[0]
    block=block.split(r"\end{frame}",1)[1]  # Omit the section divider.
    assert block.count(r"\begin{frame}")==7
    extract=TMP/"simplified_olg_theory_slides.tex"
    extract.write_text(preamble+r"\hypersetup{pdftitle={Housing Allocation across Generations}}"+"\n"+
                       r"\begin{document}"+"\n"+block+r"\end{document}"+"\n")
    engine=shutil.which("pdflatex") or "/Library/TeX/texbin/pdflatex"
    env=dict(os.environ)
    env["PATH"]=str(Path(engine).parent)+os.pathsep+env.get("PATH","")
    for src in (SOURCE,extract):
        for run in (1,2):
            with (TMP/f"{src.stem}.pass{run}.txt").open("w") as log:
                subprocess.run([engine,"-interaction=nonstopmode","-halt-on-error","-file-line-error",
                                f"-output-directory={TMP}",str(src)],cwd=ROOT/"latex",env=env,
                               stdout=log,stderr=subprocess.STDOUT,check=True)
        shutil.copy2(TMP/f"{src.stem}.pdf",PDF/f"{src.stem}.pdf")
    shutil.copy2(PDF/"september_14_presentation.pdf",SOURCE.with_suffix(".pdf"))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figures-only",action="store_true")
    args=parser.parse_args()
    for folder in (OUT,TMP,PDF):
        folder.mkdir(parents=True,exist_ok=True)
    style()
    p,Z,hh=anchor(4)
    receipt={"allocation":allocation_figure(p,hh),"transition":transition_figure(p,Z,hh)}
    receipt["parameters"]={k:float(v) for k,v in p.items()}
    receipt["sources"]={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in (Path(__file__),SOURCE,
                                     Path(__file__).with_name("verify_simplified_olg_mixed_transition.py"),
                                     Path(__file__).with_name("verify_simplified_olg_local_transition.py"))}
    if not args.figures_only:
        compile_decks()
    (OUT/"theory_slides_figure_checks.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print("Built the allocation and transition diagrams"+("." if args.figures_only else " and both slide PDFs."))
    print(f"Allocation: young gain {receipt['allocation']['young_utility_gain']:.6g}; old change {receipt['allocation']['old_utility_change']:.2g}.")
    print(f"Transition: linear equilibrium residual {receipt['transition']['map_residual']:.2g}.")


if __name__=="__main__":
    main()
