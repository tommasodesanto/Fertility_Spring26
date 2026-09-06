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

    # Both panels compare the SAME initial, impact and final equilibria.
    # Arrows connect those states schematically; they are not date-by-date
    # trajectories. The full analytical response, including overshooting, is
    # retained in the receipt and the existing supporting appendix. No static
    # P(population) or n(P) schedule is assumed. Distances are first-order
    # derivatives of levels, with symbolic level endpoints on the axes.
    fig,axes=plt.subplots(1,2,figsize=(12.5,3.8),layout="constrained")
    panels=[(population,price,np.array([n1,p1])),
            (price,fertility,np.array([p1,0.]))]
    for j,(ax,(xx,yy,endpoint)) in enumerate(zip(axes,panels)):
        impact=np.array([xx[0],yy[0]])
        ax.annotate("",xy=impact,xytext=(0,0),
                    arrowprops=dict(arrowstyle="->",color="#999999",lw=1.7,
                                    shrinkA=7,shrinkB=7,mutation_scale=14),zorder=1)
        ax.annotate("",xy=endpoint,xytext=impact,
                    arrowprops=dict(arrowstyle="->",color=RED,lw=2.,linestyle=(0,(5,3)),
                                    shrinkA=7,shrinkB=7,mutation_scale=14),zorder=1)
        ax.scatter([0],[0],color=BLUE,s=70,zorder=5)
        ax.scatter([impact[0],endpoint[0]],[impact[1],endpoint[1]],color=RED,s=70,zorder=5)
        ax.annotate(r"$A$",(0,0),xytext=(8,10),textcoords="offset points",color=BLUE,fontsize=18)
        ax.annotate(r"$I$",impact,xytext=(9,2 if j==0 else 5),textcoords="offset points",color=RED,fontsize=18)
        ax.annotate(r"$A'$",endpoint,xytext=(9,8 if j==0 else -23),textcoords="offset points",color=RED,fontsize=18)

    ax=axes[0]
    ax.set(title="(a) Housing market",xlabel=r"Adult households, $Y+O$",ylabel=r"House price, $P$")
    ax.set_xticks([0,n1],[r"$Y_0+O_0$",r"$Y_1^*+O_1^*$"])
    ax.set_yticks([0,p1],[r"$P_0^*$",r"$P_1^*$"])
    ax.set_xlim(-.20,2.90);ax.set_ylim(-.22,5.05)
    ax.text(.11,1.9,"Impact",color=GRAY,fontsize=13)
    ax.text(.88,3.50,"Population adjustment",color=RED,fontsize=13)
    ax=axes[1]
    ax.axhline(0,color=GRAY,lw=1,ls=":")
    ax.set(title="(b) House prices and fertility",xlabel=r"House price, $P$",ylabel=r"Mean fertility, $\bar n$")
    ax.set_xticks([0,p1],[r"$P_0^*$",r"$P_1^*$"])
    ax.set_yticks([0],[r"$1/\nu$"])
    ax.set_xlim(-.25,4.92);ax.set_ylim(-.020,1.38*fertility[0])
    for ax in axes:
        ax.title.set_fontsize(18)
        ax.xaxis.label.set_size(17)
        ax.yaxis.label.set_size(17)
        ax.tick_params(length=4,pad=6,labelsize=17)
    save(fig,"theory_slides_transition")
    return dict(taste_scale=p["sigma"],owner_share=float(hh["pi"]),
                first_fertility_derivative=float(fertility[0]),
                stationary_population_derivative=float(n1),
                stationary_price_derivative=float(p1),
                map_residual=float(map_error),initial_boundary_residual=float(boundary_error),
                states=states.tolist(),fertility=fertility.tolist(),population=population.tolist(),
                interpretation="Initial, impact and final equilibria from the local analytical credit response. Arrows connect those states schematically; intermediate overshooting is not drawn. The complete analytical path remains in the receipt and supporting appendix. No static one-dimensional schedule or monotonicity claim.")


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
