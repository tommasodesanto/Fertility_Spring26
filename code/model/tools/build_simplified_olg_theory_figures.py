"""Draw the local analytical transition; no nonlinear model is solved.

The cubic and jump formulas are proved in the consolidated note's transition
appendix. Roots are evaluated only to display that proved local solution.
Stationary schedules and all paths use the same first-order approximation.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/simplified_olg_amendments"


def build():
    q = beta = 0.75
    alpha = gamma = theta = kappa = nu = 1.0
    omega, W, V, chi, pi = 4.0, 6.0, 24.48, 0.6, 0.25
    E, K, discount_gap = 1 + alpha + theta, 1 + gamma + omega, 1 - q
    X = W / E
    p = (nu * theta * X - chi) / kappa
    P = p / discount_gap
    eps = kappa * p / (chi + kappa * p)
    A, B = alpha + theta * eps, alpha + theta * eps**2
    hy, ho = A * X / p, gamma * V / (K * p)
    Hpair, zeta = hy + ho, ho / hy
    R = pi * beta * K * W / (q * E * V)
    f = q / discount_gap * ((1 - pi + R) * eps - R * A / E)
    a = q / discount_gap * ((1 - pi + R) * B / A - R * A / E + zeta)
    b, c = B / A + zeta, pi * gamma / (K * discount_gap)
    coefficients = [a, f - 2*a - b + c,
                    f*(zeta-1) - eps + a + b - 2*c,
                    c - zeta*(f+eps)]
    roots = np.roots(coefficients)
    stable = roots[np.abs(roots) < 1]
    unstable = roots[np.abs(roots) > 1]
    assert len(stable) == 2 and len(unstable) == 1
    ru = float(np.real(unstable[0]))
    assert abs(np.imag(unstable[0])) < 1e-10
    assert alpha/(1+alpha) < eps < q
    assert 1 < zeta < 2*E*eps/A - 1
    assert pi*gamma/K < min(discount_gap**2, q*zeta)
    assert q*E*(V/W) > beta*K
    assert omega*discount_gap > q*gamma

    def response(g, forcing, periods=24):
        Bu = 1 + zeta / ru
        Du = a*(ru-1) + f*Bu
        zs = g/eps
        ys = (b*zs - forcing)/(1+zeta)
        z0 = (Bu*g/(ru-1) + forcing)/Du
        dz0 = ((b-c)*z0 - forcing)/a
        y1 = g - eps*z0 + f*dz0
        matrix = np.array([[1,1], stable], dtype=complex)
        cy = np.linalg.solve(matrix, [-ys, y1-ys])
        cz = np.linalg.solve(matrix, [z0-zs, z0+dz0-zs])
        t = np.arange(periods + 2)
        y = np.real(ys + cy[0]*stable[0]**t + cy[1]*stable[1]**t)
        z = np.real(zs + cz[0]*stable[0]**t + cz[1]*stable[1]**t)
        # These are checks on displayed analytical formulas, not a root proof.
        assert abs(y[0]) < 1e-11
        assert np.max(abs(np.diff(y) + eps*z[:-1] - f*np.diff(z) - g)) < 1e-10
        assert np.max(abs(y[1:-1] + zeta*y[:-2] - b*z[1:-1]
                          + a*np.diff(z)[1:] + c*np.diff(z)[:-1] + forcing)) < 1e-10
        return y, z, ys, zs

    J = A*(1+zeta)/E
    gt = q/discount_gap*(J/2-eps)
    Ct = q/discount_gap*((1+zeta)/2*(A/E+gamma/K)-(B/A+zeta))
    gy, Cy = 1/theta-1/E, eps/A-1/E
    yt, zt, yts, zts = response(gt, Ct)
    yy, zy, yys, zys = response(gy, Cy)
    dtheta, dtax, policy_date = -0.03, 0.012, 3
    ybase, zbase = dtheta*yy, dtheta*zy
    ypol, zpol = ybase.copy(), zbase.copy()
    ypol[policy_date:] += dtax*yt[:len(ypol)-policy_date]
    zpol[policy_date:] += dtax*zt[:len(zpol)-policy_date]

    def observed(y, z):
        adult = 1 + (y[:-1] + np.r_[0,y[:-2]])/2
        price = P*(1+z[:-1])
        fertility = (1+np.diff(y))/nu
        return adult, price, fertility
    nb, pb, fb = observed(ybase,zbase)
    npol, pp, fp = observed(ypol,zpol)
    ss0 = (1+dtheta*yys, P*(1+dtheta*zys), 1/nu)
    ss1 = (ss0[0]+dtax*yts, ss0[1]+P*dtax*zts, 1/nu)

    # Tangent stationary housing-clearing schedules, with equal age masses.
    a_space, b_space, c_space = alpha*X/p, gamma*V/(K*p), kappa/nu
    dh = a_space+b_space+eps*c_space
    rebate_derivative = q*P*Hpair/2
    home_theta = W*(kappa*p-alpha*chi)/(E**2*p*(chi+kappa*p))
    home_transfer = (alpha/p+kappa*theta/(chi+kappa*p))/E+gamma/(K*p)
    price_population = P*Hpair/dh
    price_taste = P*home_theta/dh
    price_tax = P*(home_transfer*rebate_derivative-dh*q/discount_gap)/dh
    fert_price = -kappa*discount_gap/(nu*(chi+kappa*p))

    plt.rcParams.update({"font.family":"DejaVu Serif", "mathtext.fontset":"dejavuserif",
                         "font.size":10, "axes.spines.top":False,
                         "axes.spines.right":False, "axes.labelsize":11})
    fig, axes = plt.subplots(1,2,figsize=(11.5,3.9))
    gray, blue, red = "#888888", "#2165a6", "#b63b3d"
    nlo, nhi = min(nb.min(),npol.min(),ss0[0],ss1[0])-.006, 1.006
    plo = min(pb.min(),pp.min(),ss0[1],ss1[1])-.06
    phi = max(P,pb.max(),pp.max())+.08
    gridn = np.linspace(nlo,nhi,180)
    gridp = np.linspace(plo,phi,180)
    schedules = [(0,0,gray,"Initial"), (dtheta,0,blue,"After the fertility decline"),
                 (dtheta,dtax,red,"With the tax reform")]
    for shock,tax,color,label in schedules:
        axes[0].plot(gridn,P+price_population*(gridn-1)+price_taste*shock+price_tax*tax,
                     color=color,lw=1.3,alpha=.7,label=label)
        axes[1].plot(gridp,1/nu+fert_price*(gridp-P)+gy*shock+gt*tax,
                     color=color,lw=1.3,alpha=.7)
    axes[0].plot(nb,pb,color=blue,lw=1.8,ls="--")
    axes[0].plot(npol[policy_date:],pp[policy_date:],color=red,lw=1.8,ls="--")
    axes[1].plot(pb,fb,color=blue,lw=1.8,ls="--")
    axes[1].plot(pp[policy_date:],fp[policy_date:],color=red,lw=1.8,ls="--")
    axes[1].axhline(1/nu,color="#aaaaaa",lw=.8,ls=":")

    def point(ax,x,y,label,color,offset):
        ax.scatter([x],[y],s=21,color=color,zorder=6)
        ax.annotate(label,(x,y),xytext=offset,textcoords="offset points",color=color,
                    fontsize=10,ha="center",va="center",zorder=7)
    def arrow(ax,start,end,color):
        ax.annotate("",end,start,arrowprops={"arrowstyle":"->","color":color,"lw":1.4},zorder=5)
    point(axes[0],1,P,r"$S_-$",gray,(9,9))
    point(axes[1],P,1,r"$S_-$",gray,(10,10))
    for ss,label,color,offset0,offset1 in [
        (ss0,r"$S_0$",blue,(-11,-11),(12,12)),
        (ss1,r"$S_1$",red,(10,-10),(-10,-12))]:
        point(axes[0],ss[0],ss[1],label,color,offset0)
        point(axes[1],ss[1],ss[2],label,color,offset1)
    i=policy_date
    point(axes[0],nb[i],pb[i],"",blue,(13,3))
    point(axes[0],npol[i],pp[i],r"$I$",red,(12,-6))
    point(axes[1],pb[i],fb[i],"",blue,(10,-10))
    point(axes[1],pp[i],fp[i],r"$I$",red,(-12,10))
    arrow(axes[0],(1,P),(nb[0],pb[0]),gray)
    arrow(axes[1],(P,1),(pb[0],fb[0]),gray)
    arrow(axes[0],(nb[i],pb[i]),(npol[i],pp[i]),red)
    arrow(axes[1],(pb[i],fb[i]),(pp[i],fp[i]),red)
    for ax,xb,yb,xp,yp in [(axes[0],nb,pb,npol,pp),(axes[1],pb,fb,pp,fp)]:
        arrow(ax,(xb[1],yb[1]),(xb[2],yb[2]),blue)
        arrow(ax,(xp[i+1],yp[i+1]),(xp[i+2],yp[i+2]),red)
    axes[0].set(xlim=(nlo,nhi),ylim=(plo,phi),xlabel=r"Adult households $N_{hh}$",ylabel=r"House price $P$")
    fmin, fmax = min(fb.min(),fp.min())-.005, max(1,fb.max(),fp.max())+.006
    axes[1].set(xlim=(plo,phi),ylim=(fmin,fmax),xlabel=r"House price $P$",ylabel=r"Mean fertility $\bar n$")
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
        ax.spines['left'].set_color('#777777');ax.spines['bottom'].set_color('#777777')
    axes[1].text(phi,1,r" $1/\nu$",va="center",fontsize=10)
    fig.legend(*axes[0].get_legend_handles_labels(),loc="upper center",ncol=3,frameon=False,bbox_to_anchor=(.5,1.01))
    fig.text(.5,.015,"Local equilibrium illustration. Solid: stationary market schedules. Dashed: first-order transition.",ha="center",fontsize=8,color='#555555')
    fig.subplots_adjust(left=.06,right=.955,bottom=.16,top=.87,wspace=.27)
    OUT.mkdir(parents=True,exist_ok=True)
    for suffix in ['pdf','png']:
        fig.savefig(OUT/f'consolidated_transition_figure.{suffix}',dpi=180)
    plt.close(fig)
    receipt={"purpose":"first-order analytical illustration, not calibration or a nonlinear transition solve",
             "endowments":"Compact positive heterogeneous w with mean W and v_i=(V/W)w_i; finite housing caps chosen slack. Only these distributional moments enter this first-order illustration.",
             "finance":"phi=q; all young finance constraints strictly binding at the reference",
             "parameters":dict(q=q,beta=beta,alpha=alpha,gamma=gamma,theta=theta,kappa=kappa,nu=nu,omega=omega,W=W,V=V,chi=chi,pi=pi),
             "taste_change":dtheta,"tax_change":dtax,"policy_date":policy_date,
             "roots":[[float(z.real),float(z.imag)] for z in roots],
             "stationary_baseline":list(ss0),"stationary_policy":list(ss1),
             "tax_impact_fertility":float(yt[1]),"taste_impact_fertility":float(yy[1])}
    (OUT/'consolidated_transition_figure.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('Saved both-panel analytical transition figure and verification receipt.')


if __name__ == '__main__':
    build()
