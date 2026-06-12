import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

alpha, beta, kappa, chi, q = 0.5, 0.4, 0.3, 0.6, 1.0
y, a, gamma, ell = 1.0, 0.20, 0.5, 0.2
w = y + a            # lifetime resources
k = a / 0.6161538027919038   # vartheta(1+gamma); set so a' = a at H=0.25

def solve_constrained(H):
    # budget (1+k)hatc + chi n = w - qH ; FOC beta/n = chi/hatc + alpha*kappa/(H-kappa n)
    lo, hi = 1e-9, H/kappa - 1e-9
    def f(n):
        hatc = (w - q*H - chi*n)/(1.0+k)
        if hatc <= 0: return -1e9
        return beta/n - chi/hatc - alpha*kappa/(H - kappa*n)
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if f(mid) > 0: lo = mid
        else: hi = mid
    n = 0.5*(lo+hi)
    hatc = (w - q*H - chi*n)/(1.0+k)
    return hatc, n

D = 1 + alpha + beta + k
hatc_u = w / D
n_u = beta * w / (D * (chi + kappa*q))
h_u = alpha * w / (D * q) + kappa * n_u

# exact equilibrium point
hatc_e, n_e = solve_constrained(0.25)
s_e = 0.25 - kappa*n_e
zeta = alpha*hatc_e/s_e - q
# slope at 0.25
num = alpha*kappa/s_e**2 - chi*q/((1+k)*hatc_e**2)
den = beta/n_e**2 + alpha*kappa**2/s_e**2 + chi**2/((1+k)*hatc_e**2)
slope = num/den
print(f"check: n={n_e:.4f} hatc={hatc_e:.4f} zeta={zeta:.4f} slope={slope:.4f} h_u={h_u:.4f} n_u={n_u:.4f}")

BLUE, RED, GRAY = "#2b5d8a", "#b03a3a", "#666666"

fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.6, 4.1))

# ---- left: marginal values ----
hs = np.linspace(0.16, 0.40, 200)
mv_y = []
for H in hs:
    hc, n = solve_constrained(H)
    mv_y.append(alpha*hc/(H - kappa*n))
mv_y = np.array(mv_y)
axL.plot(hs, mv_y, color=BLUE, lw=2.2, label="young buyer")

ho = np.linspace(0.06, 0.42, 200)
W_liq = a + (q-ell)*0.25     # 0.40
cO = W_liq - (q-ell)*ho
mv_o = gamma*cO/ho
axL.plot(ho, mv_o, color=RED, lw=2.2, ls="--", label="old incumbent")

axL.axhline(q, color=GRAY, lw=1.0, ls=":")
axL.text(0.395, q+0.03, "$q=1$", color=GRAY, ha="right", fontsize=10)

axL.plot([0.25], [q+zeta], "o", color=BLUE, ms=6)
axL.plot([1/6], [q-ell], "o", color=RED, ms=6)

axL.annotate("", xy=(0.25, q+zeta), xytext=(0.25, q),
             arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.4))
axL.text(0.256, 1.30, r"$\zeta^{O,F}=0.68$", color=BLUE, fontsize=10)
axL.annotate("", xy=(1/6, q-ell), xytext=(1/6, q),
             arrowprops=dict(arrowstyle="->", color=RED, lw=1.4))
axL.text(0.105, 0.86, r"$\ell=0.2$", color=RED, fontsize=10)

axL.text(0.205, 1.86, "move one unit:\nsurplus $0.88$", fontsize=9.5, ha="center",
         bbox=dict(boxstyle="round,pad=0.3", fc="#f5f0e6", ec="#b8a268", lw=0.8))

axL.set_xlabel("floorspace")
axL.set_ylabel("marginal value of space (goods)")
axL.set_xlim(0.05, 0.41); axL.set_ylim(0.55, 2.15)
axL.legend(frameon=False, fontsize=9.5, loc="upper right")
axL.set_title("Two-sided misallocation", fontsize=11)

# ---- right: n(H) ----
Hs = np.linspace(0.10, 0.40, 300)
ns = []
for H in Hs:
    if H < h_u:
        _, n = solve_constrained(H)
    else:
        n = n_u
    ns.append(n)
axR.plot(Hs, ns, color=BLUE, lw=2.2)
axR.plot([0.25], [n_e], "o", color=BLUE, ms=6)
axR.axvline(h_u, color=GRAY, lw=1.0, ls=":")
axR.text(h_u+0.004, 0.165, r"cap stops binding: $h^u=0.34$", rotation=90,
         color=GRAY, fontsize=9, va="bottom")

dx = 0.045
axR.plot([0.25-dx, 0.25+dx], [n_e-slope*dx, n_e+slope*dx], color=RED, lw=1.6, ls="-")
axR.text(0.252, n_e-0.012, r"$\dd n/\dd H=0.33$".replace(r"\dd","d"), color=RED, fontsize=10)

axR.set_xlabel("effective family-housing cap $H$")
axR.set_ylabel("completed fertility $n(H)$")
axR.set_xlim(0.10, 0.40); axR.set_ylim(0.14, 0.26)
axR.set_title("Fertility against the cap", fontsize=11)

for ax in (axL, axR):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

fig.tight_layout()
out = "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/figures/example_misallocation.pdf"
fig.savefig(out, bbox_inches="tight")
print("saved", out)
