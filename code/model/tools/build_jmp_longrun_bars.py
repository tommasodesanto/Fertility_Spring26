"""Long-run (stationary GE at the 2023 psi) effects of the rebated property tax and of mortgage credit (JMP slide).

Reads output/model/reconciliation_14p40_20261007/psi2023/results.json (no solve) and writes
output/model/jmp_slides_policy_functions_20261007/longrun_bars.{pdf,png}.
Run: output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/tools/build_jmp_longrun_bars.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "output/model/reconciliation_14p40_20261007/psi2023/results.json"
OUT = ROOT / "output/model/jmp_slides_policy_functions_20261007"

d = json.loads(SRC.read_text())
base = d["psi110_new_base"]
arms = [("Property tax 1.06% to 2%,\nrebated", d["psi110_new_ptax2"]), ("Financed share\n80% to 95%", d["psi110_new_phi095"])]
panels = [
    ("House price (% change)", lambda a: 100 * (a["price"] / base["price"] - 1)),
    ("Population (% change)", lambda a: 100 * (a["pop"] / base["pop"] - 1)),
    ("Ownership, ages 18-29 (pp)", lambda a: 100 * (a["own_18_29"] - base["own_18_29"])),
]
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 3, figsize=(11, 3.2), sharey=True)
colors = ["#2a7f62", "#08519c"]
for ax, (title, f) in zip(axes, panels):
    vals = [f(a) for _, a in arms]
    y = [1, 0]
    ax.barh(y, vals, color=colors, height=0.55)
    ax.axvline(0, color="0.3", lw=0.8)
    span = max(abs(v) for v in vals)
    for yi, v in zip(y, vals):
        ax.text(v + (0.04 * span if v >= 0 else -0.04 * span), yi, f"{v:+.1f}", va="center",
                ha="left" if v >= 0 else "right", fontsize=10)
    ax.set_xlim(-1.35 * span if min(vals) < 0 else -0.15 * span, 1.35 * span if max(vals) > 0 else 0.15 * span)
    ax.set_title(title, loc="left", fontsize=11)
    ax.grid(axis="x", alpha=0.3)
axes[0].set_yticks([1, 0])
axes[0].set_yticklabels([n for n, _ in arms])
fig.tight_layout()
OUT.mkdir(parents=True, exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"longrun_bars.{ext}", dpi=200, bbox_inches="tight")
print("wrote", OUT / "longrun_bars.pdf", [round(f(a), 2) for _, f in panels for _, a in arms])
