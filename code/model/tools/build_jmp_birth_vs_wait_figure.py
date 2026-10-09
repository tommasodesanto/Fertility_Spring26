"""Birth versus wait from the same state (JMP slides): rooms, consumption and ownership, childless ages 22-33.

Reads the saved branch table (no model solve):
  output/model/rental_menu_precaution_20261007/tables_branches.md  (base 14.402, fixed price, rebate held)
Writes output/model/jmp_slides_policy_functions_20261007/birth_vs_wait.{pdf,png}.
Run: output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/tools/build_jmp_birth_vs_wait_figure.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "output/model/rental_menu_precaution_20261007/tables_branches.md"
OUT = ROOT / "output/model/jmp_slides_policy_functions_20261007"

# (table group, income z) -> slide label
ROWS = [
    ("renters, b = 0", "0.47", "Renter, no savings,\nlow earnings"),
    ("renters, b = 0", "0.78", "Renter, no savings,\nmiddle earnings"),
    ("owners at the collateral floor", "0.78", "Owner at borrowing\nlimit, middle earnings"),
    ("renters, b > 0", "1.29", "Renter with savings,\nhigh earnings"),
    ("owners above the floor", "1.29", "Owner above limit,\nhigh earnings"),
]


def read_rows():
    lines = SRC.read_text().split("## Who would use")[0].splitlines()
    table = {}
    for ln in lines:
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        if len(cells) == 9 and "/" in cells[4]:
            table[(cells[0], cells[1])] = cells
    out = []
    for grp, z, label in ROWS:
        c = table[(grp, z)]  # KeyError if the source table changes
        own = [float(x) for x in c[4].split("/")]
        rooms = [float(x) for x in c[5].split("/")]
        cons = [float(x) for x in c[6].split("/")]
        out.append(dict(label=label, own=own, rooms=rooms, cons=cons))
    return out


def main():
    rows = read_rows()
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.0), sharey=True)
    y = np.arange(len(rows))[::-1]
    wait_c, birth_c = "#9ecae1", "#08519c"
    panels = [("Rooms", "rooms", None), ("Consumption", "cons", None), ("Owns at end of period", "own", (0, 1))]
    for ax, (title, key, xlim) in zip(axes, panels):
        for yi, r in zip(y, rows):
            w, b = r[key]
            ax.plot([w, b], [yi, yi], color="0.75", lw=2, zorder=1)
            ax.scatter([w], [yi], color=wait_c, s=60, zorder=2, label="wait" if yi == y[0] else None)
            ax.scatter([b], [yi], color=birth_c, s=60, zorder=3, label="have the child" if yi == y[0] else None)
        if key == "rooms":
            ax.axvline(6.0, color="0.4", ls="--", lw=1)
            ax.text(6.05, y[-1] - 0.45, "rental cap", fontsize=9, color="0.3")
        ax.set_title(title, loc="left", fontsize=12)
        ax.grid(axis="x", alpha=0.3)
        if xlim:
            ax.set_xlim(*xlim)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([r["label"] for r in rows], fontsize=10)
    axes[0].set_ylim(y[-1] - 0.7, y[0] + 0.7)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, frameon=False, fontsize=11, bbox_to_anchor=(0.6, 1.06))
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"birth_vs_wait.{ext}", dpi=200, bbox_inches="tight")
    print(f"wrote {OUT / 'birth_vs_wait.pdf'}")


if __name__ == "__main__":
    main()
