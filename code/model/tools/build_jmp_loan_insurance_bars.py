"""Bar chart for the JMP slides: loans versus insurance, from the markdown table
"Credit, transfers and price on one calibration" in output/model/rental_menu_precaution_20261007/tables_floor_credit.md (parsed at run time).
Left: completed fertility (% change, column "completed vs bench"). Right: ownership 18-29, pp change = 100*(row value - benchmark value).

Run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
     output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/tools/build_jmp_loan_insurance_bars.py
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "output/model/rental_menu_precaution_20261007/tables_floor_credit.md"
OUT = ROOT / "output/model/jmp_slides_policy_functions_20261007"
HEADING = "## Credit, transfers and price on one calibration"
# (label, prefix of the table's experiment cell), top to bottom
ROWS = [("House prices and rents \u221210%", None),   # not in the table: PRICE_DOWN below
        ("Financed share 80% to 95%", "LTV 80% -> 95%"),
        ("Renter credit, 1 year of earnings", "renter credit D = 0.52"),
        ("Renter credit, 2 years of earnings", "renter credit D = 1.03"),
        ("Insurance: transfer when earnings are low", "income-contingent transfer"),
        ("Same total cost, paid to everyone", "certain transfer, equal expected value")]
# Price and rent x0.9 (cells A6_P090 vs A6): not in tables_floor_credit.md. Computed from the saved solution arrays with
# outcomes() of output/model/rental_menu_precaution_20261007/analyze_floor_credit.py (same function as the table; it reproduces the table's
# A6 0.438 and A6_P110 0.409 ownership 18-29 and completed fertility 1.866 / 1.767): completed fertility 1.963842 vs 1.865698, ownership 0.468824.
PRICE_DOWN = dict(ceb=1.963842, ceb0=1.865698, own=0.468824)


def parse():
    lines = SRC.read_text().splitlines()
    k = lines.index(HEADING)
    tab = []
    for ln in lines[k + 1:]:
        if ln.startswith("#"):
            break
        if ln.startswith("|"):
            tab.append([c.strip() for c in ln.strip().strip("|").split("|")])
    head, body = tab[0], tab[2:]
    ic, io = head.index("completed vs bench"), head.index("ownership 18-29")
    bench = next(r for r in body if r[0] == "benchmark")
    own0 = float(bench[io])
    out = []
    for label, pref in ROWS:
        if pref is None:
            out.append((label, 100 * (PRICE_DOWN["ceb"] / PRICE_DOWN["ceb0"] - 1), 100 * (PRICE_DOWN["own"] - own0)))
            continue
        r = next(r for r in body if r[0].startswith(pref))
        out.append((label, float(r[ic].replace("%", "").replace("+", "")), 100 * (float(r[io]) - own0)))
    return out, own0


def main():
    data, own0 = parse()
    print("benchmark ownership 18-29:", own0)
    for d in data:
        print(f"{d[0]:38s} completed {d[1]:+.2f}%  ownership {d[2]:+.2f} pp")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "font.family": "sans-serif"})
    labels = [d[0] for d in data]
    y = list(range(len(data)))[::-1]       # first row on top
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, idx, ttl in [(axes[0], 1, "Completed fertility (% change)"), (axes[1], 2, "Ownership, ages 18-29 (pp change)")]:
        vals = [d[idx] for d in data]
        cols = ["#9ecae1" if v < 0 else "#2171b5" for v in vals]
        ax.barh(y, vals, color=cols, height=0.62)
        span = max(abs(v) for v in vals)
        for yy, v in zip(y, vals):
            ax.text(v + (0.02 * span if v >= 0 else -0.02 * span), yy, f"{v:+.1f}".replace("-", "\u2212"), va="center",
                    ha="left" if v >= 0 else "right", fontsize=10)
        lo, hi = min(0, min(vals)), max(0, max(vals))
        ax.set_xlim(lo - 0.18 * span, hi + 0.18 * span)
        ax.axvline(0, color="0.3", lw=0.9)
        ax.set_title(ttl, loc="left", fontsize=11)
        ax.grid(axis="x", alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(labels)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "loan_vs_insurance_bars.pdf")
    fig.savefig(OUT / "loan_vs_insurance_bars.png", dpi=200)
    print("wrote", OUT / "loan_vs_insurance_bars.pdf")


if __name__ == "__main__":
    main()
