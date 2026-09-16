#!/usr/bin/env python3
"""Compare two sandbox runs' moments.csv and parameters.csv.

    PYTHONPATH=. python sandbox/compare.py baseline kappa_h_zero

Writes compare.md and one overlay PDF of key policy plots to
output/model/sandbox/compare_<A>_vs_<B>/.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

SANDBOX_ROOT = Path(__file__).resolve().parent
MODEL_ROOT = SANDBOX_ROOT.parent
REPO_ROOT = MODEL_ROOT.parents[1]
DEFAULT_OUT_ROOT = REPO_ROOT / "output/model/sandbox"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def compare_moments(a: list[dict[str, str]], b: list[dict[str, str]]) -> list[str]:
    b_by_moment = {row["moment"]: row for row in b}
    lines = ["| Moment | A | B | B - A |", "|---|---:|---:|---:|"]
    for row in a:
        moment = row["moment"]
        other = b_by_moment.get(moment)
        try:
            va = float(row["model"])
            vb = float(other["model"]) if other else float("nan")
            diff = f"{vb - va:.6g}"
            va_s, vb_s = f"{va:.6g}", f"{vb:.6g}"
        except (TypeError, ValueError):
            va_s, vb_s, diff = row["model"], (other["model"] if other else "-"), "-"
        lines.append(f"| {moment} | {va_s} | {vb_s} | {diff} |")
    return lines


def compare_parameters(a: list[dict[str, str]], b: list[dict[str, str]]) -> list[str]:
    b_by_param = {row["parameter"]: row for row in b}
    lines = ["| Parameter | A | B | B - A |", "|---|---:|---:|---:|"]
    for row in a:
        name = row["parameter"]
        if name == "_solved_price":
            continue
        other = b_by_param.get(name)
        try:
            va = float(row["estimate"])
            vb = float(other["estimate"]) if other else float("nan")
            diff = f"{vb - va:.6g}"
            va_s, vb_s = f"{va:.6g}", f"{vb:.6g}"
        except (TypeError, ValueError):
            va_s, vb_s, diff = row["estimate"], (other["estimate"] if other else "-"), "-"
        lines.append(f"| {name} | {va_s} | {vb_s} | {diff} |")
    return lines


def overlay_policy_pdf(name_a: str, moments_a: list[dict[str, str]], name_b: str,
                        moments_b: list[dict[str, str]], out_pdf: Path) -> bool:
    """Bar-chart overlay of the two runs' scalar moments (from moments.csv).

    moments.csv holds one solved-model number per row, not raw policy-function
    arrays (the four-file output contract keeps arrays out of moments.csv), so
    this overlays the moment table rather than reconstructing policy-function
    curves. Each run's own graphs.pdf carries the full policy-function pages.
    """
    import math
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    b_by_moment = {row["moment"]: row for row in moments_b}
    labels, values_a, values_b = [], [], []
    for row in moments_a:
        try:
            va = float(row["model"])
            vb = float(b_by_moment[row["moment"]]["model"])
        except (KeyError, TypeError, ValueError):
            continue
        if not (math.isfinite(va) and math.isfinite(vb)):
            continue
        labels.append(row["moment"])
        values_a.append(va)
        values_b.append(vb)
    if not labels:
        return False

    with PdfPages(out_pdf) as pdf:
        chunk = 6
        for start in range(0, len(labels), chunk):
            sl = slice(start, start + chunk)
            fig, ax = plt.subplots(figsize=(10, 6))
            x = range(len(labels[sl]))
            ax.barh([xi - 0.2 for xi in x], values_a[sl], height=0.4, label=name_a)
            ax.barh([xi + 0.2 for xi in x], values_b[sl], height=0.4, label=name_b)
            ax.set_yticks(list(x))
            ax.set_yticklabels(labels[sl], fontsize=8)
            ax.legend()
            ax.set_title(f"Moment overlay: {name_a} vs {name_b}")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name_a")
    parser.add_argument("name_b")
    parser.add_argument("--root", default=str(DEFAULT_OUT_ROOT))
    args = parser.parse_args()

    root = Path(args.root)
    dir_a, dir_b = root / args.name_a, root / args.name_b
    for d in (dir_a, dir_b):
        if not (d / "moments.csv").exists():
            raise FileNotFoundError(f"Missing run output: {d}/moments.csv (run sandbox/run_ss.py first)")

    moments_a, moments_b = read_csv(dir_a / "moments.csv"), read_csv(dir_b / "moments.csv")
    params_a, params_b = read_csv(dir_a / "parameters.csv"), read_csv(dir_b / "parameters.csv")

    out_dir = root / f"compare_{args.name_a}_vs_{args.name_b}"
    out_dir.mkdir(parents=True, exist_ok=True)

    lines = [f"# Compare {args.name_a} vs {args.name_b}", "", "## Moments", ""]
    lines += compare_moments(moments_a, moments_b)
    lines += ["", "## Parameters", ""]
    lines += compare_parameters(params_a, params_b)
    (out_dir / "compare.md").write_text("\n".join(lines) + "\n")

    made_pdf = overlay_policy_pdf(args.name_a, moments_a, args.name_b, moments_b, out_dir / "overlay.pdf")
    print(f"Wrote {out_dir}/compare.md" + (f" and {out_dir}/overlay.pdf" if made_pdf else " (overlay.pdf skipped: no comparable finite moments)"))


if __name__ == "__main__":
    main()
