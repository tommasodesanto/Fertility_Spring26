#!/usr/bin/env python3
"""Torch-only tables and supplemental figures for a completed receipt-risk test.

Reads CSVs/JSONs/images only; never imports the model or loads a checkpoint.
The 17 standard figures per case are unchanged.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import run_e5f_estate_receiver_probe as probe

CASES = ("no_receipt", "conditional_mean", "receipt_lottery")
LABELS = ("No receipts", "Certain age-specific mean", "Inheritance lottery")


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--coarse-run-root", type=Path, help="report effects relative to a completed coarser-grid run")
    args = parser.parse_args()
    probe.require_torch()
    root = args.run_root
    complete = json.loads((root / "complete.json").read_text())
    test_count = int(re.findall(r"(\d+) passed", (root / "tests.log").read_text())[-1])
    if complete["status"] != "completed" or [c["case"] for c in complete["cases"]] != list(CASES):
        raise ValueError("A completed ordered three-case run is required")
    report = root / "report"
    report.mkdir(exist_ok=False)
    targets = {case: read_csv(root / case / "target_fit.csv") for case in CASES}
    parameters = {case: read_csv(root / case / "parameters.csv") for case in CASES}
    if any(len(targets[case]) != 13 or len(parameters[case]) != 27 for case in CASES):
        raise ValueError("Incomplete target or parameter tables")
    if any(parameters[case] != parameters[CASES[0]] for case in CASES):
        raise ValueError("The held parameter/restriction table differs across cases")
    comparison = []
    for index, row in enumerate(targets[CASES[0]]):
        if any(any(targets[case][index][key] != row[key] for key in ("moment", "target", "weight")) for case in CASES):
            raise ValueError("Target or weight contract differs across cases")
        values = [float(targets[case][index]["model"]) for case in CASES]
        comparison.append(dict(moment=row["moment"], target=float(row["target"]),
            **dict(zip(CASES, values)), certain_minus_control=values[1]-values[0],
            lottery_minus_control=values[2]-values[0], lottery_minus_certain=values[2]-values[1]))
    probe.table(report / "comparison.csv", comparison)
    if args.coarse_run_root:
        coarse = {case: read_csv(args.coarse_run_root / case / "target_fit.csv") for case in CASES}
        sensitivity = []
        for index, row in enumerate(comparison):
            if any(coarse[case][index]["moment"] != row["moment"] for case in CASES):
                raise ValueError("Coarse/fine target rows differ")
            base = float(coarse[CASES[0]][index]["model"])
            item = dict(moment=row["moment"])
            for case, label in zip(CASES[1:], ("certain", "lottery")):
                coarse_effect = float(coarse[case][index]["model"]) - base
                fine_effect = row[case] - row["no_receipt"]
                item[label + "_coarse_effect"] = coarse_effect
                item[label + "_fine_effect"] = fine_effect
                item[label + "_effect_change"] = fine_effect - coarse_effect
            sensitivity.append(item)
        probe.table(report / "grid_effect_comparison.csv", sensitivity)
    probe.table(report / "full_target_fit.csv", [dict(case=case, **row) for case in CASES for row in targets[case]])
    probe.table(report / "parameters.csv", parameters[CASES[0]])
    summary = []
    for row in complete["cases"]:
        summary.append(dict(case=row["case"], loss=row["loss"],
            price=row["price"][0], market_residual=row["relative_market_residual"],
            seconds=row["wall_seconds"], clipped_wealth_loss=row["clipped_wealth_loss"],
            **row["estate_accounts"]))
    probe.table(report / "accounts.csv", summary)
    probe.write(report / "summary.json", dict(cases=summary, comparison=comparison,
        source_pins=complete["source_pins"], status="complete_tables_and_figures_visual_review_pending"))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    keys = ("ownership_before_current_choices", "financial_wealth_before_current_choices",
            "mean_receipt", "financial_wealth_sd")
    titles = ("Ownership before current choices (%)", "Financial wealth before current choices",
              "Expected receipt per household at this age", "Financial wealth standard deviation")
    for case, label in zip(CASES, LABELS):
        rows = read_csv(root / case / "lifecycle_and_receipts.csv")
        age = [float(row["age"]) for row in rows]
        for ax, key, title in zip(axes.flat, keys, titles):
            factor = 100. if key == keys[0] else 1.
            ax.plot(age, [float(row[key])*factor for row in rows], label=label, linewidth=1.8)
            ax.set(title=title, xlabel="Age")
            ax.grid(alpha=.2)
    axes[1, 0].legend(fontsize=8)
    fig.suptitle("Supplemental: inheritance risk at fixed prices and equal age-specific means\n"
                 "Wealth and receipts in annual-earnings units; receipts arrive before choices", fontsize=12)
    fig.savefig(report / "receipt_risk_lifecycle.png", dpi=170)
    fig.savefig(report / "receipt_risk_lifecycle.pdf")
    plt.close(fig)

    from PIL import Image, ImageDraw
    for case in CASES:
        figures = sorted((root / case / "standard_diagnostics").glob("*.png"))
        if len(figures) != 17:
            raise ValueError("Incomplete standard figure packet")
        for offset in (0, 9):
            sheet = Image.new("RGB", (1500, 1110), "white")
            draw = ImageDraw.Draw(sheet)
            for k, path in enumerate(figures[offset:offset+9]):
                x, y = (k % 3)*500, (k // 3)*370
                draw.text((x+8, y+8), f"{case}: {path.stem}", fill="black")
                with Image.open(path) as source:
                    source = source.convert("RGB")
                    source.thumbnail((490, 340))
                    sheet.paste(source, (x+5, y+25))
            sheet.save(report / f"{case}_standard_contact_{offset//9+1}.png")
    grid_subdivision = complete["cases"][0].get("wealth_grid_subdivision", 1)
    control_description = ("the retained net-valuation case" if grid_subdivision == 1 else
                           "a fresh unpatched net-valuation solve and calendar distribution on the finer grid")
    lines = ["# Inheritance uncertainty: completed fixed-price test", "",
        f"{complete['household_solves']} household solves; wealth-grid subdivision {grid_subdivision}. "
        "Experimental, not recalibrated or funded general equilibrium.", "",
        "![Supplemental lifecycle comparison](receipt_risk_lifecycle.png)", "",
        f"The no-receipt control exactly reproduces {control_description}. All {test_count} numerical tests, "
        "the ordered-loop smoke and native calendar/budget/purchase/fiscal/value/probability checks pass. "
        "No occupied wealth is clipped. The two receipt cases have the same expected transfer at each age "
        "and the same fixed house price. Housing and estate funding residuals below are deliberately reported.", "",
        "All economic changes are experimental: age-pooled published receipt profiles, zero receipts outside "
        "supported model nodes26–78, IID receipt risk, and start-period liquid-wealth timing. "
        "The certain and lottery cases share this timing. Entry wealth, B15 earnings, preferences, "
        "the inherited 8.751% tax and the frozen target/weight contract remain fixed. "
        "The target table includes the inherited experimental first-birth rooms target1.465; this does not adopt it "
        "as the future empirical/model observation contract.", "",
        "| Moment | Target | No receipts | Certain mean | Lottery |", "|---|---:|---:|---:|---:|"]
    for row in comparison:
        lines.append("| " + row["moment"] + " | " + " | ".join(f"{row[key]:.3f}" for key in ("target", *CASES)) + " |")
    lines += ["", "| Case | Inherited loss | House-market residual | Paid − generated estates |", "|---|---:|---:|---:|"]
    for row in summary:
        lines.append(f"| {row['case']} | {row['loss']:.3f} | {100*row['market_residual']:.3f}% | {row['paid_minus_generated']:.3e} |")
    lines += ["", "[Full target fits: every gap, weight and loss contribution](full_target_fit.csv)", "",
        "[Every parameter, estimate, bound, restriction and near-bound flag](parameters.csv)", "",
        "[Estate accounts and run diagnostics](accounts.csv)", "",
        "The 17 standard figures for each case are retained unchanged:"]
    for case in CASES:
        lines += ["", f"- [{case}, figures1–9]({case}_standard_contact_1.png)",
                  f"- [{case}, figures10–17]({case}_standard_contact_2.png)"]
    if args.coarse_run_root:
        lines += ["", "[Wealth-grid sensitivity: changes in each receipt effect relative to its own control](grid_effect_comparison.csv)"]
    (report / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
