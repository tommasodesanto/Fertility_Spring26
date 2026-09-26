#!/usr/bin/env python3
"""Re-export the existing estate diagnostic figures with compact income legends.

Torch only. Uses saved solutions; no equilibrium or household solve is called.
Frozen source and original figures remain untouched.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pickle
from pathlib import Path
from unittest.mock import patch

import run_e5f_estate_receiver_probe as probe


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    probe.require_torch()
    if args.output.exists():
        raise FileExistsError(args.output)
    complete = json.loads((args.run_root / "complete.json").read_text())
    if complete["status"] != "exact_three_case_loop_complete":
        raise RuntimeError("completed numerical comparison required")
    args.output.mkdir(parents=True)
    _, _, _, _, _, runtime, _ = probe.verified_runtime(args.output)
    adapter = probe.load_adapter()
    adapter.install(runtime["model"], args.output / "estate_adapter")
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    import numpy as np

    original_legend, original_save = Axes.legend, Figure.savefig
    checks = []

    def curves(ax):
        digest = hashlib.sha256()
        for line in ax.get_lines():
            for values in (line.get_xdata(), line.get_ydata()):
                digest.update(np.asarray(values).tobytes())
        return digest.hexdigest()

    def legend(ax, *positional, **keywords):
        handles, labels = ax.get_legend_handles_labels()
        if not positional and labels and all(label.startswith("z=") for label in labels):
            before = curves(ax)
            labels = [f"z={float(label[2:]):.3f}" for label in labels]
            keywords.update(fontsize=8, handlelength=1.5, columnspacing=.6)
            result = original_legend(ax, handles, labels, **keywords)
            if curves(ax) != before:
                raise RuntimeError("legend formatting altered plotted data")
            ax._estate_legend_formatted = True
            return result
        return original_legend(ax, *positional, **keywords)

    def save(fig, path, *positional, **keywords):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for ax in fig.axes:
            if not getattr(ax, "_estate_legend_formatted", False):
                continue
            box = ax.get_window_extent(renderer)
            legend_box = ax.get_legend().get_window_extent(renderer)
            if (legend_box.x0 < box.x0-1 or legend_box.y0 < box.y0-1
                    or legend_box.x1 > box.x1+1 or legend_box.y1 > box.y1+1):
                raise RuntimeError("formatted legend exceeds its axes: " + str(path))
            checks.append(dict(figure=Path(path).name, curve_sha256=curves(ax),
                               legend_within_axes=True))
        return original_save(fig, path, *positional, **keywords)

    cases = {}
    for case in probe.CASES:
        source = args.run_root / case
        output = args.output / case
        output.mkdir()
        checkpoint = source / "initial_state.pkl.gz"
        with gzip.open(checkpoint, "rb") as stream:
            packet = pickle.load(stream)
        before = len(checks)
        with patch.object(Axes, "legend", legend), patch.object(Figure, "savefig", save):
            runtime["audit"].standard_diagnostics(packet, output, validate_production_young=False)
        expected = sorted(p.name for p in (source/"standard_diagnostics").glob("*.png"))
        actual = sorted(p.name for p in (output/"standard_diagnostics").glob("*.png"))
        if expected != actual or len(actual) != 17 or len(checks) == before:
            raise RuntimeError("figure set changed or no income legends were checked")
        cases[case] = dict(checkpoint_sha256=probe.sha(checkpoint), figure_count=17,
                           formatted_legends=len(checks)-before)
    probe.write(args.output / "figure_formatting.json", dict(
        status="complete_same_figures_legend_formatting", native_ge_solves=0,
        changes="income legend labels to three decimals; font 8, handle length 1.5, column spacing .6",
        cases=cases, checks=checks, original_figures_untouched=True))


if __name__ == "__main__":
    main()
