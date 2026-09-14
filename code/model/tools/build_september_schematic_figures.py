"""Build the September housing/fertility schematics with preference labels.

Reuse the canonical curves and endpoints, adding the author's fixed-price
fertility response between A and B in the right panel only. Write PDF slides
and PNG previews without modifying the shared canonical builder.
"""

from pathlib import Path
import importlib.util
import tempfile

import matplotlib.pyplot as plt
from matplotlib.text import Annotation


ROOT = Path(__file__).resolve().parents[3]
SOURCE = Path(__file__).with_name("build_fertility_population_housing_transition_figures.py")
OUTPUT = ROOT / "latex" / "figures"


def _load_source():
    spec = importlib.util.spec_from_file_location("canonical_schematics", SOURCE)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load canonical builder: {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    source = _load_source()
    values = source.illustration()
    source.verify(values)
    OUTPUT.mkdir(parents=True, exist_ok=True)

    original_savefig = plt.Figure.savefig
    original_close = source.plt.close
    with tempfile.TemporaryDirectory(prefix="september-schematic-") as temp_dir:
        source.OUTDIR = Path(temp_dir)

        def savefig_with_psi_labels(figure, path, *args, **kwargs):
            if Path(path).suffix.lower() != ".pdf":
                return
            for text in figure.findobj(match=plt.Text):
                text.set_text(text.get_text().replace("v_0", r"\psi_0").replace("v_1", r"\psi_1"))
                if text.get_text() == r"$I$":
                    text.set_text(r"$B$")
                elif text.get_text() == r"$A'$":
                    text.set_text(r"$C$")
            stage = Path(path).stem.rsplit("_", 1)[-1]
            # Label replacement fertility as n-bar rather than the normalized 1.
            if stage == "initial":
                figure.axes[1].set_yticks([values["replacement"]], labels=[r"$\bar n$"])
            else:
                figure.axes[1].set_yticks(
                    [values["impact_fertility"], values["replacement"]],
                    labels=[r"$\widetilde n$", r"$\bar n$"],
                )
            if stage in {"impact", "adjustment"}:
                child_axis = figure.axes[1]
                point_a = (values["old_price"], values["replacement"])
                point_b = (values["impact_price"], values["impact_fertility"])
                for label in child_axis.texts:
                    if label.get_text() == r"$B$":
                        label.set_position((point_b[0] + 0.012, point_b[1] - 0.085))
                    elif label.get_text() == r"$C$":
                        label.set_position((values["new_price"] - 0.052,
                                            values["replacement"] + 0.105))
                direct_fertility = float(source.child_demand(
                    values["old_price"], values["new_theta"],
                    values["child_cost"], values["child_space"],
                ))
                point_aprime = (values["old_price"], direct_fertility)
                assert direct_fertility < point_b[1] < point_a[1]
                # Replace only the former direct A-to-I arrow in the right panel.
                arrows = [text for text in child_axis.texts
                          if isinstance(text, Annotation) and tuple(text.xy) == point_b]
                assert len(arrows) == 1
                arrows[0].remove()
                child_axis.scatter(*point_aprime, color="#a3473f", s=40, zorder=6)
                child_axis.text(point_aprime[0] + 0.012, point_aprime[1] - 0.085,
                                r"$A'$", color="#a3473f")
                arrow_color = "#555555" if stage == "impact" else "#999999"
                source.add_arrow(child_axis, point_a, point_aprime, color=arrow_color)
                source.add_arrow(child_axis, point_aprime, point_b, color=arrow_color)
            stem = Path(path).stem.replace("housing_fertility_stage_", "september_housing_fertility_stage_")
            original_savefig(figure, OUTPUT / f"{stem}.pdf", *args, **kwargs)
            original_savefig(figure, OUTPUT / f"{stem}.png", dpi=180, pad_inches=0.05)

        source.plt.Figure.savefig = savefig_with_psi_labels
        source.plt.close = lambda *args, **kwargs: None
        try:
            for stage in ("initial", "impact", "adjustment"):
                source.build_housing_fertility_stage_figure(
                    values,
                    stage=stage,
                    output_stem=f"housing_fertility_stage_{stage}",
                )
        finally:
            source.plt.Figure.savefig = original_savefig
            source.plt.close = original_close
            plt.close("all")


if __name__ == "__main__":
    main()
