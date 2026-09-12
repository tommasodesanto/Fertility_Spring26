"""Build the September housing/fertility schematics with preference labels.

This wrapper reuses the canonical schematic builder without changing its source
or numerical construction.  It intercepts the builder's save calls, relabels
only the legend text from ``v`` to ``\\psi``, and writes the three requested
PDFs.
"""

from pathlib import Path
import importlib.util
import tempfile

import matplotlib.pyplot as plt


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
            stem = Path(path).stem.replace("housing_fertility_stage_", "september_housing_fertility_stage_")
            original_savefig(figure, OUTPUT / f"{stem}.pdf", *args, **kwargs)

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
