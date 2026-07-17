"""Re-solve a saved record and regenerate benchmark_lifecycle_compare.png (no title)."""

from __future__ import annotations
import sys, json, argparse
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "code/model"))
sys.path.insert(0, str(REPO / "code/model/tools"))

from make_full_dispersion_diagnostic_note import build_setup, solve_record, make_grid  # noqa: E402
from make_direct_fit_slide_note import plot_benchmark_lifecycle_compare  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--hR-max", type=float, default=8.0)
    args = parser.parse_args()

    record = json.loads(args.record.read_text())
    setup = build_setup(record, args.hR_max)
    sol, P, _ = solve_record(record, setup)
    b_grid = make_grid(P)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    plot_benchmark_lifecycle_compare(sol, P, b_grid, args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
