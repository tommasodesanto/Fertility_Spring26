#!/usr/bin/env python3
"""Collect strict, exactly repeated E5 fixed-kappa fertility probes."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from intergen_eqscale_seq_optimized.collect_e1 import eligible_tight


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULTS_ROOT = ROOT / "output/model/eqscale_seq_e5_probe_20260725/production"
DEFAULT_OUTDIR = ROOT / "output/model/eqscale_seq_e5_probe_20260725/report"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    paths = sorted(args.results_root.glob("chain_*/summary.json"))
    if not paths:
        raise RuntimeError(f"no E5 probe chain summaries under {args.results_root}")
    frontier: list[dict[str, Any]] = []
    chain_rows: list[dict[str, Any]] = []
    for path in paths:
        summary = json.loads(path.read_text(encoding="utf-8"))
        meta = summary.get("metadata") or {}
        valid = (meta.get("arm") == "E5_PROBE_KE" and int(meta.get("free_parameter_count", -1)) == 9
                 and int(meta.get("target_count", -1)) == 12 and "probe_fixed_kappa_fert" in meta)
        if not valid:
            raise RuntimeError(f"{path} does not satisfy the E5 fixed-kappa probe contract")
        tight = eligible_tight(summary)
        chain_rows.append({"summary_path": str(path), "seed": meta.get("seed"), "eligible": tight is not None,
                           "probe_fixed_kappa_fert": meta["probe_fixed_kappa_fert"],
                           "search_cases": summary.get("n_cases_completed"), "search_strict": summary.get("n_strict"),
                           "tight_loss": tight.get("rank_loss") if tight else None,
                           "tight_residual": tight.get("market_residual") if tight else None})
        if tight is None:
            continue
        models = {row["moment"]: row.get("model") for row in tight.get("target_fit", [])}
        row = {"probe_fixed_kappa_fert": meta["probe_fixed_kappa_fert"], "tight_loss": tight["rank_loss"],
               **models}
        theta = tight.get("theta") or {}
        for domain_row in meta.get("active_domain") or []:
            name = domain_row["name"]
            row[name] = theta.get("beta" if name == "beta_annual" else name)
        psi = float(theta.get("psi_child", float("nan")))
        row["psi_at_bound"] = abs(psi) > 0.98 * float(meta["psi_bound"])
        frontier.append(row)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "frontier.csv", frontier)
    write_csv(args.outdir / "chain_summary.csv", chain_rows)


if __name__ == "__main__":
    main()
