#!/usr/bin/env python3
"""Select the primary A3 winner without choosing the external LTV variant by loss."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--primary-ltv", type=float, default=0.4)
    parser.add_argument("--primary-theta1", type=float, default=0.25)
    args = parser.parse_args()

    candidates: list[dict[str, Any]] = []
    for summary_path in sorted(args.results_root.glob("task_*/summary.json")):
        payload = json.loads(summary_path.read_text())
        meta = dict(payload.get("metadata", {}))
        # Only the twice re-solved tight winner is eligible; the loose-search
        # record is retained for optimizer diagnostics but never reported as a
        # calibration result.
        best = payload.get("best_tight")
        mechanism = dict(meta.get("mechanism", {}))
        fixed = dict(meta.get("fixed_parameters", {}))
        if not isinstance(best, dict) or meta.get("arm") != "A3":
            continue
        if not math.isclose(float(mechanism.get("owner_ltv_terminal_share", math.nan)), args.primary_ltv):
            continue
        if not math.isclose(float(fixed.get("theta1", math.nan)), args.primary_theta1):
            continue
        if not bool(best.get("strict_converged", False)):
            continue
        best = dict(best)
        best["source_summary"] = str(summary_path)
        candidates.append(best)

    if not candidates:
        raise SystemExit("no strict primary A3 candidate was produced")
    accepted = [row for row in candidates if bool(row.get("lifecycle", {}).get("hard_acceptance_pass", False))]
    pool = accepted or candidates
    winner = min(pool, key=lambda row: float(row.get("rank_loss", math.inf)))
    out = {
        "selection_rule": {
            "arm": "A3",
            "primary_ltv_terminal_multiplier": args.primary_ltv,
            "primary_theta1": args.primary_theta1,
            "external_variants_not_loss_selected": True,
            "hard_acceptance_required_if_any_chain_passes": True,
            "winner_passes_hard_acceptance": bool(winner.get("lifecycle", {}).get("hard_acceptance_pass", False)),
        },
        "candidate_count": len(candidates),
        "accepted_candidate_count": len(accepted),
        "best": winner,
        "theta": winner["theta"],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True))
    print(args.out)


if __name__ == "__main__":
    main()
