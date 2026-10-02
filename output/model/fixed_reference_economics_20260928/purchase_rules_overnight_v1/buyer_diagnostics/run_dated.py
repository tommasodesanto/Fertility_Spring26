"""Financial-access addendum from an authenticated dated diagnostic packet."""
from __future__ import annotations

import argparse
import gzip
import json
import pickle
import sys
from pathlib import Path

import numpy as np

from financial_access import matched_first_birth_access, write_compact_summary

ROOT = Path(__file__).resolve().parents[5]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine-root", type=Path, required=True,
                    help="Matching isolated engine root containing refactor_lab")
    ap.add_argument("--packet", type=Path, required=True,
                    help="Trusted selected dated diagnostic_packet.pkl.gz")
    ap.add_argument("--rule", choices=("hard", "quarter"), required=True)
    ap.add_argument("--observed-phi", choices=(0.8, 1.0), type=float, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit("Refusing existing dated access output")
    sys.path[:0] = [str(args.engine_root.resolve()), str(ROOT / "code/model/tools"), str(ROOT / "code/model"),
                    str(ROOT / "code/model/experiments/transition_readiness")]
    # The native policy bundle class is defined here; import before unpickling.
    import run_dynamic_population_transition  # noqa: F401
    with gzip.open(args.packet, "rb") as stream:
        packet = pickle.load(stream)
    required = {"parameters", "b_grid", "evaluation", "shared", "period"}
    if not isinstance(packet, dict) or not required.issubset(packet):
        raise RuntimeError("Dated packet lacks exact native objects")
    P = packet["parameters"]
    grid = np.asarray(packet["b_grid"], float)
    ev = packet["evaluation"]
    policy = ev.policy
    price = np.asarray(policy.price, float)
    if price.ndim != 1 or np.any(price <= 0) or not np.isfinite(price).all():
        raise RuntimeError("Invalid dated price")
    result = matched_first_birth_access(
        P, packet["shared"], grid, price,
        np.asarray(ev.g_pre, float), np.asarray(policy.fert_probs, float),
        rule=args.rule,
        observed_tenure_probs=np.asarray(policy.tenure_probs, float),
        observed_phi=args.observed_phi)
    write_compact_summary(args.out, result)
    audit = result["observed_choice_audit"]
    if audit["share_outside_map"] is not None and audit["share_outside_map"] > 1e-7:
        raise RuntimeError("Dated observed owner choices fall outside pointwise financial map")
    receipt = args.out.with_suffix(".receipt.json")
    receipt.write_text(json.dumps({
        "status": "dated_financial_access_diagnostic_passed",
        "period": int(packet["period"]),
        "rule": args.rule,
        "observed_phi": args.observed_phi,
        "packet": str(args.packet.resolve()),
        "summary": str(args.out.resolve()),
        "no_model_solve": True,
        "financial_access_only": True,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
