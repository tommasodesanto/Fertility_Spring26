#!/usr/bin/env python3
"""Collect and adjudicate the four E1 promotion-task checkpoints."""
from __future__ import annotations
import argparse, json
from pathlib import Path
def main() -> None:
    p=argparse.ArgumentParser(); p.add_argument("--run-root",type=Path,required=True); p.add_argument("--output",type=Path,default=None); a=p.parse_args(); records={}
    for task in ("demand-map","root-case","strict-repeat","throughput"):
        path=a.run_root / f"{task}.json"
        if path.exists(): records[task]=json.loads(path.read_text())
    demand=records.get("demand-map",{}).get("points",[]); root=records.get("root-case",{}); repeat=records.get("strict-repeat",{}); throughput=records.get("throughput",{})
    summary={"completed_tasks":sorted(records),"demand_map_pass":bool(demand) and all(x["sign_equal"] for x in demand),"root_case_pass":bool(root) and root.get("forced_wrong_direction_exercised") and root["ordinary"]["market_residual"] <= 2.5e-5 and root["forced_wrong_direction"]["market_residual"] <= 2.5e-5 and root["forced_wrong_direction"].get("directional_fallback_used"),"strict_repeat_pass":bool(repeat) and repeat.get("price_bitwise_equal") and repeat.get("loss_bitwise_equal") and repeat.get("moments_bitwise_equal"),"throughput_pass":len(throughput.get("candidates",[])) in (3,10) and "ge_loose_total_seconds" in throughput,"records":records}
    out=a.output or a.run_root / "promotion_summary.json"; out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n"); print(out)
if __name__ == "__main__": main()
