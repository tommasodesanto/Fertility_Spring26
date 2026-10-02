"""Authenticate one selected overnight fit, then tabulate saved policies only."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from financial_access import matched_first_birth_access, write_compact_summary

PACKET = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PACKET / "mechanism"))
import selected_runtime  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", choices=("hard", "quarter"), required=True)
    ap.add_argument("--completed", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit("Refusing existing diagnostic directory")
    args.out.mkdir(parents=True)
    rt = selected_runtime.construct(args.arm, args.completed, args.out / "runtime_auth")
    _, _, _, repeat = selected_runtime.authenticate_selected(args.arm, args.completed)
    stage = repeat / "stage/solution_arrays.npz"
    from refactor_lab.engine.parameters import get_fecundity_by_age
    with np.load(stage, allow_pickle=False) as saved:
        required = ("distribution.g_pre", "g_beginning_distribution", "fert_probs",
                    "tenure_probs", "b_grid", "p_eq")
        missing = [key for key in required if key not in saved]
        if missing:
            raise RuntimeError("Selected stage lacks required arrays: " + ", ".join(missing))
        gpre = saved["distribution.g_pre"]
        gpost = saved["g_beginning_distribution"]
        fert = saved["fert_probs"]
        tenure = saved["tenure_probs"]
        np.testing.assert_array_equal(saved["b_grid"], rt.grid)
        np.testing.assert_allclose(saved["p_eq"], [rt.reference_price], atol=1e-12, rtol=0)
        if gpre.shape != gpost.shape or abs(float(gpre.sum() - gpost.sum())) > 1e-8:
            raise RuntimeError("Pre/post fertility mass shape or total differs")
    P = rt.P
    SD = rt.model.precompute_shared(P, rt.grid)
    fec = get_fecundity_by_age(P)
    if len(fec) != gpre.shape[3]:
        raise RuntimeError("Selected fecundity vector length mismatch")
    ltv_out = args.out / "buyer_net_closing_ratio.json"
    cmd = [sys.executable, str(Path(__file__).with_name("summarize_saved_buyers.py")),
           "--stage", str(stage), "--pre-birth", str(stage),
           "--fecundity-by-age", ",".join(format(float(x), ".17g") for x in fec),
           "--out", str(ltv_out),
           "--owner-rungs", ",".join(format(float(x), ".17g") for x in P.H_own),
           "--sale-cost", format(float(P.psi), ".17g"),
           "--age-start", format(float(P.age_start), ".17g"),
           "--age-step", format(float(P.da), ".17g")]
    if args.arm == "hard":
        cmd.extend(("--hard-phi", ".8"))
    subprocess.run(cmd, check=True)
    access = matched_first_birth_access(
        P, SD, rt.grid, np.array([rt.reference_price]), gpre, fert,
        rule=args.arm, observed_tenure_probs=tenure, observed_phi=0.8)
    write_compact_summary(args.out / "matched_first_birth_financial_access.json", access)
    audit = access["observed_choice_audit"]
    if audit["share_outside_map"] is not None and audit["share_outside_map"] > 1e-7:
        raise RuntimeError("Observed first-birth owner choices fall outside financial feasibility map")
    (args.out / "completed.json").write_text(json.dumps({
        "status": "saved_policy_buyer_diagnostics_passed",
        "arm": args.arm,
        "selected_completed": str(args.completed.resolve()),
        "stage": str(stage.resolve()),
        "ltv_summary": str(ltv_out.resolve()),
        "access_summary": str((args.out / "matched_first_birth_financial_access.json").resolve()),
        "no_model_solve": True,
        "gross_mortgage_ltv_identified": False,
        "continuation_value_feasibility_not_checked": True,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
