#!/usr/bin/env python3
"""Three-case cap-eight model sensitivity, never a production calibration.

Use the certified 239-node baseline and the exact inherited 2023 measure. Only
the continuous rental upper bound changes; the owner choice set is unchanged.
"""
from __future__ import annotations
import argparse
import copy
import json
import threading
import time
from pathlib import Path
import numpy as np
import run_e5f_global_saving_quantification as saving

audit = saving.audit
baseline, calendar, transition = saving.baseline, saving.calendar, saving.transition
FINE_BASELINE_SHA = "04c48693586124c97f2f33cdb4a5a199db3e9605ff92335eae0fcf7a13fd8cb7"
AUDIT_SHA = "ba0f94f52705f43fd0bf7a18ea8c2053f9897674abeb20b030eaf37c57de8623"
SAVING_SHA = "0e89d55a6cd3d0432dcff39b24683a7f6233ac0820a03b47a94eda8e2a924586"
CAP = 8.0
MARKET_TOL = 1e-5


def parameter_record(P):
    return {k:calendar.jsonable(v) for k,v in vars(P).items() if not k.startswith("_")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint",type=Path,required=True)
    parser.add_argument("--outdir",type=Path,required=True)
    parser.add_argument("--stage",choices=("smoke","policies"),required=True)
    parser.add_argument("--smoke-summary",type=Path)
    parser.add_argument("--smoke-summary-sha256")
    args = parser.parse_args();out = args.outdir.resolve()
    if out.exists() and any(out.iterdir()): raise FileExistsError(out)
    out.mkdir(parents=True,exist_ok=True)
    for path,expected in ((args.checkpoint,FINE_BASELINE_SHA),(Path(audit.__file__),AUDIT_SHA),(Path(saving.__file__),SAVING_SHA)):
        if audit.digest(path) != expected: raise RuntimeError(f"Pinned input differs: {path}")
    if args.stage == "policies":
        if not args.smoke_summary or audit.digest(args.smoke_summary) != args.smoke_summary_sha256:
            raise RuntimeError("Exact-loop cap-eight smoke receipt missing or mismatched")
        receipt = json.loads(args.smoke_summary.read_text())
        if receipt.get("status") != "complete" or receipt.get("stage") != "smoke" or receipt.get("driver_sha256") != audit.digest(__file__):
            raise RuntimeError("This driver has no passed cap-eight smoke")
    transition.configure_sequential_model()
    saving.install_sequential_operators()
    packet = audit.load_checkpoint(args.checkpoint)
    packet["checkpoint_sha256"] = FINE_BASELINE_SHA
    _,configured = transition.configure_sequential_model()
    live = baseline.calibration.code_fingerprint_contract(configured)["bundle_sha256"]
    if live != audit.BUNDLE or packet["input_hashes"]["code_bundle_sha256"] != live:
        raise RuntimeError("Frozen production code differs")
    if len(packet["b_grid"]) != 239 or float(packet["parameters"].hR_max) != 6.0:
        raise RuntimeError("The certified fine-grid cap-six baseline is required")
    before = parameter_record(packet["parameters"])
    packet["parameters"] = copy.deepcopy(packet["parameters"])
    packet["parameters"].hR_max = CAP
    after = parameter_record(packet["parameters"])
    changed = [k for k in set(before)|set(after) if before.get(k) != after.get(k)]
    if changed != ["hR_max"]: raise RuntimeError(f"Unexpected economic parameter change: {changed}")
    contract = dict(stage=args.stage,driver_sha256=audit.digest(__file__),source_bundle_sha256=live,
        audit_helper_sha256=AUDIT_SHA,global_saving_helper_sha256=SAVING_SHA,
        inherited_fine_baseline_checkpoint_sha256=FINE_BASELINE_SHA,
        smoke_summary_sha256=args.smoke_summary_sha256,operator_contract=saving.operator_contract(),
        economic_parameter_changes={"hR_max":dict(reference=6.0,diagnostic=CAP)},
        estimated_parameters="all eleven unchanged at task_010; no re-estimation",
        external_restrictions="supply elasticity 0.63 and tenure taste dispersion 0.005 unchanged",
        rental_cap_classification="explicit diagnostic economic override; eight rooms is not an empirical estimate",
        owner_choices="unchanged 2,4,6,8,10 room options and unchanged taste distribution",
        population="unchanged exact inherited 2023 pre-choice measure on nested 239-node grid",
        entry_and_geography="inherited queue and single-market closure unchanged; no forward population path evaluated",
        fiscal="existing baseline, supply +20%, and dependent-child LTV 95% contracts; no tax/rebate changes",
        outstanding="history and target fit under cap eight; upper wealth bound, long-run and welfare effects; other housing menus",
        production_unchanged=True)
    audit.save_json(out/"run_contract.json",contract)
    audit.save_json(out/"parameter_override_receipt.json",dict(before=before,after=after,changed=changed))
    audit.save_json(out/"latest_completed_case.json",dict(status="no_case_complete"))
    audit.save_json(out/"best_so_far.json",dict(status="diagnostic_model_variant_no_calibration_selection"))
    started = time.perf_counter();stop = threading.Event()
    def heartbeat():
        while not stop.wait(60):
            audit.save_json(out/"runtime_heartbeat.json",dict(elapsed_seconds=time.perf_counter()-started,
                utc=time.strftime("%Y-%m-%d %H:%M:%S UTC",time.gmtime())))
    threading.Thread(target=heartbeat,daemon=True).start()
    original_solve = calendar.solve_policy
    def traced_solve(price,P,bg,shared,counter):
        t = time.perf_counter();result = original_solve(price,P,bg,shared,counter)
        audit.progress(out,"bellman_complete",rental_cap=float(P.hR_max),nodes=len(bg),
            price=np.asarray(price).tolist(),seconds=time.perf_counter()-t,count=counter.bellman)
        return result
    calendar.solve_policy = traced_solve
    names = ("baseline",) if args.stage == "smoke" else ("supply-plus-20","dependent-child-ltv95")
    rows = []
    try:
        for name in names:
            row = saving.run_case(packet,"global",name,True,out,market_tol=MARKET_TOL)
            case_dir = Path(row["diagnostic_directory"])
            arrays = json.loads((case_dir/"policy_array_summary.json").read_text())
            gates = dict(requested_market=row["market_residual"] <= MARKET_TOL,
                occupied_value_monotonicity=arrays["occupied_negative_steps"] == 0,
                probabilities=all(p["nonfinite"] == 0 and p["minimum"] >= -1e-7 and p["maximum"] <= 1+1e-7 for p in arrays["probabilities"].values()),
                standard_graphs=len(list((case_dir/"standard_diagnostics").glob("*.png"))) == 17)
            row.update(rental_cap=CAP,nodes=239,all_gates_passed=all(gates.values()))
            audit.save_json(case_dir/"cap_sensitivity_receipt.json",dict(gates=gates,row=row,
                interpretation="Cap-eight economic variant, not a numerical repair or a calibrated model. Budget excess remains separately reported."))
            audit.save_json(out/"latest_completed_case.json",row)
            audit.save_json(out/"best_so_far.json",dict(status="diagnostic_model_variant_no_calibration_selection",latest_admissible_case=row if all(gates.values()) else None))
            if not all(gates.values()): raise RuntimeError(f"Cap sensitivity failed; evidence retained: {gates}")
            rows.append(row);baseline.write_csv(out/"case_summary.csv",rows)
        audit.save_json(out/"summary.json",dict(contract,status="complete",all_case_gates_passed=True,rows=rows,elapsed_seconds=time.perf_counter()-started))
        audit.progress(out,"stage_complete",stage=args.stage,cases=len(rows),elapsed_seconds=time.perf_counter()-started)
    except Exception as exc:
        audit.save_json(out/"failure.json",dict(contract,status="failed",error=repr(exc),completed_rows=rows))
        raise
    finally:
        stop.set()


if __name__ == "__main__":
    main()
