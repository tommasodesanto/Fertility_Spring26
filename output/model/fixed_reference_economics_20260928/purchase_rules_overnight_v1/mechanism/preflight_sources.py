"""Import and bind each actual isolated engine without a household solve."""
from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace
import numpy as np

from selected_runtime import PACKET, ROOT, OLD, TRANSITION, load, read, require, sha


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--arm",choices=("hard","quarter"),required=True)
    args=parser.parse_args()
    for rel,digest in read(PACKET/"engine_pins.json").items():
        require(sha(PACKET/rel)==digest,"Engine source drift: "+rel)
    for rel,digest in read(PACKET/"source_pins.json").items():
        require(sha(ROOT/rel)==digest,"Integration source drift: "+rel)
    engine_root=PACKET/"engines"/args.arm
    sys.path.insert(0,str(engine_root))
    from small_credit_lab.engine import solver
    from refactor_lab.engine import solver as check_solver
    require(Path(solver.__file__).resolve().is_relative_to(engine_root)
            and Path(check_solver.__file__).resolve().is_relative_to(engine_root),
            "Imported engine escaped isolated arm")
    sys.path[:0]=[str(OLD),str(ROOT/"code/model/tools"),str(ROOT/"code/model"),str(TRANSITION)]
    import inputs
    import runner
    from floor_runtime import FloorRuntime
    sys.path.insert(0,str(runner.BASE))
    base=load("purchase_base_preflight",runner.BASE/"run_comparison.py")
    ge=load("purchase_ge_preflight",OLD.parent/"utility_calibration_round1_v1"/"phase_b_pilot.py")
    runner.install_reporter_on_authored(base.authored)
    P,grid=inputs.proposal("floor_s0")
    P,entry=inputs.entry(P,grid,"nonnegative_mean")
    P.experimental_purchase_saving_fraction=.25 if args.arm=="quarter" else 1.
    with tempfile.TemporaryDirectory() as tmp:
        ctx=base.authored.context_from_bundle(SimpleNamespace(
            bundle=ROOT/"output/model/publication_refactor_20260929/local_export_v1/inputs",
            reference_root=ROOT,out=Path(tmp)))
        base.authored.authenticate_frozen(ctx)
        rt=FloorRuntime()
        rt.model=solver;rt.P=P;rt.grid=grid;rt.ctx=ctx
        rt.rt=ctx["prepared"].rt;rt.pf=rt.rt["primitive"].pf
        rt.install_observer_adapters(Path(tmp)/"observers")
        with rt.native_bindings():
            cohort=rt.pf.calendar.entrant_cohort(np.array([1.]),P,grid)
            np.testing.assert_allclose(
                cohort.sum(axis=(1,2,4,5)),
                P.fixed_reference_entry_conditional*P.z_weights[None,:],
                rtol=0,atol=2e-16)
        print(dict(status="source_and_native_callback_import_passed_zero_solves",
                   arm=args.arm,engine=str(Path(solver.__file__).resolve()),
                   observer=str(Path(ge.__file__).resolve()),grid=[int(P.Nb),int(P.Nz)],
                   native_policy_calls=0,entry_arm=entry["arm"]))


if __name__=="__main__":main()
