"""Execute exact prepared reference/J restoration, forbidding every model solve."""
import argparse,json,sys,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--plan',required=True);p.add_argument('--output',required=True);a=p.parse_args()
root=Path(__file__).resolve().parents[7];sys.path.insert(0,str(root/'code/model/experiments/transition_readiness'))
import one_shock_floor as c
from floor_runtime import FloorRuntime
plan=json.loads(Path(a.plan).read_text());c.preflight(plan);out=Path(a.output)
rt=FloorRuntime.from_handoff(plan['handoff'],out/'runtime');adapter=c.NativeAdapter(rt,plan);controller=c.Controller(plan,adapter,out)
try:
    with rt.native_budget(time.monotonic()+240,0):
        controller.prepare()
    assert rt.total_native_calls==0 and controller.policy_calls==0
    assert rt.reference_verified and controller.seed['matrix'].shape==(24,24)
    c.write(out/'restore_receipt.json',dict(status='actual_reference_and_measured_J_restored',native_calls=0,reference_verified=True,matrix_shape=[24,24],identity=rt.identity(),original_preparation=plan['prepared_native_inputs'],controller=plan['source_files']['controller'],runtime=plan['source_files']['runtime']))
    print(json.dumps(dict(status='actual_reference_and_measured_J_restored',native_calls=0)))
except BaseException as exc:
    c.write(out/'restore_failure.json',dict(status='FAILED',exception=type(exc).__name__,reason=str(exc),native_calls=rt.total_native_calls))
    raise
