#!/usr/bin/env python3
"""Zero-lifecycle checks of the actual three-arm controller with injected evaluation."""
import importlib.util,json,tempfile,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('share_A_controller',HERE/'runner.py');runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
def main():
    seen=[]
    def evaluator(state,arm,out,deadline):
        assert deadline<=state['end'];seen.append(arm)
        value=dict(lifecycle_seconds=0.,cohort_summary={'completed_fertility':float(len(seen))},renewal_residual_reported_not_imposed=.1,relative_market_residual_reported_not_imposed=.2)
        runner.write(out/'closure.json',value);return value
    with tempfile.TemporaryDirectory() as folder:
        end=time.time()+60
        records=runner.run(Path(folder)/'success',end,state={'end':end},evaluator=evaluator,mock=True)
        assert seen==list(runner.ARMS) and len(records)==3
        receipt=json.loads((Path(folder)/'success/completed.json').read_text());assert receipt['lifecycle_solves']==0 and receipt['case_attempts']==3
        def fail(*args):raise RuntimeError('injected accounting failure')
        try:runner.run(Path(folder)/'failure',end,state={},evaluator=fail,mock=True)
        except RuntimeError:pass
        else:raise AssertionError('Fatal case was swallowed')
        failed=json.loads((Path(folder)/'failure/completed.json').read_text());assert failed['case_attempts']==1 and failed['status']=='failed_no_retry'
    print(json.dumps(dict(status='passed',lifecycle_solves=0,checks=['exact_three_arm_order','zero_mock_lifecycle','600s_case_1200s_total_deadline','checkpoint_each_case','fatal_gate_stops_without_retry'])))
if __name__=='__main__':main()
