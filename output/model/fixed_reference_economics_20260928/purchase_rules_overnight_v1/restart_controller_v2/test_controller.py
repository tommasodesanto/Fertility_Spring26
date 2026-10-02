"""Zero-model-solve loop checks for the restart controller."""
import json
import sys
import tempfile
import time
from pathlib import Path

PACKET=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PACKET))
sys.argv=['run_psi.py','--chain','0']
import run_psi
from restart_controller_v2.controller import install_censored_optimizer,require_launcher_scaffold


def main():
    install_censored_optimizer(run_psi,8)
    with tempfile.TemporaryDirectory() as tmp:
        scaffold=Path(tmp)/'launch'
        scaffold.mkdir()
        (scaffold/'numba_cache').mkdir()
        (scaffold/'matplotlib').mkdir()
        (scaffold/'search.log').write_text('apptainer startup text\n')
        require_launcher_scaffold(scaffold)
        (scaffold/'restart_contract.json').write_text('{}\n')
        try: require_launcher_scaffold(scaffold)
        except AssertionError: pass
        else: raise AssertionError('Duplicate scientific results were accepted')
        (scaffold/'restart_contract.json').unlink()
        (scaffold/'numba_cache'/'stale-cache').write_text('x')
        try: require_launcher_scaffold(scaffold)
        except AssertionError: pass
        else: raise AssertionError('Nonempty launcher cache was accepted')
        out=Path(tmp)/'continue'
        n=[0]
        def evaluate(label,point,deadline):
            n[0]+=1
            if n[0]==2:
                return dict(status='budget_exhausted',reason='uncomputed_bounded_budget',lifecycle_solves=31)
            return dict(status='passed',residual=[float(n[0])]+[0.]*9,
                        lifecycle_solves=1,report=str(out/label))
        receipt=run_psi.optimize(out,{'x':.5},{'x':(0.,1.)},('x',),evaluate,time.time()+3600,toy=True)
        cases=json.loads((out/'cases.json').read_text())
        assert len(cases)>=3 and cases[0]['status']=='passed' and cases[1]['optimizer_only_censor']
        assert cases[1]['loss'] is None and cases[1]['objective'] is None
        assert cases[1]['computed_valid_loss'] is False
        assert cases[2]['status']=='passed' and receipt['selected']['status']=='passed'
        assert receipt['search_stop_reason']!='native_evaluation_budget_exhausted'
        late=run_psi.optimize(Path(tmp)/'clock',{'x':.5},{'x':(0.,1.)},('x',),evaluate,time.time()+1,toy=True)
        assert late['objective_calls']==0 and late['search_stop_reason']=='four_hour_actual_start_final_reserve'
        print(json.dumps(dict(status='passed_zero_model_solves',cases=len(cases),censored=1,
                              late_stop=late['search_stop_reason'])))


if __name__=='__main__': main()
