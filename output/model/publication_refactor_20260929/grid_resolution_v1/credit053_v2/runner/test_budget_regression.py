"""Mocked controller regression; no Bellman, KFE or native observer calls."""
import csv,json,tempfile,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import run_comparison as runner
ge=runner.ge

class MockBudget:
    stage_deadline_seconds=300
    def __init__(self,out):self.out=out;self.deadline_epoch=time.time()+2400;self.used_lifecycle=1;self.max_lifecycle=20
    @property
    def remaining_lifecycle(self):return self.max_lifecycle-self.used_lifecycle
    def claim_lifecycle(self,label):
        if self.used_lifecycle>=20:raise RuntimeError('20-call cap')
        self.used_lifecycle+=1

def main():
    prior_solve,prior_observe=ge.solve_fixed_price,ge.observe_price
    with tempfile.TemporaryDirectory(dir=runner.HERE) as tmp:
        out=Path(tmp)
        def stage(q):return dict(price=np.array([q]),sol=SimpleNamespace(a=np.array([q])),sd=SimpleNamespace(b=np.array([1.])),case_deadline_epoch=time.time()+300)
        for converge in (True,False):
            folder=out/str(converge);folder.mkdir();budget=MockBudget(folder);labels=[];roots=[]
            context=dict(fp=SimpleNamespace(write=runner.write),prepared=object(),manifest={},objective={},runtime=object(),reference={},P=SimpleNamespace(psi_child=.1355551166583114),q_ref=1.,out=folder,deadline_epoch=budget.deadline_epoch)
            def solve(ctx,d,q,b,label,path):
                assert d==.53;b.claim_lifecycle(label);labels.append(label)
                if label.startswith('root_'):roots.append(label)
                return stage(q)
            def observe(ctx,live,label,*,final=False):
                q=float(live['price'][0]);res=(0. if converge and (len(roots)>=8 or final or 'repeat' in label) else (-.01 if q<1.02 else .01))
                if final:
                    directory=folder/'phase_b_ge'/label;directory.mkdir(parents=True)
                    for name,count,key in [('target_fit.csv',14,'moment'),('parameters.csv',31,'parameter')]:
                        with (directory/name).open('w',newline='') as f:
                            w=csv.DictWriter(f,fieldnames=[key,'value']);w.writeheader();w.writerows({key:str(i),'value':'1'} for i in range(count))
                return dict(price=q,renewal_residual=res,population_scale=1.03)
            ge.solve_fixed_price,ge.observe_price=solve,observe
            try:result=ge.run_phase_b(context,dict(selected_d_bar=.53,selected_live=stage(1.)),budget)
            finally:ge.solve_fixed_price,ge.observe_price=prior_solve,prior_observe
            if converge:
                assert result['status']=='passed' and len(roots)>=8 and budget.used_lifecycle>6 and 'selected_repeat' in labels
                success=dict(mock_counted_calls=budget.used_lifecycle,phase_b_new_calls=len(labels),root_trials=len(roots),repeat_retained=True)
            else:
                assert result['status']=='uncomputed_bounded_budget' and budget.used_lifecycle==19 and budget.remaining_lifecycle==1 and 'selected_repeat' not in labels
                failure=dict(status=result['status'],mock_counted_calls=budget.used_lifecycle,repeat_reserve=budget.remaining_lifecycle)
            budget.used_lifecycle=20
            try:budget.claim_lifecycle('excess')
            except RuntimeError:pass
            else:raise AssertionError('20-call cap was not enforced')
        actual=runner.ArmBudget(out/'actual_budget',time.time()+2400,smoke=False)
        for i in range(20):actual.claim_lifecycle('root_%02d'%i)
        try:actual.claim_lifecycle('excess')
        except RuntimeError:pass
        else:raise AssertionError('Actual ArmBudget did not reject call21')
    receipt=dict(status='passed_mock_budget_regression',lifecycle_solves=0,extended_convergence=success,bounded_nonconvergence=failure,actual_20_call_cap_rejected_21=True,case_seconds=300,total_seconds=2400)
    (runner.HERE/'budget_regression_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
