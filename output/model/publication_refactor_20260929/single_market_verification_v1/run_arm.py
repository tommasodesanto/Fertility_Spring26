"""One fixed-D14 full closure arm; copied reviewed loop, six-case/1200s cap."""
from pathlib import Path
from types import SimpleNamespace
import argparse, json, os, platform, sys, time
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE/'source')]
import driver as authored
import single_price
import phase_b_ge as ge
from small_credit_lab.contract import validate_contract

class ReplayBudget(authored.Budget):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs);self.max_lifecycle=6
        self.progress("initialized", max_lifecycle=6, stage_deadline_seconds=300)
    def claim_lifecycle(self,label):
        if self.smoke:raise RuntimeError('Smoke cannot claim lifecycle work')
        if self.used_lifecycle>=6:raise RuntimeError('Six-lifecycle cap reached')
        if time.time()+1>=self.deadline_epoch:raise RuntimeError('1200-second global deadline reached')
        super().claim_lifecycle(label)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['local-preflight','smoke','full']);ap.add_argument('--reference-root',type=Path,required=True);ap.add_argument('--bundle',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--deadline-epoch',type=float,required=True)
    args=ap.parse_args();started=time.time();deadline=min(args.deadline_epoch,started+1200)
    if args.out.exists():raise SystemExit('Refusing existing output')
    args.out.mkdir(parents=True)
    budget=ReplayBudget(args.out,deadline,smoke=args.mode!='full')
    try:
        context=authored.context_from_bundle(args);validate_contract(context['P']);context['deadline_epoch']=deadline
        context['selected_d_bar']=.14;context['reference_psi']=float(context['P'].psi_child)
        authored.write_json(args.out/'input_identity.json',context['loaded'].identity)
        authored.write_json(args.out/'runtime.json',dict(python=sys.version,numpy=np.__version__,node=platform.node(),platform=platform.platform(),threads={k:os.getenv(k) for k in ['OMP_NUM_THREADS','NUMBA_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']},global_seconds=1200,deadline_epoch=deadline,D=.14,wealth_nodes=len(context['b_grid']),income_nodes=len(context['P'].z_grid)))
        if args.mode=='local-preflight':
            context.update(fp=SimpleNamespace(write=authored.write_json),prepared=object(),manifest={},objective={},runtime=object(),reference={})
        else:authored.authenticate_frozen(context)
        if args.mode!='full':
            (args.out/'phase_b_ge').mkdir()
            # Actual GE root loop with a synthetic seed and mocked model calls.
            # The full seed call matches the verified grid-control driver.
            assert budget.used_lifecycle==0
            result=ge.smoke_phase_b(context,budget)
            assert result['lifecycle_solves']==0 and budget.used_lifecycle==0
            # Verify cap and deadline rejection without calling the model.
            gate=ReplayBudget(args.out/'gate_checks',deadline,smoke=False);gate.used_lifecycle=6
            try:gate.claim_lifecycle('must_reject');raise AssertionError('case cap was not enforced')
            except RuntimeError:pass
            gate.used_lifecycle=0;gate.deadline_epoch=time.time()-1
            try:gate.claim_lifecycle('must_reject');raise AssertionError('deadline was not enforced')
            except RuntimeError:pass
        else:
            seed=single_price.solve_fixed_price(context,.14,context['q_ref'],budget,'phase_a_fixed_d14',args.out/'seed')
            result=ge.run_phase_b(context,dict(selected_d_bar=.14,selected_live=seed),budget)
            if result['status']!='passed':raise RuntimeError('Full closure not certified within budget')
            if budget.used_lifecycle!=6:raise RuntimeError('Matched reference six-case loop changed')
        authored.write_json(args.out/'completed.json',dict(result=result,lifecycle_solves=budget.used_lifecycle,input_identity=context['loaded'].identity,budget_gates="six-case and expired-deadline rejection passed" if args.mode!='full' else "enforced",preflight_observer='mocked_local_not_authenticated' if args.mode=='local-preflight' else 'frozen_authenticated',elapsed_seconds=time.time()-started))
    except BaseException as exc:
        authored.write_json(args.out/'failure.json',dict(status='failed',error=str(exc),type=type(exc).__name__,lifecycle_solves=budget.used_lifecycle,no_auto_retry=True));raise
if __name__=='__main__':main()
