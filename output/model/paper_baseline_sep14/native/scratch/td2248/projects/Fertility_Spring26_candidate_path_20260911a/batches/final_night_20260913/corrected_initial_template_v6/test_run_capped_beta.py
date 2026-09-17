import contextlib, copy, io, math, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import run_capped_beta as c


def restrictions():
    d={n:dict(lower=.01,upper=20.) for n in c.ALL_NAMES}
    d['beta_annual']=dict(lower=.94,upper=.99)
    return d

def item(beta=.99): return dict(parameters={n:(beta if n=='beta_annual' else 1.) for n in c.ALL_NAMES},initial_psi=.2)

def score(q):
    # beta's optimum is strictly inside the cap; psi/fixed rows catch accidental parameter leakage.
    v=np.array([q['parameters'][n]-(.97 if n=='beta_annual' else 1.) for n in c.ALL_NAMES])
    r=np.r_[v,np.zeros(3)]
    fits=[dict(restriction_id='normalization',target=2.1,model=2.1,gap=0.,scored=False)]
    fits += [dict(restriction_id=str(i),target=0.,model=float(x),gap=float(x),actual_weight=1.,loss_contribution=float(x*x),scored=True) for i,x in enumerate(r)]
    rows=[dict(parameter=n,estimate=q['parameters'][n],lower=.01,upper=20.,status='estimated') for n in c.ALL_NAMES]
    rows += [dict(parameter='psi_child',estimate=q['initial_psi']+.01,status='derived'),dict(parameter='fixed_x',estimate=3.,status='fixed')]
    return dict(contract_sha256=c.APPROVED_OBJECTIVE,free_parameter_count=9,loss=float(r@r),target_fit=fits,parameters=rows,normalization=dict(target=2.1,psi_child=q['initial_psi']+.01,completed_fertility=2.1,absolute_gap=0.,status='verified'))

def plan(seconds=1e6):
    q=item(); return dict(resume_proposal=dict(q,case_id='pinned'),workers=18,rounds=3,maximum_search_cases=90,search_seconds=seconds)

class CappedBetaTests(unittest.TestCase):
    def execute(self, mutate=None, seconds=1e6):
        calls=[]; p=plan(seconds); seed=score(p['resume_proposal'])
        def worker(q):
            calls.append(copy.deepcopy(q)); s=score(q)
            ans=dict(case_id=q['case_id'],proposal=q,status='verified',score=s,loss=s['loss'],output='/synthetic',receipt_status='verified_scored_candidate')
            if q.get('repetitions')==2: ans['second_signature_equal']=True; ans['repetitions']=2
            if mutate: ans=mutate(ans)
            return ans
        with tempfile.TemporaryDirectory() as t,contextlib.redirect_stdout(io.StringIO()): summary,best=c.search(p,seed,restrictions(),worker,Path(t))
        return summary,best,calls

    def test_three_round_free_beta_search(self):
        s,b,calls=self.execute(); self.assertEqual(s['status'],'completed'); self.assertLessEqual(s['search_attempted_cases'],90); self.assertEqual(len(calls),s['search_attempted_cases']+3)
        self.assertLess(b['score']['parameters'][0]['estimate'],.99)
        self.assertEqual(calls[-1]['parameters'].keys(), item()['parameters'].keys())

    def test_bound_beta_has_only_inward_derivative_and_all_j_columns(self):
        seed=dict(case_id='x',proposal=item(),score=score(item()))
        ps=c.derivative_proposals(seed,restrictions(),0); beta=[x for x in ps if x['column']==0]
        self.assertEqual(len(beta),1); self.assertLess(beta[0]['parameters']['beta_annual'],.99)
        J=c.jacobian([dict(status='verified',proposal=dict(x,step=x['step']),score=score(x)) for x in ps],seed)
        self.assertEqual(J.shape,(12,9))

    def test_above_cap_rejected_before_worker(self):
        q=item(1.); self.assertRaises(ValueError,c.validate_proposal,q,restrictions())

    def test_seed_and_final_repetition_failures_stop(self):
        s,_,calls=self.execute(lambda q: dict(q,receipt_status='bad') if q['case_id']=='seed_exact_repeat_01' else q)
        self.assertEqual(s['status'],'stopped_for_review'); self.assertEqual(len(calls),1)
        s,_,_=self.execute(lambda q: dict(q,second_signature_equal=False) if q.get('repetitions')==2 else q)
        self.assertEqual(s['status'],'stopped_for_review')

    def test_no_time_keeps_verified_seed_and_final_repeats(self):
        s,_,calls=self.execute(seconds=1); self.assertEqual(s['status'],'completed'); self.assertEqual(s['stop_reason'],'stage_deadline_or_case_budget'); self.assertEqual(len(calls),3)

    def test_no_fabricated_jacobian_column(self):
        with self.assertRaises(RuntimeError): c.derivative_column([],np.zeros(12))

    def test_deadline_is_rechecked_after_derivative_stage(self):
        now=[0.]
        def advance(q):
            if q['case_id'].startswith('r0_d'):now[0]=7800.
            return q
        with patch.object(c.time,'monotonic',side_effect=lambda:now[0]):
            s,_,calls=self.execute(advance,seconds=7800)
        self.assertEqual(s['status'],'completed')
        self.assertEqual(s['stop_reason'],'stage_deadline_or_case_budget')
        self.assertFalse(any('_joint_' in q['case_id'] for q in calls))
        self.assertTrue(s['selected_exact_repetitions_verified'])

    def test_output_beta_is_free_with_actual_cap_and_near_bound(self):
        q=item(.9899);b=dict(score=score(q),output='/synthetic')
        with tempfile.TemporaryDirectory() as t:
            c.save_tables(Path(t),b,restrictions())
            s=c.read(Path(t)/'selected_capped_score.json')
        beta=next(r for r in s['parameters'] if r['parameter']=='beta_annual')
        self.assertEqual(s['free_parameter_count'],9)
        self.assertEqual(beta['status'],'estimated_capped')
        self.assertEqual((beta['lower'],beta['upper']),(.94,.99))
        self.assertTrue(beta['near_bound'])

    def test_mass_rejection_is_excluded_then_recurrence_stops(self):
        def one(q):
            return dict(q,status='rejected_mass_gate') if q['case_id']=='r0_d1_-1' else q
        s,_,_=self.execute(one); self.assertEqual(s['status'],'completed')
        def two(q):
            return dict(q,status='rejected_mass_gate') if q['case_id'] in ('r0_d1_-1','r0_d2_-1') else q
        s,_,_=self.execute(two); self.assertEqual(s['status'],'stopped_for_review')

if __name__=='__main__': unittest.main()
