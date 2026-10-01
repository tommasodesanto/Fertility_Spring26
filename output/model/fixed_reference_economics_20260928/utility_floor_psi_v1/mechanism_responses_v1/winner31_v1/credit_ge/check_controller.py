"""Zero-lifecycle tests of signed q0 routing, budgets, and selected repeats."""
import math,shutil,tempfile,time
from pathlib import Path
from types import SimpleNamespace
import run_ge as ge

plan_path=Path(__file__).parent/'plan.json'
plan=ge.verify_plan(plan_path,zeroLC=True)
winner=ge.ROOT/'output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/ROOT'

def check(root_factor,expected_direction,expected_first):
    with tempfile.TemporaryDirectory(prefix='winner31_ge_mock_') as temporary:
        folder=Path(temporary);q0=folder/'fixed_q0';q0.mkdir()
        for name in ('target_fit.csv','parameters.csv'):shutil.copyfile(winner/name,q0/name)
        ge.write(q0/'receipt.json',dict(status='completed_support_limited_diagnostic',regime='lifetime_repayment_only',price_factor=1.0,source_binding_sha256=ge.sha(ge.MECHANISM/'source_binding.json'),target_fit_sha256=ge.sha(q0/'target_fit.csv'),parameters_sha256=ge.sha(q0/'parameters.csv')))
        ge.write(q0/'closure.json',dict(candidate_base_loss=31.28400725566496,candidate_price_q0=.719168368828958,price=.719168368828958,candidate_psi_child=.17156192800028292,grid_nodes=120,standard_plot_count=17,support_diagnostic=dict(status='occupied_support_pass_unoccupied_alternatives_unverified'),natural_support_certified=False,renewal_residual_reported_not_imposed=math.log(root_factor)))
        args=SimpleNamespace(out=folder/'ge',plan=plan_path,deadline_epoch=time.time()+2400,q0_credit_case=q0)
        seen=[]
        def fake(name,factor,role,case_deadline,out):
            assert case_deadline<=args.deadline_epoch and case_deadline<=time.time()+301
            if role!='repeat':assert case_deadline<=args.deadline_epoch-300+1e-6
            seen.append((name,factor,role))
            receipt=dict(status='passed_support_limited_case',lifecycle_solves=1,plan_sha256=ge.sha(args.plan),price=plan['candidate_q0']*factor,renewal_residual=math.log(root_factor/factor),closure=dict(population_scale=1.),checkpoint=dict(path='mock'),standard_plot_count=17,standard_plot_sha256={f'p{i}.png':'mock' for i in range(17)},native_scaled_step=dict(status='passed'),repeat_check=dict(status='passed',standard_plot_hashes_exact=17))
            (out/name).mkdir();ge.write(out/name/'receipt.json',receipt)
            return receipt
        result=ge.controller(args,plan,evaluate=fake,mock=True)
        assert result['status']=='completed_support_limited_stationary_ge_diagnostic'
        assert result['direction']==expected_direction and seen[0][0]==expected_first
        assert seen[-1][0]=='selected_repeat' and len(seen)<=11
        assert result['total_lifecycle_attempts']==len(seen)
        assert all(.8<=factor<=1.35 for _,factor,_ in seen)
        if expected_direction!='root_at_q0':assert not any(name=='q0_selected' for name,_,_ in seen)
        return len(seen)

counts={
    'lower':check(.87,'lower','lower_95'),
    'upper':check(1.13,'upper','upper_105'),
    'root_at_q0':check(1.,'root_at_q0','q0_selected')}
with tempfile.TemporaryDirectory(prefix='winner31_ge_bad_identity_') as temporary:
    folder=Path(temporary);q0=folder/'fixed_q0';q0.mkdir()
    shutil.copyfile(winner/'parameters.csv',q0/'parameters.csv')
    target_text=(winner/'target_fit.csv').read_text().replace('cps_childlessness','wrong_target_name',1)
    (q0/'target_fit.csv').write_text(target_text)
    ge.write(q0/'receipt.json',dict(status='completed_support_limited_diagnostic',regime='lifetime_repayment_only',price_factor=1.0,source_binding_sha256=ge.sha(ge.MECHANISM/'source_binding.json'),target_fit_sha256=ge.sha(q0/'target_fit.csv'),parameters_sha256=ge.sha(q0/'parameters.csv')))
    ge.write(q0/'closure.json',dict(candidate_base_loss=31.28400725566496,candidate_price_q0=.719168368828958,price=.719168368828958,candidate_psi_child=.17156192800028292,grid_nodes=120,standard_plot_count=17,support_diagnostic=dict(status='occupied_support_pass_unoccupied_alternatives_unverified'),natural_support_certified=False,renewal_residual_reported_not_imposed=-.1))
    args=SimpleNamespace(out=folder/'ge',plan=plan_path,deadline_epoch=time.time()+2400,q0_credit_case=q0)
    try:ge.controller(args,plan,evaluate=lambda *_: (_ for _ in ()).throw(AssertionError('lifecycle dispatched')),mock=True)
    except RuntimeError as error:assert 'target identity' in str(error)
    else:raise AssertionError('Altered target identity was accepted')
print('PASS zero-LC signed directions, reviewed domains, fresh q0/root repeats, 11-new-attempt and 300/2400 deadlines',counts)
