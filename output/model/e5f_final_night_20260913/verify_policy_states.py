"""Torch read-only native comparison of saved policy inputs; no model solve."""
import argparse,gzip,hashlib,json,os,pickle,sys,time
from pathlib import Path
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[name]='1'
import numpy as np

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--case-dir',type=Path,required=True)
    ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--helper',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    read=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    manifest=read(args.manifest);plan=read(manifest['prior_plan']);root=Path(plan['source_root'])
    sys.path[:0]=[str(args.helper),str(root/'code/model/tools'),str(root/'code/model')]
    import e5f_rebated_surprises as rebated
    from run_e5f_final_rebated_history import verify_pins,_same_state,save,clean
    verify_pins(manifest['file_sha256']);verify_pins(plan['file_sha256']);rebated._runtime()
    args.out.parent.mkdir(parents=True,exist_ok=True);assert not args.out.exists()
    started=time.time();receipt=dict(status='running',no_model_solve=True)
    try:
        case=args.case_dir;contract=read(case/'contract_receipt.json');fits=read(case/'realized_fit.json')
        assert contract['manifest_sha256']==sha(args.manifest) and len(fits)==4
        conditional=contract.get('conditional_history_count')
        if conditional is not None:assert conditional==6 and contract['history_refitted'] is False
        complete='conditioning_history_complete.json' if conditional is not None else 'finite_history_complete.json'
        assert read(case/complete)['realized']==fits
        state=case/'realized_state_2023.pkl.gz'
        with gzip.open(state,'rb') as f:inherited=pickle.load(f)
        assert inherited.year==2023
        heads=inherited.households.g_pre;parameters=[];grids=[];rules=[];quantities=[];prices=[];root_hashes={};sources=[state,args.manifest,case/'contract_receipt.json',Path(__file__)]
        for name,tax in [('baseline_rebate',.01),('tax2_rebate',.02)]:
            folder=case/'policies'/name;r=read(folder/'root_receipt.json');s=read(folder/'summary.json')
            assert s['finite_converged'] and s['annual_tax']==tax and r['converged'] and r['finite_horizon_market_fiscal_converged']
            assert r['start_year']==2023 and r['count']==contract['count'] and r['case']==contract['case'] and r['psi']==fits[-1]['psi']
            assert r['final']['mapping_valid'] and r['final_reproduction_max_abs']<=2e-10
            snapshot=folder/'first_period_diagnostics.pkl.gz'
            with gzip.open(snapshot,'rb') as f:pack=pickle.load(f)
            evaluation=pack['evaluation'];P=pack['parameters'];np.testing.assert_array_equal(evaluation.g_pre,heads)
            assert evaluation.g_pre.dtype==heads.dtype
            assert P.tau_H==4*tax and P.tau_pay==.179 and P.psi_child==fits[-1]['psi']
            assert P.user_cost_rate==P.R_gross+P.delta+P.tau_H-1.
            grids.append(pack['b_grid']);rules.append(pack['supply_rule']);prices.append(float(evaluation.policy.price[0]));quantities.append(float(evaluation.supply_by_loc[0]))
            np.testing.assert_array_equal(pack['supply_rule'].quantity(evaluation.policy.price),evaluation.supply_by_loc)
            parameters.append(P);root_hashes[name]=sha(folder/'root_receipt.json');sources.extend([folder/'summary.json',folder/'root_receipt.json',snapshot])
        np.testing.assert_array_equal(grids[0],grids[1]);assert _same_state(rules[0],rules[1]), 'Different housing supply rules'
        left,right=map(vars,parameters);assert left.keys()==right.keys()
        differences=[key for key in left if not _same_state(left[key],right[key])]
        allowed={'tau_H','user_cost_rate','pension','pension_by_loc','property_tax_lump_sum_transfer','income'}
        assert set(differences)<=allowed,('Unexpected parameter change',differences)
        np.testing.assert_array_equal(parameters[0].income[:,:parameters[0].J_R],parameters[1].income[:,:parameters[1].J_R])
        receipt.update(status='passed',conditional_history_count=conditional,common_initial_g_pre=True,common_grid=True,common_supply_rule=True,supply_rule_fields=clean(vars(rules[0])),
            prices=prices,housing_supply=quantities,implied_supply_elasticity=float(np.log(quantities[1]/quantities[0])/np.log(prices[1]/prices[0])),initial_g_pre_sha256=hashlib.sha256(heads.tobytes()).hexdigest(),
            initial_shape=list(heads.shape),initial_dtype=str(heads.dtype),case_contract_sha256=sha(case/'contract_receipt.json'),
            policy_root_sha256=root_hashes,allowed_parameter_differences=sorted(allowed),actual_parameter_differences=differences,
            all_other_parameter_fields_exact=True,worker_income_exact=True,input_sha256={str(p):sha(p) for p in sources},
            horizon_verified=False,production_eligible=False,elapsed_seconds=time.time()-started)
        verify_pins(manifest['file_sha256']);save(args.out,receipt);print(json.dumps(dict(status='passed',differences=differences)))
    except BaseException as exc:
        receipt.update(status='failed',error_type=type(exc).__name__,error=str(exc),elapsed_seconds=time.time()-started)
        save(args.out,receipt);raise
if __name__=='__main__':main()
