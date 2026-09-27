#!/usr/bin/env python3
"""Compare Torch daytime smokes to the selected Mac case; explicitly promote."""
import argparse,csv,hashlib,json,pathlib
p=argparse.ArgumentParser();p.add_argument('--root',type=pathlib.Path,required=True);p.add_argument('--approve',action='store_true');a=p.parse_args();root=a.root

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text())
def write(p,x):
    with p.open('x') as f: json.dump(x,f,sort_keys=True,indent=2);f.write('\n')
def rows(p,key):
    with p.open() as f: return {r[key]:r for r in csv.DictReader(f)}
contract=root/'contract.json';c=read(contract);smoke=root/'run/smoke/complete.json';r=read(smoke)
assert sha(contract)=='dc5d49ecf05d81b00ee5c515fca222b6de7935fa72e64d7345bca5d48e3fdb04'
assert r['status']=='exact_loop_smoke_passed' and r['contract_sha256']==sha(contract)
assert len(r['records'])==2 and all(x['status']=='success' for x in r['records'])
comparisons=[]
for record in r['records']:
    case=pathlib.Path(record['case_path']);receipt=read(case/'receipt.json')
    assert receipt['point']==c['initial_point']
    assert receipt['case_checkpoint_sha256']==sha(case/'initial_state.pkl.gz')
    assert len(list((case/'standard_diagnostics').glob('*.png')))==17
    for file,key,numeric in [('target_fit.csv','moment',['target','model','gap','weight','loss_contribution']),('parameters.csv','parameter',['estimate','lower','upper'])]:
        old=rows(root/'input'/('selected_'+file),key);new=rows(case/file,key)
        assert old.keys()==new.keys()
        for name in old:
            for field in numeric:
                x,y=old[name][field],new[name][field]
                if x==y: continue
                assert x and y,(name,field)
                difference=float(y)-float(x)
                assert abs(difference)<=1e-5,(name,field,difference)
                comparisons.append(dict(case=record['case'],table=file,row=name,field=field,old=float(x),new=float(y),difference=difference))
    assert receipt['market_residual']<=2e-4
    assert abs(receipt['normalization']['completed_fertility']-2.1)<=5e-4
result=dict(status='cross_host_smokes_passed',rows_checked={'targets_per_case':14,'parameters_per_case':31},tolerance=1e-5,differences=comparisons,smoke_sha256=sha(smoke),smoke_contract_sha256=sha(contract))
if not (root/'cross_host_comparison.json').exists():write(root/'cross_host_comparison.json',result)
if a.approve:
    # The caller has inspected smoke evidence; only approval metadata may change.
    production=root/'production_contract.json';c['status']='approved_production';c['production_blockers']=[]
    c['verified_smoke']={'path':str(smoke),'sha256':sha(smoke)}
    write(production,c)
    write(root/'run/acceptance.json',dict(status='accepted_for_bounded_search',contract_path=str(production),contract_sha256=sha(production),smoke_receipt=str(smoke),smoke_sha256=sha(smoke),search_output=str(root/'run/search'),authority='Lead acceptance under Tommaso daytime continuation authorization'))
print(json.dumps(result,indent=2))
