import hashlib,json,time
from pathlib import Path
from experiments.ces_normalized_shares import adapter
root=Path('/work/deployment')
packet=json.loads((root/'followup_tools/historical_start_point.json').read_text())
out=Path('/work/results/run');out.mkdir(parents=True,exist_ok=False)
end=time.time()+1500
(out/'start.json').write_text(json.dumps(dict(packet,deadline_epoch=end,maximum_lifecycle_solves=32,diagnostic_only=True,no_search=True),indent=2)+'\n')
P,grid=adapter.load_inputs(packet['parameters'])
evaluate=adapter.make_evaluator(out,'historical_start',P,grid,end,target_fingerprint=packet['target_fingerprint'],weight_fingerprint=packet['weight_fingerprint'])
result=evaluate('historical_start',packet['parameters'],end)
(out/'completed.json').write_text(json.dumps(dict(result,diagnostic_only=True,no_search=True,source=packet['source'],source_parameters_sha256=packet['parameter_csv_sha256']),indent=2)+'\n')
print(json.dumps({k:result.get(k) for k in ('status','reason','loss','price','lifecycle_solves','report')}))
