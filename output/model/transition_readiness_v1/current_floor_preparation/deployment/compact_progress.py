import json,time
from pathlib import Path
root=Path('/scratch/td2248/projects/transition_readiness_v1/current_floor/results/smoke_v4')
def read(p):
 try:return json.loads(p.read_text())
 except (OSError,ValueError):return None
summary=dict(captured_epoch=time.time(),label='smoke_v4',heartbeat=read(root/'heartbeat.json'),launcher_terminal=read(root/'launcher_terminal.json'),failure=read(root/'failure.json'),reference_reconstruction=read(root/'native_reference/reference_reconstruction.json'),seed_receipt=read(root/'seed_receipt.json'),maps=[])
if summary['reference_reconstruction']:
 summary['reference_reconstruction']={k:v for k,v in summary['reference_reconstruction'].items() if k in ('status','policy_calls','lifecycle_calls','checkpoint_sha256')}
for p in sorted(root.glob('seed/map_*/mapping.json'))+sorted(root.glob('six_date_smoke/**/mapping.json')):
 r=read(p)
 if not r:continue
 def maximum(key):
  v=r.get(key,[])
  return max((abs(float(x)) for x in v),default=None)
 row=dict(path=str(p.relative_to(root)),mtime=p.stat().st_mtime,policy_calls=r.get('policy_calls'),total_native_calls=r.get('total_native_calls'),accounting_valid=r.get('accounting_valid'),gates=r.get('gates'),market_maximum_residual=maximum('market_residual'),fiscal_maximum_residual=maximum('fiscal_residual'))
 if row['market_maximum_residual'] is not None and row['fiscal_maximum_residual'] is not None:
  row['maximum_gate_ratio']=max(row['market_maximum_residual']/2e-4,row['fiscal_maximum_residual']/2e-5)
 summary['maps'].append(row)
summary['latest_completed_map']=max(summary['maps'],key=lambda r:r['mtime']) if summary['maps'] else None
summary['best_completed_map']=min((r for r in summary['maps'] if 'maximum_gate_ratio' in r),key=lambda r:r['maximum_gate_ratio'],default=None)
summary['mapping_count']=len(summary['maps'])
endpoint=read(root/'six_date_smoke/endpoint/root.json')
summary['endpoint_root']={k:endpoint.get(k) for k in ('converged','status','evaluations','elapsed_seconds','history')} if endpoint else None
summary['effective_status']='FAILED' if summary['failure'] else ('terminal' if summary['launcher_terminal'] else 'running')
summary['reference_partials']=[]
for p in sorted(root.glob('native_reference/repeat_*/native_solve_unverified.json')):
 r=read(p)
 if r:summary['reference_partials'].append(dict(path=str(p.relative_to(root)),mtime=p.stat().st_mtime,status=r.get('status'),policy_calls=r.get('policy_calls'),total_native_calls=r.get('total_native_calls')))
print(json.dumps(summary,indent=2))
