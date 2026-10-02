"""Read compact actualfit_v1 evidence; never import the model or write its run."""
import json,time
from pathlib import Path
root=Path('/scratch/td2248/projects/transition_readiness_v1/normalized_resumed_fit_v1/results/normalized_resumed_fit_v1')
def read(p):
 try:return json.loads(p.read_text())
 except (OSError,ValueError):return None
def compact(r):
 if isinstance(r,dict):
  return {k:compact(v) for k,v in r.items() if k not in ('identity','reference_identity','source_files','target_contract','preparation','source_pins','plot_sha256','reports','native_reply','matrix','diagnostic_packets')}
 if isinstance(r,list):return [compact(v) for v in r] if len(r)<=32 else dict(omitted_list_length=len(r))
 return r
r=dict(captured_epoch=time.time(),run=str(root))
for name in ('launcher_start.json','launcher_terminal.json','heartbeat.json','failure.json','native_input_reuse.json','latest_completed.json','best_so_far.json','latest_exploratory_state.json','complete.json'):
 x=read(root/name)
 if x is not None:r[name]=compact(x)
r['candidate_receipts']=[]
for p in sorted(root.glob('candidate_*/complete.json')):
 x=read(p)
 if x is not None:r['candidate_receipts'].append(dict(path=str(p),summary=compact(x)))
r['checkpoint_receipts']=[]
for p in sorted(root.glob('candidate_*/state_2023_checkpoint/checkpoint_receipt.json')):
 x=read(p)
 if x is not None:r['checkpoint_receipts'].append(dict(path=str(p),summary=compact(x)))
r['status']='FAILED' if 'failure.json' in r else ('COMPLETE' if 'complete.json' in r else ('RUNNING' if 'launcher_start.json' in r else 'NOT_STARTED'))
print(json.dumps(r,indent=2))
