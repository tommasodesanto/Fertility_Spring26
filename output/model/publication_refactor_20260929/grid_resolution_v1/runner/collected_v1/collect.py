"""Read-only compact remote collection; excludes arrays, checkpoint and PNG bytes."""
import base64,hashlib,json,shlex,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
REMOTE='/scratch/td2248/projects/grid_resolution_120x9_v1/results'
script=r'''
import base64,hashlib,json
from pathlib import Path
root=Path('/scratch/td2248/projects/grid_resolution_120x9_v1/results')
patterns=['preflight/completed.json','full/failure.json','full/control_160x15/completed.json','full/control_160x15/effective_input_contract.json','full/control_160x15/external_workflow_timing.json','full/control_160x15/frozen_observer_identity.json','full/control_160x15/seed/summary.json','full/control_160x15/phase_b_ge/selected.json','full/control_160x15/phase_b_ge/*/closure.json','full/control_160x15/phase_b_ge/selected_root/target_fit.csv','full/control_160x15/phase_b_ge/selected_root/parameters.csv','full/control_160x15/phase_b_ge/selected_root/gates.json','full/control_160x15/phase_b_ge/selected_repeat_final/target_fit.csv','full/control_160x15/phase_b_ge/selected_repeat_final/parameters.csv','full/control_160x15/phase_b_ge/selected_repeat_final/gates.json','full/proposal_120x9/effective_input_contract.json','full/proposal_120x9/failure.json','full/proposal_120x9/frozen_observer_identity.json','full/proposal_120x9/latest.json','full.log']
paths=set(p for pattern in patterns for p in root.glob(pattern) if p.is_file())
files={str(p.relative_to(root)):dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size,base64=base64.b64encode(p.read_bytes()).decode()) for p in sorted(paths)}
plots={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for final in ['selected_root','selected_repeat_final'] for p in sorted((root/'full/control_160x15/phase_b_ge'/final/'standard_diagnostics').glob('*.png'))}
print(json.dumps(dict(files=files,control_plot_hashes=plots)))
'''
result=subprocess.run(['ssh','torch','/share/apps/anaconda3/2025.06/bin/python -c '+shlex.quote(script)],capture_output=True,text=True,check=True)
packet=json.loads(result.stdout);receipts={}
for rel,entry in packet['files'].items():
    blob=base64.b64decode(entry.pop('base64'));assert hashlib.sha256(blob).hexdigest()==entry['sha256']
    target=HERE/rel;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(blob);receipts[rel]=entry
(HERE/'remote_hash_receipt.json').write_text(json.dumps(dict(remote_root=REMOTE,files=receipts,control_plot_hashes=packet['control_plot_hashes'],remote_writes=False,arrays_downloaded=False),indent=2)+'\n')
print(json.dumps(dict(collected_files=len(receipts),collected_bytes=sum(e['bytes'] for e in receipts.values()),control_plot_hashes=len(packet['control_plot_hashes']))))
