from pathlib import Path
import hashlib,json,os,subprocess,sys
cases=json.loads(Path('contracts/cases.json').read_text())
index=int(os.environ['SLURM_ARRAY_TASK_ID']); case=cases[index]
assert index==case['index'] and 0<=index<19
path=Path('contracts')/case['contract']
assert hashlib.sha256(path.read_bytes()).hexdigest()==case['contract_sha256']
subprocess.run([sys.executable,'code/model/tools/run_e5f_initial_revision_probe.py','--contract',str(path),'--contract-sha256',case['contract_sha256'],'--case','new_balanced','--output','output/panel/'+case['case_id']],check=True)
