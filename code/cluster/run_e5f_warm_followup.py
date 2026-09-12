"""After a successful warm terminal, reproduce it and run a short native forecast."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--arm',type=int,required=True);a=parser.parse_args()
    root=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a')
    batch=root/'batches/night_warm_followup_20260912';warm=root/'batches/night_warm_terminal_20260912'
    result=warm/f'arm_{a.arm}';receipt_path=result/'root_receipt.json'
    if not (result/'summary.json').exists() or read(result/'summary.json')['status']!='passed_terminal_root_diagnostic':
        (batch/f'skipped_{a.arm}.json').write_text(json.dumps(dict(reason='Warm terminal did not pass; no dependent forecast launched'))+'\n');return
    receipt=read(receipt_path);assert receipt['converged'];f=receipt['final']
    plan=read(root/'batches/night_surprise_frontier_20260912/plan.json')
    plan.update(output_root=str(batch/'results'),seed_steps=[-.01414,-.0175],
        terminal_driver=str(warm/'warm_terminal_driver.py'),
        terminal_contract_overrides=dict(schema='e5f_candidate_terminal_warm_start_v1',diagnostic_root_start=dict(asset_price=f['prices'][0],pension_period=f['fiscal_values'][0],source_receipt=dict(path=str(receipt_path),sha256=sha(receipt_path)))),
        diagnostic_question='Use a verified same-preference terminal as numerical start, reproduce it, and inspect the first short forecast. No historical fit or policy promotion.')
    plan['file_sha256'].update({str(p):sha(p) for p in [*batch.glob('*.py'),warm/'warm_terminal_driver.py',receipt_path,result/'summary.json']})
    path=batch/f'plan_{a.arm}.json';path.write_text(json.dumps(plan,indent=2)+'\n')
    env=dict(os.environ,PYTHONPATH=f'{batch}:{root}/code/model/tools:{root}/code/model')
    subprocess.run([sys.executable,'-B','-m','unittest','test_e5f_surprise_overnight','test_e5f_successive_surprises','-q'],env=env,check=True,timeout=180)
    subprocess.run([sys.executable,'-B',str(batch/'run_e5f_successive_surprises_overnight.py'),'--plan',str(path),'--arm',str(a.arm)],env=env,check=True,timeout=5460)
if __name__=='__main__':main()
