"""Independent shock branch; wait for exact-loop smoke before corresponding stage."""
import json, sys, time
from pathlib import Path
import run_candidate_path as path

FIRST=path.HERE/'results/delta_m005'

def wait_for_gate(stage,key,seconds):
 deadline=time.monotonic()+seconds
 while time.monotonic()<deadline:
  source=FIRST/stage/'summary.json'
  if source.exists():
   summary=json.loads(source.read_text())
   if not summary.get(key):raise RuntimeError(f'Initial {stage} smoke not accepted; no dependent stage launched')
   print(f'Accepted {stage} smoke gate',flush=True);return
  if (FIRST/'failure.json').exists():raise RuntimeError('Initial smoke failed; dependent branch withheld')
  time.sleep(30)
 raise TimeoutError(f'Initial {stage} gate did not appear within bounded wait')

if __name__=='__main__':
 wait_for_gate('terminal','endpoint_numerically_verified',1900)
 original=path.run_stage
 def gated(folder,name,contract,driver):
  if name=='history_6':wait_for_gate('history_6','finite_horizon_market_fiscal_converged',1900)
  return original(folder,name,contract,driver)
 path.run_stage=gated
 path.main()
