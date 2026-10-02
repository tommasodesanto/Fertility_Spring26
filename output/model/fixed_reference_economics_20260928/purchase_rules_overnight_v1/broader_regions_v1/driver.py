"""Inject one pinned broad start into the unchanged 80%-financing model runner."""
import argparse
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

SOURCE=Path(__file__).resolve().parent
PACKET=Path(os.environ.get('REGIONS_PACKET_ROOT',SOURCE.parent))
CENSOR=Path(os.environ.get('REGIONS_REVIEWED_CENSOR',PACKET/'restart_controller_v2/controller.py'))

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,obj):Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--slot',type=int,choices=range(16),required=True)
 ap.add_argument('--stage',choices=('init','search','postcheck'),required=True)
 ap.add_argument('--root',type=Path,required=True)
 ap.add_argument('--deadline-epoch',type=float,required=True)
 a=ap.parse_args()
 design=json.loads((SOURCE/'design.json').read_text())
 assert sha(PACKET/'plan.json')==design['original_plan_sha256']
 assert sha(CENSOR)==design['reviewed_censor_source_sha256']
 plan=json.loads((PACKET/'plan.json').read_text())
 assert plan['target_fingerprint']==design['original_target_contract_sha256']
 assert plan['weight_fingerprint']==design['original_weight_fingerprint']
 assert plan['bounds']==design['full_approved_bounds']
 assert plan['free_coordinates']==design['free_coordinates']
 assert a.deadline_epoch<=design['absolute_deadline_epoch']
 assert time.time()<a.deadline_epoch
 start=design['starts'][a.slot]
 chain=start['original_chain_index']
 assert chain==(a.slot if a.slot<8 else a.slot+16)
 assert start['arm']==('hard' if a.slot<8 else 'quarter')
 for k,val in start['parameters'].items():
  lo,hi=plan['bounds'][k];assert lo<=val<=hi
 stage_dir=a.root/a.stage
 assert not stage_dir.exists(), 'Refusing existing broad-region stage'
 a.root.mkdir(parents=True,exist_ok=True)
 sys.path.insert(0,str(PACKET))
 sys.argv=['run_psi.py','--chain',str(chain),'--out',str(stage_dir),
           '--deadline-epoch',str(a.deadline_epoch)]
 if a.stage=='init':sys.argv.append('--initialize-only')
 elif a.stage=='search':sys.argv.append('--fast-objective')
 else:sys.argv.extend(['--verify-only',str(a.root/'search/search_completed.json')])
 import run_psi as original
 assert original.CONFIG['nearby_starts'][chain]['purchase_rule']==start['arm']
 assert original.CONFIG['profiles'][start['weight_profile']]=={}
 original.CONFIG['nearby_starts'][chain]['parameters']=dict(start['parameters'])
 original.CONFIG['nearby_starts'][chain]['profile']=start['weight_profile']
 original.CONFIG['nearby_starts'][chain]['initial_price']=start['initial_price']
 if a.stage=='search':
  spec=importlib.util.spec_from_file_location('reviewed_censor',CENSOR)
  helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
  optimizer_sha=helper.install_censored_optimizer(original,design['maximum_objective_calls_per_chain'])
 else:optimizer_sha=None
 receipt=dict(slot=a.slot,arm=start['arm'],original_chain_index=chain,stage=a.stage,
              design_sha256=sha(SOURCE/'design.json'),runner_sha256=sha(PACKET/'run_psi.py'),
              reviewed_censor_sha256=sha(CENSOR),optimizer_source_sha256=optimizer_sha,
              deadline_epoch=a.deadline_epoch,absolute_deadline_epoch=design['absolute_deadline_epoch'],
              max_objective_calls=design['maximum_objective_calls_per_chain'],
              final_reserve_seconds=design['final_reserve_seconds'],parameters=start['parameters'],
              original_plan_sha256=design['original_plan_sha256'],
              target_contract_sha256=design['original_target_contract_sha256'],
              weight_fingerprint=design['original_weight_fingerprint'],
              experimental_initialization_only=True)
 write(a.root/f'{a.stage}_region_contract.json',receipt)
 original.main()

if __name__=='__main__':main()
