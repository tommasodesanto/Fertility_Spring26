"""Fail closed on external gate, source and incumbent identity; zero native solves."""
from pathlib import Path
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def check(gate,manifest,incumbent,own,gate_plan):
 assert gate['status']=='incumbent_replay_passed'
 assert gate['source_manifest_sha256']==own['companion_manifest_sha256']
 assert manifest['incumbent_sha256']==sha(HERE/'incumbent.json')==own['incumbent_sha256']
 assert manifest['target_fingerprint']==incumbent['target_fingerprint']==gate_plan['target_fingerprint']
 assert manifest['weight_fingerprint']==incumbent['weight_fingerprint']==gate_plan['weight_fingerprint']
 assert gate['old_h_P_bounds']==[.1,2.3] and gate['extended_h_P_bounds']==[.1,2.6]
 assert gate['all_other_coordinates_fixed'] is True
 old,ext=gate['old'],gate['extended']
 for arm in (old,ext):
  assert arm['status']=='passed' and arm['lifecycle_solves']>0
  assert arm['repeat']['status']=='exact_full_ge_repeat_passed'
  assert arm['repeat']['target_rows']==14 and arm['repeat']['parameter_rows']==31 and len(arm['repeat']['standard_plot_hashes'])==17
  assert len(arm['target_fit'])==14 and len(arm['parameters'])==31
 assert old['target_fit']==ext['target_fit']==incumbent['postcheck_target_fit']
 assert old['parameters']==incumbent['postcheck_parameters']
 for a,b in zip(old['parameters'],ext['parameters']):
  assert a['parameter']==b['parameter'] and a['estimate']==b['estimate']
  if a['parameter']!='h_P':assert a==b
  else:
   assert b['lower']=='0.1' and b['upper']=='2.6'
   assert {k:v for k,v in a.items() if k not in ('upper','near_bound')}=={k:v for k,v in b.items() if k not in ('upper','near_bound')}
 assert old['price']==ext['price']==incumbent['selected_price']
 assert old['H0_derived']==ext['H0_derived']
 assert abs(old['base_loss']-incumbent['base_loss'])<1e-10
 return {'status':'authenticated_external_gate','native_solves_this_check':0}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--gate-root',type=Path,required=True);ap.add_argument('--gate-job',required=True);args=ap.parse_args()
 assert args.gate_job.isdigit()
 own=json.loads((HERE/'manifest.json').read_text());inc=json.loads((HERE/'incumbent.json').read_text());plan=json.loads((HERE/'plan.json').read_text())
 for rel,digest in own['sha256'].items():assert sha(ROOT/rel)==digest,rel
 remote=args.gate_root;packet=Path('output/model/fixed_reference_economics_20260928/normalized_floor_extension_v2');comp=remote/'source'/packet
 manifest=json.loads((comp/'manifest.json').read_text());assert sha(comp/'manifest.json')==own['companion_manifest_sha256']
 assert sha(comp/'incumbent.json')==sha(HERE/'incumbent.json')
 assert sha(comp/'run.py')==manifest['runner_sha256']
 v2root=remote.parent/'normalized_calibration_v2'
 for rel,digest in manifest['sha256'].items():
  folder=remote if 'normalized_floor_extension_v2/' in rel else v2root
  assert sha(folder/'source'/rel)==digest,rel
 # Exact complete source-pins bytes and every inherited pin must agree with v2.
 assert sha(v2root/'source'/packet.parent/'normalized_calibration_v2/normalized_objective.py')==sha(HERE/'normalized_objective.py')
 vp=v2root/'source'/packet.parent/'normalized_calibration_v2/source_pins.json'
 assert sha(vp)==sha(HERE/'source_pins.json')
 for rel,digest in json.loads(vp.read_text()).items():assert sha(v2root/'source'/rel)==digest,rel
 gate=json.loads((remote/'results/gate0/run/completed.json').read_text())
 terminal=json.loads((remote/'results/gate0/launcher_terminal.json').read_text())
 assert terminal['exit_code']==0 and str(terminal['slurm_job_id']).split('_')[0]==args.gate_job
 print(json.dumps(check(gate,manifest,inc,own,plan)))
if __name__=='__main__':main()
