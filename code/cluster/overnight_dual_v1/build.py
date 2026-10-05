"""Build two isolated derivative stages from authenticated prior archives."""
import hashlib,json,tarfile,copy,re
from pathlib import Path
H=Path(__file__).resolve().parent;R=H.parents[2]
OUT=R/'output/model/overnight_dual_20261004_v1'
PA=R/'output/model/experiments/birth_count_choice/estate_a_binary_continuation_20261004_v2/deployment'
PB=R/'output/model/fixed_reference_economics_20260928/soft_timing_continuation_20261003_v1/deployment'
RA='/scratch/td2248/projects/estate_birth_overnight_20261004_v1'
RB='/scratch/td2248/projects/soft_timing_overnight_20261004_v1'
DA='code/model/experiments/birth_count_choice/cluster_calibrate.py'
DB='code/cluster/soft_timing_calibration/continuation/calibrate.py'
BP='output/model/fixed_reference_economics_20260928/soft_timing_continuation_20261003_v1/start_plan.json'
def sha(b):return hashlib.sha256(b).hexdigest()
def read_source(archive,rel):
 with tarfile.open(archive) as t:return t.extractfile('source/'+rel).read()
def write(p,b):p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
def save_json(p,x):write(p,(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n').encode())
def replace_once(s,a,b):
 assert s.count(a)==1,(a,s.count(a));return s.replace(a,b,1)
def exception_block(call,indent):
 i=indent
 return (i+'try:result = '+call+'\n'+i+'except RuntimeError as exc:\n'
  +i+'    msg=str(exc)\n'
  +i+'    if msg == "native GE acceptance failed: uncomputed_price_unbracketed":\n'
  +i+'        result=dict(status="inadmissible_numerical",reason=msg,error_type=type(exc).__name__,rejection_kind="price_unbracketed_diagnostic_caps")\n'
  +i+'    elif msg == "native GE acceptance failed: uncomputed_bounded_budget":\n'
  +i+'        result=dict(status="inadmissible_numerical",reason=msg,error_type=type(exc).__name__,rejection_kind="native_bounded_budget")\n'
  +i+'    elif type(exc).__name__ == "InheritedDistributionInfeasible" and getattr(exc,"classification",None)=="inherited_distribution_infeasible" and hasattr(exc,"audit"):\n'
  +i+'        result=dict(status="inadmissible_numerical",reason=msg,error_type=type(exc).__name__,rejection_kind="typed_inherited_distribution_gate")\n'
  +i+'    elif type(exc).__name__ == "InfeasibleThetaError" and hasattr(exc,"stage"):\n'
  +i+'        result=dict(status="inadmissible_numerical",reason=msg,error_type=type(exc).__name__,rejection_kind="typed_dead_node_gate")\n'
  +i+'    else:raise')
def prepare_a():
 d=OUT/'a/deployment';d.mkdir(parents=True,exist_ok=True)
 inv=json.loads((PA/'inventory.json').read_text());source=read_source(PA/'stage.tar.gz',DA).decode()
 a=source.index('def checked_plan(');b=source.index('def experimental_parameter_rows(',a)
 source=source[:a]+(H/'a_checked_plan.py.txt').read_text()+source[b:]
 source=replace_once(source,'        ensure_native_launch_budget(time.time(),deadline-RESERVE)',
                     '        if not args.mock_smoke:ensure_native_launch_budget(time.time(),deadline-RESERVE)')
 source=replace_once(source,"        if best is not None and best['loss'] < 13.:",
                     "        if not args.mock_smoke and best is not None and best['loss'] < 13.:")
 old="                raise\n        row=dict(label=label,parameters=point,**result)"
 add=("                if message == 'native GE acceptance failed: uncomputed_price_unbracketed':\n"
      "                    result=dict(status='inadmissible_numerical',reason=message,error_type=type(exc).__name__,rejection_kind='price_unbracketed_diagnostic_caps')\n"
      "                elif type(exc).__name__ == 'InheritedDistributionInfeasible' and getattr(exc,'classification',None)=='inherited_distribution_infeasible' and hasattr(exc,'audit'):\n"
      "                    result=dict(status='inadmissible_numerical',reason=message,error_type=type(exc).__name__,rejection_kind='typed_inherited_distribution_gate')\n"
      "                elif type(exc).__name__ == 'InfeasibleThetaError' and hasattr(exc,'stage'):\n"
      "                    result=dict(status='inadmissible_numerical',reason=message,error_type=type(exc).__name__,rejection_kind='typed_dead_node_gate')\n"
      "                else:raise\n        row=dict(label=label,parameters=point,**result)")
 source=replace_once(source,old,add)
 write(d/'overlays'/DA,source.encode());inv['files'][DA]=sha(source.encode());inv['derivative_driver_sha256']=inv['files'][DA]
 plan=json.loads((PA/'starts.json').read_text());plan['overnight_parent_starts_sha256']=sha((PA/'starts.json').read_bytes())
 plan['starts']=plan['starts'][:10];plan['start_provenance']=plan['start_provenance'][:10]
 plan['overnight_canceled_sources']=[]
 for n in (1,2):
  name=f'canceled_chain{n}_best_so_far.json';blob=(PA/name).read_bytes();write(d/'control'/name,blob)
  best=json.loads(blob)['best'];assert best['status']=='passed' and best['target_fingerprint']==inv['target_fingerprint'] and best['weight_fingerprint']==inv['weight_fingerprint']
  plan['starts'].append(best['parameters']);plan['start_provenance'].append(dict(index=10+len(plan['overnight_canceled_sources']),kind='canceled_saved_feasible_provisional',source=name))
  plan['overnight_canceled_sources'].append(dict(filename=name,sha256=sha(blob),loss=best['loss']))
 assert len(plan['starts'])==12 and len({json.dumps(x,sort_keys=True) for x in plan['starts']})==12
 save_json(d/'control/starts.json',plan);inv['remote_root']=RA;inv['task_count']=12;inv['task_wall_seconds']=43200
 save_json(d/'inventory.json',inv)
 l=(H.parent/'estate_birth_calibration/binary_continuation_v2/launch_torch.sh').read_text()
 l=l.replace('/scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2',RA)
 l=l.replace('preflight|smoke|production','preflight|mock|production').replace('^(0|[1-9]|1[0-9])$','^(0|[1-9]|1[01])$')
 l=l.replace('"$python" "$remote/verify_starts.py" "$remote/control/starts.json"\n','')
 a=l.index('if [[ "$mode" == production ]]');b=l.index('"$python" - "$inputs"',a);l=l[:a]+l[b:]
 l=l.replace('binds=(--bind "$frozen:$repo:ro" --bind "$remote:/work/deployment:ro"','binds=(--bind "$frozen:$repo:ro" --bind "$remote:/work/deployment:ro" --bind "/scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2:/work/parent_v2:ro"')
 l=l.replace('else wall_seconds=5400; fi','else wall_seconds=300; fi')
 l=l.replace('[[ "$mode" == smoke ]] && options+=(--smoke)','[[ "$mode" == mock ]] && options+=(--mock-smoke)')
 l=l.replace('binds+=(--bind "$out:/work/results:rw")','binds+=(--bind "$out:/work/results:rw")\napptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" /work/deployment/verify_stage.py --container')
 write(d/'launch_torch.sh',l.encode());return inv
def prepare_b():
 d=OUT/'b/deployment';d.mkdir(parents=True,exist_ok=True)
 inv=json.loads((PB/'inventory.json').read_text());source=read_source(PB/'stage.tar.gz',DB).decode()
 source=replace_once(source,'start + 21600','start + 43200')
 source=replace_once(source,'        if time.time() >= deadline - RESERVE:',
                     '        if time.time() >= deadline - RESERVE - 800:')
 source=replace_once(source,'        result = evaluate(label, point, deadline - RESERVE)',exception_block('evaluate(label, point, deadline - RESERVE)','        '))
 write(d/'overlays'/DB,source.encode());inv['files'][DB]=sha(source.encode())
 plan=json.loads((PB.parent/'start_plan.json').read_text())
 for n in (6,8):
  rel=f'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_{n}/run/completed.json'
  blob=(R/rel).read_bytes();receipt=json.loads(blob);assert receipt['status']=='selected_numerically_verified' and receipt['arm']=='alternative'
  write(d/'overlays'/rel,blob);inv['files'][rel]=sha(blob)
  plan['starts'].append(receipt['selected']['parameters']);plan['start_provenance'].append(dict(group='verified_alternative',native_loss=receipt['native_loss'],source_chain=n,source_path=rel,source_sha256=sha(blob),source_status=receipt['status']))
 assert len(plan['starts'])==12 and len({json.dumps(x,sort_keys=True) for x in plan['starts']})==12
 plan['chain_count']=12;plan['wall_seconds_per_chain']=43200
 payload=(json.dumps(plan,indent=2,sort_keys=True,allow_nan=False)+'\n').encode();write(d/'overlays'/BP,payload);inv['files'][BP]=sha(payload);inv['start_plan_sha256']=sha(payload)
 inv['remote_root']=RB;save_json(d/'inventory.json',inv)
 l=(H.parent/'soft_timing_calibration/continuation/launch_torch.sh').read_text()
 l=l.replace('/scratch/td2248/projects/soft_timing_continuation_20261003_v1',RB)
 l=l.replace('#SBATCH --time=06:00:00','#SBATCH --time=12:00:00').replace('^(0|[1-9])$','^(0|[1-9]|1[01])$').replace('^[0-9]$','^(0|[1-9]|1[01])$')
 l=l.replace('wall_seconds=21600; else wall_seconds=5400','wall_seconds=43200; else wall_seconds=300')
 l=l.replace('binds+=(--bind "$out:/work/results:rw")','binds+=(--bind "$out:/work/results:rw")\napptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" /work/deployment/verify_stage.py --container')
 write(d/'launch_torch.sh',l.encode());return inv
def finish(arm,inv):
 d=OUT/arm/'deployment';remote=RA if arm=='a' else RB
 verify=(f'''import argparse,hashlib,json\nfrom pathlib import Path\nROOT=Path({remote!r});REPO=Path({str(R)!r})\ndef sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()\nif __name__=="__main__":\n p=argparse.ArgumentParser();p.add_argument("--host",action="store_true");p.add_argument("--container",action="store_true");a=p.parse_args();assert a.host!=a.container\n root=Path("/work/deployment") if a.container else ROOT;m=json.loads((root/"manifest.json").read_text());inv=json.loads((root/"inventory.json").read_text());source=REPO if a.container else root/"source"\n assert sha(root/"inventory.json")==m["inventory_sha256"]\n for rel,digest in inv["files"].items():assert sha(source/rel)==digest,rel\n for rel,digest in m["entrypoints"].items():assert sha(root/rel)==digest,rel\n print(json.dumps(dict(status="passed_zero_solves",arm={arm!r},files=len(inv["files"]),target=inv["target_fingerprint"],weight=inv["weight_fingerprint"])))\n''').encode()
 write(d/'verify_stage.py',verify)
 manifest=dict(arm=arm,parent_inventory_sha256=sha((PA if arm=='a' else PB).joinpath('inventory.json').read_bytes()),inventory_sha256=sha((d/'inventory.json').read_bytes()),entrypoints={n:sha((d/n).read_bytes()) for n in ('launch_torch.sh','verify_stage.py')},target_fingerprint=inv['target_fingerprint'],weight_fingerprint=inv['weight_fingerprint'],source_changes=list(str(p.relative_to(d/'overlays')) for p in (d/'overlays').rglob('*') if p.is_file()))
 save_json(d/'manifest.json',manifest)
 print(arm,remote,'manifest',sha((d/'manifest.json').read_bytes()),'source_changes',manifest['source_changes'])
if __name__=='__main__':
 OUT.mkdir(parents=True,exist_ok=True)
 a=prepare_a();b=prepare_b();finish('a',a);finish('b',b)
