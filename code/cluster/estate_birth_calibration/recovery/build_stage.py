"""Build isolated recovery derivative from pinned failed-parent stage, without solves."""
import gzip, hashlib, io, json, tarfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
PARENT=ROOT/'output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/stage.tar.gz'
OUT=ROOT/'output/model/experiments/birth_count_choice/estate_a_recovery_20261004_v1/deployment'
REMOTE='/scratch/td2248/projects/estate_birth_recovery_20261004_v1'
DRIVER='code/model/experiments/birth_count_choice/cluster_calibrate.py'
PARENT_SHA='d14a39bcbb55067060ca492943f1e0067000c10733a48c3a3aff3b06e61e2afe'
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
    with tarfile.open(PARENT) as tar:
        pin=tar.extractfile('inventory.json').read();assert sha(pin)==PARENT_SHA
        old=json.loads(pin)
        source={n.removeprefix('source/'):tar.extractfile(n).read() for n in tar.getnames() if n.startswith('source/')}
    assert {n:sha(b) for n,b in source.items()}==old['files']
    orig=source[DRIVER].decode();begin=orig.index('def checked_plan(');end=orig.index('def experimental_parameter_rows(',begin)
    replacement=(HERE/'checked_plan.py.txt').read_text();new=orig[:begin]+replacement+orig[end:]
    def change(before,after):
        nonlocal new
        assert new.count(before)==1,before
        new=new.replace(before,after,1)
    change('class BudgetStop(Exception): pass\n\n', '''class BudgetStop(Exception): pass

def ensure_native_launch_budget(now, search_deadline):
    if now >= search_deadline - 800:
        raise BudgetStop('native_700s_minimum_plus_100s_guard')

def authenticated_budget_error(message, remaining, case_exists):
    return (message in ('native GE acceptance failed: uncomputed_bounded_budget',
                        'No time reserve for selected reporting and exact repeat',
                        'No time for exact repeat') and remaining <= 700 and case_exists)

''')
    change('start+21600','start+43200')
    change('wall_seconds=21600,deterministic_seed','wall_seconds=43200,deterministic_seed')
    change("require(0<=args.chain<5,'Matched chain outside 0..4');cap=plan['arms'][args.arm];seed=plan['starts'][args.chain]", "require(plan['recovery_arm']==args.arm, 'Recovery arm/plan mismatch')\n    require(0<=args.chain<len(plan['starts']),'Recovery chain outside plan');cap=plan['arms'][args.arm];seed=plan['starts'][args.chain]")
    change('starts_count=5,','starts_count=len(plan[\'starts\']),')
    change("provisional_seed_source_sha256=plan['provisional_seed_source_sha256'],provisional_seed=plan['provisional_seed']", "provisional_seed_source_sha256=plan['provisional_seed_source_sha256'],provisional_seed=plan['provisional_seed']")
    change("if time.time()>=deadline-RESERVE:raise BudgetStop('final_native_reserve_reached')", "ensure_native_launch_budget(time.time(),deadline-RESERVE)")
    change('else:result=evaluate(label,point,deadline-RESERVE)', '''else:
            try:result=evaluate(label,point,deadline-RESERVE)
            except RuntimeError as exc:
                message=str(exc);remaining=deadline-RESERVE-time.time()
                if authenticated_budget_error(message,remaining,(out/label).is_dir()):
                    write(out/'native_budget_stop.json',dict(reason='authenticated_native_time_exhaustion',
                          original_error=message,remaining_search_seconds=remaining,case_label=label,
                          search_deadline_epoch=deadline-RESERVE,source='pinned_native_engine'))
                    raise BudgetStop('native_evaluation_budget_exhausted') from exc
                raise''')
    assert new.count('def checked_plan(')==1
    source[DRIVER]=new.encode()
    scripts={p.name:p.read_bytes() for p in HERE.iterdir() if p.is_file() and p.suffix in ('.py','.sh') and p.name!='build_stage.py'}
    inv=dict(parent_inventory_sha256=PARENT_SHA,parent_driver_sha256=old['files'][DRIVER],derivative_driver_sha256=sha(source[DRIVER]),derivative_scope='checked_plan; chain count; 12-hour cap; narrow search budget guard/translation',files={n:sha(b) for n,b in sorted(source.items())},entrypoints={n:sha(b) for n,b in sorted(scripts.items())},target_fingerprint=old['target_fingerprint'],weight_fingerprint=old['weight_fingerprint'],selected_source_sha256=old['selected_source_sha256'],remote_root=REMOTE,parent_job_id='19127370',no_auto_retry=True)
    OUT.mkdir(parents=True,exist_ok=True);ib=(json.dumps(inv,sort_keys=True,indent=2)+'\n').encode();(OUT/'inventory.json').write_bytes(ib)
    entries={'source/'+n:b for n,b in source.items()};entries.update(scripts);entries['inventory.json']=ib
    with (OUT/'stage.tar.gz').open('wb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as gz:
        with tarfile.open(fileobj=gz,mode='w') as tar:
            for n,b in sorted(entries.items()):
                info=tarfile.TarInfo(n);info.size=len(b);info.mtime=0;info.mode=0o755 if n.endswith('.sh') else 0o644;tar.addfile(info,io.BytesIO(b))
    receipt=dict(status='prepared_zero_solves',archive=str(OUT/'stage.tar.gz'),archive_sha256=sha((OUT/'stage.tar.gz').read_bytes()),inventory_sha256=sha(ib),source_files=len(source),remote_root=REMOTE)
    (OUT/'stage_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
