"""Build a pinned derivative of the verified October 4 recovery stage."""
import gzip, hashlib, io, json, tarfile
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
PARENT_ARCHIVE=ROOT/'output/model/experiments/birth_count_choice/estate_a_recovery_20261004_v1/deployment/stage.tar.gz'
OUT=ROOT/'output/model/experiments/birth_count_choice/estate_a_binary_continuation_20261004_v2/deployment'
DRIVER='code/model/experiments/birth_count_choice/cluster_calibrate.py'
RECOVERY_SHA='974a14e243da6a2ad0572bb9825b47ab349828f9144cbbad69e74b40a4408b22'
REMOTE='/scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2'

def sha(b): return hashlib.sha256(b).hexdigest()

def main():
    with tarfile.open(PARENT_ARCHIVE) as tar:
        ib=tar.extractfile('inventory.json').read(); assert sha(ib)==RECOVERY_SHA,'Recovery inventory drift'
        old=json.loads(ib)
        source={n.removeprefix('source/'):tar.extractfile(n).read() for n in tar.getnames() if n.startswith('source/')}
    assert {n:sha(b) for n,b in source.items()}==old['files'],'Recovery source archive drift'
    driver=source[DRIVER].decode()
    a=driver.index('def checked_plan(');b=driver.index('def experimental_parameter_rows(',a)
    replacement=(HERE/'checked_plan.py.txt').read_text()
    driver=driver[:a]+replacement+driver[b:]
    oldcheck="require(plan['recovery_arm']==args.arm, 'Recovery arm/plan mismatch')"
    assert driver.count(oldcheck)==1
    driver=driver.replace(oldcheck,"require(plan['continuation_arm']==args.arm, 'Continuation arm/plan mismatch')",1)
    # Stop only the optimizer after a recorded full-GE loss crosses the threshold.
    needle="heartbeat('case_completed',objective_calls=calls,completed_full_ge=len(cases),best_loss=best['loss'] if best else None)\n        return row['objective']"
    assert driver.count(needle)==1
    driver=driver.replace(needle,"heartbeat('case_completed',objective_calls=calls,completed_full_ge=len(cases),best_loss=best['loss'] if best else None)\n        if best is not None and best['loss'] < 13.:\n            raise BudgetStop('provisional_loss_below_13_fresh_native_gate_required')\n        return row['objective']",1)
    source[DRIVER]=driver.encode()
    scripts={p.name:p.read_bytes() for p in HERE.iterdir() if p.is_file() and p.suffix in ('.py','.sh','.txt') and p.name!='build_stage.py'}
    inv=dict(parent_inventory_sha256=RECOVERY_SHA,parent_job_id='19141024',
      parent_driver_sha256=old['files'][DRIVER],derivative_driver_sha256=sha(source[DRIVER]),
      derivative_scope='checked plan authentication; 20-chain continuation plan; provisional <13 early search stop only before unchanged native selected-point repeat gate',
      files={n:sha(b) for n,b in sorted(source.items())},entrypoints={n:sha(b) for n,b in sorted(scripts.items())},
      target_fingerprint=old['target_fingerprint'],weight_fingerprint=old['weight_fingerprint'],
      selected_source_sha256=old['selected_source_sha256'],remote_root=REMOTE,no_auto_retry=True,
      objective_max_calls_per_task=500,task_wall_seconds=43200,native_reserve_seconds=1800,task_count=20)
    OUT.mkdir(parents=True,exist_ok=True);ib=(json.dumps(inv,sort_keys=True,indent=2)+'\n').encode()
    (OUT/'inventory.json').write_bytes(ib)
    entries={'source/'+n:b for n,b in source.items()};entries.update(scripts);entries['inventory.json']=ib
    with (OUT/'stage.tar.gz').open('wb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as gz:
      with tarfile.open(fileobj=gz,mode='w') as tar:
        for n,b in sorted(entries.items()):
          ti=tarfile.TarInfo(n);ti.size=len(b);ti.mtime=0;ti.mode=0o755 if n.endswith('.sh') else 0o644;tar.addfile(ti,io.BytesIO(b))
    receipt=dict(status='prepared_zero_solves_no_submission',archive=str(OUT/'stage.tar.gz'),archive_sha256=sha((OUT/'stage.tar.gz').read_bytes()),
      inventory_sha256=sha(ib),source_files=len(source),parent_inventory_sha256=RECOVERY_SHA,remote_root=REMOTE)
    (OUT/'stage_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))

if __name__=='__main__': main()
