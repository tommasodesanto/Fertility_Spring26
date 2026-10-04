"""Fail-closed plan check for the 20 one-birth continuation starts."""
import argparse,hashlib,json,math,sys
from pathlib import Path
from prepare_starts import audit
RECOVERY='974a14e243da6a2ad0572bb9825b47ab349828f9144cbbad69e74b40a4408b22'
PARENT='d14a39bcbb55067060ca492943f1e0067000c10733a48c3a3aff3b06e61e2afe'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def verify(path,recovery,parent,receipt_path,audit_sources=True):
    p=json.loads(Path(path).read_text());assert len(p['starts'])==20 and len(p['start_provenance'])==20
    assert p['recovery_inventory_sha256']==RECOVERY and p['parent_inventory_sha256']==PARENT
    assert p['continuation_arm']=='binary' and p['arms']=={'binary':1}
    for point in p['starts']:
        assert set(point)==set(p['bounds'])
        assert all(math.isfinite(v) and p['bounds'][k][0]<=v<=p['bounds'][k][1] for k,v in point.items())
    cats=[x['kind'] for x in p['start_provenance']]
    assert cats.count('verified_recovery_selected')==5 and cats.count('parent_saved_feasible_checkpoint')==5 and cats.count('bounded_perturbation')==10
    assert min(p['pairwise_normalized_distances'])>0.005
    receipt=json.loads(Path(receipt_path).read_text())
    assert receipt['starts_sha256']==sha(path) and receipt['count']==20
    if audit_sources: assert audit(p,recovery,parent)['status']=='all_20_starts_reauthenticated'
    return dict(status='passed',sha256=sha(path),starts=20,minimum_normalized_pairwise_distance=min(p['pairwise_normalized_distances']),source_check='passed' if audit_sources else 'performed_during_remote_staging')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('plan',type=Path);ap.add_argument('--recovery',type=Path,default=Path('/scratch/td2248/projects/estate_birth_recovery_20261004_v1'));ap.add_argument('--parent',type=Path,default=Path('/scratch/td2248/projects/estate_birth_calibration_20261003_v3'));ap.add_argument('--receipt',type=Path);ap.add_argument('--structure-only',action='store_true')
    a=ap.parse_args();receipt=a.receipt or a.plan.with_name('starts_receipt.json');print(json.dumps(verify(a.plan,a.recovery,a.parent,receipt,not a.structure_only),indent=2))
