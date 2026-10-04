"""Write immutable entrypoint pins and verified control bytes before Torch staging."""
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
OUT=ROOT/'output/model/experiments/birth_count_choice/estate_a_global_search_20261004_v1/deployment'
V2=ROOT/'output/model/experiments/birth_count_choice/estate_a_binary_continuation_20261004_v2/deployment'
RECOVERY=ROOT/'output/model/experiments/birth_count_choice/estate_a_recovery_20261004_v1/collection/binary/provenance/completed.json'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    assert sha(V2/'inventory.json')=='f27c6c71b967d5850c00107f9f6f765950631739aec21b8668928efae6a40cbe'
    inv=json.loads((V2/'inventory.json').read_text())
    assert inv['target_fingerprint']=='c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70'
    assert inv['weight_fingerprint']=='f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4'
    incumbent=json.loads(RECOVERY.read_text())
    assert incumbent['status']=='selected_numerically_verified' and incumbent['native_loss']==21.275413361071312
    assert incumbent['target_fingerprint']==inv['target_fingerprint'] and incumbent['weight_fingerprint']==inv['weight_fingerprint']
    OUT.mkdir(parents=True,exist_ok=True)
    for name in ('stage_manifest.json','incumbent.json','parent_starts.sha256'):
        assert not (OUT/name).exists(),name+' already staged'
    (OUT/'incumbent.json').write_bytes(RECOVERY.read_bytes())
    (OUT/'parent_starts.sha256').write_text(sha(V2/'starts.json')+'\n')
    names=['prepare_plan.py','explore.py','verify_stage.py','launch_torch.sh','submit_once.py']
    m=dict(stage='estate_a_global_search_20261004_v1',parent_inventory_sha256=sha(V2/'inventory.json'),
           parent_starts_sha256=sha(V2/'starts.json'),incumbent_sha256=sha(RECOVERY),
           entrypoints={n:sha(HERE/n) for n in names},target_fingerprint=inv['target_fingerprint'],
           weight_fingerprint=inv['weight_fingerprint'],no_auto_extension=True)
    (OUT/'stage_manifest.json').write_text(json.dumps(m,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(status='prepared_zero_solves',manifest_sha256=sha(OUT/'stage_manifest.json'),
                          parent_inventory_sha256=m['parent_inventory_sha256'],parent_starts_sha256=m['parent_starts_sha256'])))
if __name__=='__main__':main()
