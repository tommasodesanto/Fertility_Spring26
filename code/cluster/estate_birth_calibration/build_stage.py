"""Bundle authenticated parent dependencies plus matched estate-A source inputs."""
import gzip,hashlib,io,json,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
PACKETS=Path('output/model/fixed_reference_economics_20260928')
OUT=ROOT/'output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3'
PARENT=ROOT/PACKETS/'soft_timing_calibration_20261002_v1/deployment/attempt2/stage.tar.gz'
PARENT_SHA='f7a8fec4ff370fd3690c0d0068ca595b75a17dd8aaac3bd47f6009ef73ecd68b'
STARTS='output/model/experiments/birth_count_choice/estate_a_calibration_v1/start_plan.json'
ANCHOR=str(PACKETS/'soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json')
REMOTE='/scratch/td2248/projects/estate_birth_calibration_20261003_v3'
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
    assert sha(PARENT.read_bytes())==PARENT_SHA,'Parent archive drift'
    with tarfile.open(PARENT) as archive:
        prior=json.load(archive.extractfile('inventory.json'))
        source={n.removeprefix('source/'):archive.extractfile(n).read() for n in archive.getnames() if n.startswith('source/')}
    assert {k:sha(v) for k,v in source.items()}==prior['files'],'Parent inventory drift'
    # Recover byte-exact parent dependencies; add current experiment sources and portable caller inputs.
    for directory in ('code/model/experiments/birth_count_choice','code/model/production/reference_inputs'):
        for path in (ROOT/directory).rglob('*'):
            if path.is_file() and '__pycache__' not in path.parts and path.suffix not in ('.pyc','.nbc','.nbi'):
                source[str(path.relative_to(ROOT))]=path.read_bytes()
    for rel in (STARTS,ANCHOR,str(PACKETS/'alternative_wealth_cluster_20261003_v1/pause_20261003.json')):source[rel]=(ROOT/rel).read_bytes()
    overlay='output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/local_runtime'
    for rel in (overlay+'/bootstrap.py',overlay+'/frozen_sources/e5f_exact_policy_cache.py',overlay+'/frozen_sources/test_e5f_exact_policy_cache.py'):
        source[rel]=(ROOT/rel).read_bytes()
    expected={'e5f_exact_policy_cache.py':'d51bbd13026288db0194b44b420bb49ff6970588f40a906242978d49038a5d6f','test_e5f_exact_policy_cache.py':'50784af81c5e6fd64d345d504c8c7fa81209b05a6a655789a224022485208ae5'}
    for name,digest in expected.items():assert sha(source[overlay+'/frozen_sources/'+name])==digest,'Frozen overlay drift: '+name
    table=json.loads(source[STARTS]);assert table['source_checkpoint_sha256']==sha(source[ANCHOR])
    assert sha(source[table['provisional_seed_source']])==table['provisional_seed_source_sha256']
    assert len(table['starts'])==5 and table['bounds']['beta_annual']==[.93,.99]
    entrypoints={p.name:sha(p.read_bytes()) for p in HERE.iterdir() if p.is_file() and p.suffix in ('.py','.sh') and p.name!='build_stage.py'}
    inventory=dict(files={k:sha(v) for k,v in sorted(source.items())},entrypoints=entrypoints,
      parent_archive_sha256=PARENT_SHA,target_fingerprint=table['target_fingerprint'],weight_fingerprint=table['weight_fingerprint'],
      selected_source_sha256=sha(source[ANCHOR]),start_plan_sha256=sha(source[STARTS]),remote_root=REMOTE,
      source_prefix='source/',no_cache_or_results=True,arms=table['arms'],starts_count=5,chains=10)
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'inventory.json').write_text(json.dumps(inventory,indent=2,sort_keys=True)+'\n')
    entries={'source/'+k:v for k,v in source.items()};entries['inventory.json']=(OUT/'inventory.json').read_bytes()
    entries.update({n:(HERE/n).read_bytes() for n in entrypoints})
    with (OUT/'stage.tar.gz').open('wb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as gz:
        with tarfile.open(fileobj=gz,mode='w') as archive:
            for name,blob in sorted(entries.items()):
                info=tarfile.TarInfo(name);info.size=len(blob);info.mtime=0;info.mode=0o755 if name.endswith('.sh') else 0o644
                archive.addfile(info,io.BytesIO(blob))
    receipt=dict(status='prepared_no_submission',archive=str(OUT/'stage.tar.gz'),sha256=sha((OUT/'stage.tar.gz').read_bytes()),
       source_files=len(source),start_plan_sha256=inventory['start_plan_sha256'],remote_root=REMOTE,
       target_fingerprint=table['target_fingerprint'],weight_fingerprint=table['weight_fingerprint'])
    (OUT/'stage_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
