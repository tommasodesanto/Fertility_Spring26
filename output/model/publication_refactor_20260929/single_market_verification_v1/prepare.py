"""Immutable single-arm materialization; no lifecycle calculation or remote operation."""
from pathlib import Path
import ast, hashlib, json, shutil, subprocess, sys, tarfile
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
MODEL=ROOT/'code/model/experiments/stationary_single_market'
PRIOR=ROOT/'output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed'
STAGE=HERE/'source'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if STAGE.exists():raise SystemExit('Refusing source-stage overwrite')
STAGE.mkdir()
pins={}
for name in ['driver.py','phase_a.py','phase_b_ge.py','single_price.py']:
    source=PRIOR/name;pins[str(source.relative_to(ROOT))]=sha(source);shutil.copyfile(source,STAGE/name)
package=STAGE/'source/small_credit_lab';package.mkdir(parents=True)
# Runtime files only; preserve reviewed module organization including ..contract.
for name in ['__init__.py','inputs.py','credit.py','contract.py']:
    source=MODEL/name;pins[str(source.relative_to(ROOT))]=sha(source);shutil.copyfile(source,package/name)
(package/'engine').mkdir()
for source in sorted((MODEL/'engine').glob('*.py')):
    pins[str(source.relative_to(ROOT))]=sha(source);shutil.copyfile(source,package/'engine'/source.name)
for name in ['joint_nested.py','two_shock_choice.py','fertility_nested.py']:
    assert not (package/'engine'/name).exists(),name
for p in STAGE.rglob('*.py'):ast.parse(p.read_text())
shutil.copyfile(HERE/'run_arm.py',STAGE/'run_arm.py')
shutil.copyfile(HERE/'compare.py',STAGE/'compare.py')
(HERE/'preparation.json').write_text(json.dumps(dict(status='prepared_zero_solves',model_commit='1e24652d',source_sha256=pins,deleted_engine_modules=['joint_nested.py','two_shock_choice.py','fertility_nested.py'],bundle_sha256='427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7',global_seconds=1200,case_seconds=300,lifecycle_cap=6,root_reporting_repeat_reserve_seconds=700,unsecured_credit_limit=.14,baseline_job=18879780,baseline='/scratch/td2248/projects/grid_resolution_120x9_v1/results/full/control_160x15',node='cs713',economics='unchanged D14 diagnostic relative to matched grid control; fixed preferences and psi; actual birth renewal and absolute housing supply'),indent=2)+'\n')
# Preserve driver modules byte-for-byte; budget/controller adaptation is in run_arm.
assert all(sha(STAGE/name)==sha(PRIOR/name) for name in ['driver.py','phase_a.py','phase_b_ge.py','single_price.py'])
assert all(sha(ROOT/rel)==digest for rel,digest in pins.items())
for rel,digest in pins.items():
    if rel.startswith('code/model/experiments/stationary_single_market/'):
        blob=subprocess.check_output(['git','show','1e24652d:'+rel])
        assert hashlib.sha256(blob).hexdigest()==digest,rel
with (STAGE/'source.sha256').open('w') as f:
    for p in sorted(STAGE.rglob('*')):
        if p.is_file() and p.name!='source.sha256':f.write(f'{sha(p)}  {p.relative_to(STAGE)}\n')
with tarfile.open(HERE/'source.tar.gz','w:gz') as archive:archive.add(STAGE,arcname='source')
(HERE/'source_archive.sha256').write_text(sha(HERE/'source.tar.gz')+'  source.tar.gz\n')
print('prepared',len(pins),'pinned source files; zero solves')
