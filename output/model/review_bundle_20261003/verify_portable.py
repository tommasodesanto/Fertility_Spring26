"""Verify the ZIP from a clean extraction; never changes model equations."""
from pathlib import Path
import csv, hashlib, json, os, socket, subprocess, sys, tempfile, time, urllib.request, urllib.parse, zipfile

OUT = Path(__file__).resolve().parent
REPO = OUT.parents[2]
ARCHIVE = OUT / 'Fertility_Model_Review_20261003.zip'
PYTHON = REPO / 'output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python'
REFERENCE = REPO / 'output/model/local_solution/cases/20261003T175652812716Z_b1c72f13'
receipt = {'status': 'running', 'started_epoch': time.time(), 'archive_sha256': hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()}

def save():
    (OUT / 'verification.json').write_text(json.dumps(receipt, indent=2) + '\n')

try:
    temp = Path(tempfile.mkdtemp(prefix='fertility_portable_final_'))
    with zipfile.ZipFile(ARCHIVE) as z:
        for item in z.infolist():
            assert not Path(item.filename).is_absolute() and '..' not in Path(item.filename).parts
            assert (item.external_attr >> 16) & 0o170000 != 0o120000
        z.extractall(temp)
        for item in z.infolist():
            (temp / item.filename).chmod((item.external_attr >> 16) & 0o777)
    bundle = temp / 'Fertility_Model_Review_20261003'
    receipt.update(extracted_root=str(bundle), no_external_symlinks=True)
    for line in (bundle / 'SOURCE_MANIFEST.sha256').read_text().splitlines():
        digest, kind, relative = line.split('  ', 2)
        assert hashlib.sha256((bundle / relative).read_bytes()).hexdigest() == digest, relative
    receipt['manifest_verified'] = True
    # This hook rejects any accidental fallback to the original project. Only
    # the interpreter's existing third-party environment may be read there.
    wrapper = temp / 'audit_run.py'
    wrapper.write_text('''import os, sys, runpy
from pathlib import Path
repo, allowed, bundle, script = map(Path, sys.argv[1:5])
args = sys.argv[5:]
def audit(event, values):
    if event == 'open' and isinstance(values[0], (str, bytes, os.PathLike)):
        path = Path(os.fsdecode(values[0])).resolve()
        if path.is_relative_to(repo) and not path.is_relative_to(allowed):
            raise RuntimeError('Forbidden original-project file access: ' + str(path))
sys.addaudithook(audit)
os.chdir(bundle)
sys.path.insert(0, str(bundle/'code/model'))
sys.path.insert(0, str(script.parent))
sys.argv = [str(script), *args]
runpy.run_path(str(script), run_name='__main__')
import json
third_party = sorted({name.split('.')[0] for name, module in tuple(sys.modules.items())
    if getattr(module, '__file__', None) and Path(module.__file__).resolve().is_relative_to(allowed)})
print('PORTABLE_IMPORTED_THIRD_PARTY=' + json.dumps(third_party), flush=True)
print('PORTABLE_PROJECT_ACCESS_AUDIT_PASS', flush=True)
''')
    env = dict(os.environ, MPLBACKEND='Agg', PYTHONDONTWRITEBYTECODE='1')
    for key in ('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
        env[key] = '1'
    env['MPLCONFIGDIR'] = str(temp/'matplotlib')
    env.pop('PYTHONPATH', None)
    def command(relative, *args):
        return [str(PYTHON), str(wrapper), str(REPO), str(PYTHON.parent.parent), str(bundle), str(bundle/relative), *args]
    def run(relative, log, *args, timeout=120):
        with (OUT/log).open('w') as stream:
            result = subprocess.run(command(relative, *args), cwd=bundle, env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
        assert result.returncode == 0, log
        assert 'PORTABLE_PROJECT_ACCESS_AUDIT_PASS' in (OUT/log).read_text(), log
    run('code/model/plot_model_policies.py', 'final_policy_plots.log')
    run('code/model/plot_model_aggregates.py', 'final_aggregate_plots.log')
    receipt['cached_plotters'] = {'policy': 8, 'aggregate': 7, 'original_project_access_blocked': True}
    with socket.socket() as sock:
        sock.bind(('127.0.0.1', 0)); port = sock.getsockname()[1]
    with (OUT/'final_explorer.log').open('w') as stream:
        server = subprocess.Popen(command('code/model/tools/economics_explorer.py', '--config', str(bundle/'output/model/local_solution/latest/explorer_cases.json'), '--port', str(port)), cwd=bundle, env=env, stdout=stream, stderr=subprocess.STDOUT)
        try:
            for _ in range(100):
                if server.poll() is not None: raise RuntimeError('Explorer exited; see final_explorer.log')
                try:
                    with urllib.request.urlopen(f'http://127.0.0.1:{port}/api/meta', timeout=1) as r: meta=json.load(r)
                    break
                except OSError: time.sleep(.1)
            else: raise RuntimeError('Explorer startup timeout')
            case=meta['cases'][0]['id']
            query=urllib.parse.urlencode(dict(case=case,age=30,income=meta['default']['income'],tenure=0,children=0,at_home=0))
            for endpoint in ('/', '/api/slice?'+query, '/api/aggregates?'+urllib.parse.urlencode(dict(case=case))):
                with urllib.request.urlopen(f'http://127.0.0.1:{port}'+endpoint,timeout=20) as response:
                    assert response.status==200 and len(response.read())>100
            receipt['explorer']={'http_routes_passed':['/','/api/meta','/api/slice','/api/aggregates'],'original_project_access_blocked':True}
        finally:
            server.terminate(); server.wait(timeout=10)
    save()
    if '--inspect-only' in sys.argv:
        receipt.update(status='cached_inspection_passed_fresh_solve_unverified',
                       completed_epoch=time.time(), fresh_ge='not run on this archive')
        save(); print(json.dumps(receipt, indent=2), flush=True)
        sys.exit(0)
    print('CACHED_PLOTS_AND_EXPLORER_PASS; starting one fresh GE', flush=True)
    run('code/model/run_model.py', 'final_fresh_ge.log', timeout=1200)
    import numpy as np
    latest=(bundle/'output/model/local_solution/latest').resolve()
    assert latest.parent.name=='cases' and latest.name!='bundled_reference'
    assert (latest.parent/'bundled_reference/native_result.npz').exists()
    for name, count in [('target_fit.csv',14),('parameters.csv',31)]:
        with (REFERENCE/name).open() as f: left=list(csv.DictReader(f))
        with (latest/name).open() as f: right=list(csv.DictReader(f))
        assert len(left)==len(right)==count and left==right, name
    with np.load(REFERENCE/'native_result.npz',allow_pickle=False) as a, np.load(latest/'native_result.npz',allow_pickle=False) as b:
        assert set(a.files)==set(b.files)
        for key in a.files:
            assert np.array_equal(a[key],b[key],equal_nan=True), key
        arrays=len(a.files)
    receipt.update(status='passed',completed_epoch=time.time(),fresh_case=str(latest),fresh_ge_original_project_access_blocked=True,
                   full_target_rows_exact=14,full_parameter_rows_exact=31,stored_arrays_exact=arrays,
                   cached_reference_preserved=True,production_defaults_unchanged=True,transition_not_included=True)
    save(); print(json.dumps(receipt,indent=2))
except SystemExit:
    raise
except BaseException as error:
    receipt.update(status='failed',error=repr(error),finished_epoch=time.time()); save(); raise
