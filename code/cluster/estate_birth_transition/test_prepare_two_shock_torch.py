"""Filesystem-only packaging checks; these tests never import a native model."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('packaging',HERE/'prepare_two_shock_torch.py')
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)


class PackagingTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name).resolve()
    def tearDown(self):self.tmp.cleanup()
    def put(self,path,blob=b'bytes'):
        path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(blob);return p.sha(path)
    def test_safe_relative_paths(self):
        for bad in ('/tmp/x','../x','a/../../x','.'):
            with self.assertRaises(ValueError):p.safe_rel(bad)
    def test_authentication_rejects_drift_and_symlinks(self):
        digest=self.put(self.root/'x')
        p.authenticate_files(self.root,{'x':digest})
        (self.root/'link').symlink_to(self.root/'x')
        with self.assertRaises(ValueError):p.authenticate_files(self.root,{'link':digest})
        (self.root/'x').write_bytes(b'drift')
        with self.assertRaises(ValueError):p.authenticate_files(self.root,{'x':digest})
    def test_exact_copy_is_readonly(self):
        digest=self.put(self.root/'src/a/file')
        p.copy_exact(self.root/'src',self.root/'dest',{'a/file':digest})
        self.assertEqual(p.sha(self.root/'dest/a/file'),digest)
        self.assertEqual((self.root/'dest/a/file').stat().st_mode&0o222,0)
    def test_relocate_plan_only_declared_paths(self):
        old=self.root/'canonical';new=self.root/'frozen/source'
        pin=dict(path=str(old/'x'),sha256=self.put(new/'x'))
        plan=dict(source_files=dict(runtime=pin),handoff=pin,target_contract=dict(annual=pin,blocks=pin,rows=[1]),identity={'stable':1},fit={'tolerance':.005})
        plan=json.loads(json.dumps(plan))
        result=p.relocate_plan(plan,old,new)
        self.assertEqual(result['identity'],plan['identity']);self.assertEqual(result['fit'],plan['fit'])
        self.assertEqual(result['source_files']['runtime']['path'],str(new/'x'))
        self.assertEqual(plan['handoff']['path'],str(old/'x'))
        (new/'x').write_bytes(b'changed')
        with self.assertRaises(ValueError):p.relocate_plan(plan,old,new)
    def test_relocate_overlay_preserves_original_namespace(self):
        canonical=self.root/'canonical';snapshot=self.root/'old/files';destination=self.root/'fresh'
        files={f'code/{n}.py':'abc' for n in range(1241)}
        metadata=dict(schema='exact_read_only_source_overlay_v1',manifest_sha256=p.EXPECTED_MANIFEST,
          all_source_hashes_verified=True,original_root=str(canonical),snapshot_root=str(snapshot),
          mapping={str(canonical/rel):str(snapshot/rel) for rel in files},manifest_file_count=1241)
        result=p.relocate_overlay(metadata,{'files':files},destination,canonical)
        self.assertEqual(result['original_root'],str(canonical));self.assertEqual(result['snapshot_root'],str(destination/'frozen/source_overlay/files'))
        self.assertEqual(result['manifest_path'],str(destination/'frozen/source'/p.MANIFEST_REL))
        self.assertEqual(metadata['snapshot_root'],str(snapshot))
        metadata['mapping']['extra']='x'
        with self.assertRaises(ValueError):p.relocate_overlay(metadata,{'files':files},destination,canonical)
    def test_old_overlay_authentication_resolves_legacy_symlink(self):
        canonical=self.root/'canonical';snapshot=self.root/'old/files';destination=self.root/'fresh'
        archived=canonical/'calibration_archive/legacy';archived.mkdir(parents=True)
        alias=canonical/'code/model/legacy';alias.parent.mkdir(parents=True)
        alias.symlink_to(archived,target_is_directory=True)
        files={f'code/model/legacy/{n}.py':'abc' for n in range(1241)}
        metadata=dict(schema='exact_read_only_source_overlay_v1',manifest_sha256=p.EXPECTED_MANIFEST,
          all_source_hashes_verified=True,original_root=str(canonical),snapshot_root=str(snapshot),
          mapping={str((canonical/rel).resolve()):str((snapshot/rel).resolve()) for rel in files},manifest_file_count=1241)
        result=p.relocate_overlay(metadata,{'files':files},destination,canonical)
        # The old recorded map must authenticate its resolved alias; the fresh
        # map retains lexical canonical names for resolution in its own mount.
        rel='code/model/legacy/0.py'
        self.assertIn(str(archived/'0.py'),metadata['mapping'])
        self.assertEqual(result['mapping'][str(canonical/rel)],str(destination/'frozen/source_overlay/files'/rel))
        wrong=dict(metadata);wrong['mapping']=dict(metadata['mapping'])
        wrong['mapping'][str(archived/'0.py')]=str(snapshot/'wrong.py')
        with self.assertRaisesRegex(ValueError,'Overlay path mapping differs'):
            p.relocate_overlay(wrong,{'files':files},destination,canonical)
    def test_existing_destination_fails_before_reads(self):
        destination=self.root/'existing';destination.mkdir()
        args=argparse.Namespace(repo=self.root,destination=destination,remote_root='/scratch/new',python='python')
        with self.assertRaisesRegex(ValueError,'Fresh destination'):p.build(args)
        self.assertEqual(list(destination.iterdir()),[])
    def test_inventory_identity_failure_leaves_no_destination(self):
        repo=self.root/'repo';base=repo/p.BASE_REL
        self.put(base/'deployment_v9/inventory.json',b'{}')
        args=argparse.Namespace(repo=repo,destination=self.root/'fresh',remote_root='/scratch/new',python='python')
        with self.assertRaisesRegex(ValueError,'inventory SHA'):p.build(args)
        self.assertFalse(Path(args.destination).exists())
    def test_prepare_uses_fresh_process_and_staged_driver(self):
        destination=self.root/'fresh'
        with patch.object(p.subprocess,'run') as run:p.prepare_plans(destination,'/python')
        args=run.call_args[0][0]
        self.assertEqual(args[0],'/python');self.assertEqual(args[1],'-c')
        self.assertEqual(args[-2],str(destination/'frozen/source'/Path(p.DRIVER_REL).parent))
        self.assertIn('import two_shock as d',args[2]);self.assertIn('d.prepare_manifest',args[2])
        self.assertEqual(run.call_args.kwargs['timeout'],120)
        self.assertEqual(run.call_args.kwargs['env']['NUMBA_NUM_THREADS'],'1')
    def test_launcher_shell_syntax_and_bind_namespace(self):
        stage=self.root/'space package';remote=Path('/scratch/example');canonical=self.root/'canonical'
        text=p.launcher_text(stage,remote,canonical);script=self.root/'launcher.sh';script.write_text(text)
        subprocess.run(['bash','-n',str(script)],check=True)
        self.assertIn('--bind "$remote:$local_root:ro"',text)
        self.assertIn('--bind "$remote/frozen/source:$repo:ro"',text)
        self.assertIn('--bind "$remote/results:$local_root/results:rw"',text)
        self.assertIn('runtime=r.build_runtime(plan=plan,output=out,smoke=True)',text)
        self.assertIn('runtime.rt.total_native_calls==0',text)
        self.assertIn('#SBATCH --cpus-per-task=1',text);self.assertIn('#SBATCH --mem=24G',text)
        self.assertIn('--run --smoke-receipt-pin "$smoke_pin"',text)
        self.assertNotIn('#SBATCH --array',text)
        self.assertIn('mkdir "$claim" ||',text);self.assertIn('mkdir "$job" ||',text);self.assertIn('kill -KILL -- "-$child_pid"',text)
    def test_launcher_budgets_and_single_threads(self):
        text=p.launcher_text(self.root/'stage',Path('/scratch/new'),self.root/'repo')
        self.assertIn('preflight) seconds=600;; smoke) seconds=1800;; fit) seconds=21600',text)
        for value in ('NUMBA_NUM_THREADS=1','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1','VECLIB_MAXIMUM_THREADS=1','NUMEXPR_NUM_THREADS=1'):
            self.assertIn(value,text)
        self.assertIn('host_verification.json',text);self.assertIn('container_verification.json',text)
    def test_background_bounded_command_retains_heredoc(self):
        text=p.launcher_text(self.root/'stage',Path('/scratch/new'),self.root/'repo')
        function=text[text.index('run_bounded() {'):text.index('run_bounded "$python"')]
        fake=self.root/'bin';fake.mkdir()
        for name,body in [('setsid','exec "$@"'),('timeout','shift; shift; shift; exec "$@"')]:
            script=fake/name;script.write_text('#!/bin/bash\n'+body+'\n');script.chmod(0o755)
        target=self.root/'captured'
        shell='PATH='+str(fake)+':$PATH\ndeadline_epoch=$(($(date +%s)+30))\n'+function+"\nrun_bounded bash -c 'cat > \"$1\"' bash "+str(target)+" <<'DATA'\nconstructor input\nDATA\n"
        subprocess.run(['bash','-c',shell],check=True)
        self.assertEqual(target.read_text(),'constructor input\n')
    def test_global_mode_claim_refuses_different_name(self):
        text=p.launcher_text(self.root/'stage',Path('/scratch/new'),self.root/'repo')
        first=text.index('if [[ "$mode" != preflight ]]; then\n claim=')
        claim=text[first:text.index('out="$remote/results/',first)]
        remote=self.root/'remote';(remote/'jobs').mkdir(parents=True)
        command='set -euo pipefail\nremote='+str(remote)+'\nmode=smoke\nname=$1\nSLURM_JOB_ID=1\n'+claim
        subprocess.run(['bash','-c',command,'bash','first'],check=True,capture_output=True)
        second=subprocess.run(['bash','-c',command,'bash','second'],capture_output=True,text=True)
        self.assertEqual(second.returncode,2)
        self.assertIn('Refusing duplicate or unknown',second.stdout)
        self.assertEqual((remote/'jobs/smoke.claim/run_name').read_text(),'first\n')
        self.assertFalse((remote/'jobs/smoke_second').exists())
    def test_verify_rejects_uninventoried_control_drift(self):
        stage=self.root/'stage';control='inputs/fit_manifest.json'
        digest=self.put(stage/control,b'{}')
        p.dump(stage/'inventory.json',dict(schema=p.SCHEMA,files={control:digest}))
        self.put(stage/control,b'{"tampered":true}')
        with self.assertRaisesRegex(ValueError,'Authenticated file mismatch'):p.verify(stage)

if __name__=='__main__':unittest.main()
