"""Filesystem-only packaging checks; these tests never import a native model."""
import argparse
import hashlib
import importlib.util
import importlib.machinery
import json
import os
import sys
import time
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
    def test_reference_graph_follows_explicit_constructor_pins(self):
        root=self.root/'repo';case=root/'original_case';pair=root/'portable/pair';export=root/'export';docs={}
        def item(path,value=None):
            record=dict(path=str(path),sha256=hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest())
            if value is not None:docs[str(path)]=value
            return record
        artifact=item(export/'table.csv')
        source=item(root/'native_inventory.json',dict(files={'code/native.py':'native-sha'}))
        base=item(root/'base.json',dict(reference_root=str(pair),parent_lock={'path':'/obsolete/scratch','sha256':'unused'}))
        native=item(root/'native.json',dict(files={'tool':item(root/'tool.py')},base_contract=base,
             objective=item(root/'objective.json',{}),source_manifest=source,source_root=str(root/'native_source')))
        contract=item(root/'contract.json',dict(files={'builder':item(root/'builder.py')},reference_case=str(case)))
        manifest=dict(contract=contract,objective=item(root/'objective.json',{}),source_manifest=item(root/'source_manifest.json',{}),
           source_contract=item(root/'source_contract.json',{}),native_ancestry_contract=native,
           local_export=str(export),artifact_hashes={'table.csv':artifact['sha256']},checkpoint={'sha256':'reference-checkpoint'})
        docs[str(pair/'inputs/launch_lock.json')]=dict(runtime_file_sha256={'run_pair.py':'pair-driver'},
          objective_sha256='pair-objective',proposal_bank_sha256='bank',source_manifest_sha256='pair-source',ancestor_sha256='ancestor',tax_driver_sha256='tax')
        docs[str(pair/'inputs/objective.json')]={};docs[str(pair/'inputs/proposal_bank.json')]={}
        docs[str(pair/'inputs/source_manifest.json')]=dict(source_root=str(pair/'source'),files={'code/pair.py':'pair-code'})
        files,provenance=p.reference_input_pins(manifest,root,lambda pin:docs[pin['path']])
        self.assertEqual(files['contract.json'],contract['sha256'])
        self.assertEqual(files['export/initial_state.pkl.gz'],'reference-checkpoint')
        self.assertEqual(files['native_source/code/native.py'],'native-sha')
        self.assertEqual(files['portable/pair/source/code/pair.py'],'pair-code')
        self.assertEqual(files['original_case/initial_state.pkl.gz'],p.CASE_CHECKPOINT_SHA)
        for name,digest in p.CASE_LEAVES.items():
            self.assertEqual(files['original_case/'+name],digest)
            self.assertIn(p.REVIEW_MANIFEST_SHA,provenance['original_case/'+name][0])
        self.assertNotIn('/obsolete/scratch',files)
    def test_reference_resolver_rejects_missing_and_drifted_bytes(self):
        repo=self.root/'repo';frozen=self.root/'frozen';overlay=self.root/'overlay'
        manifest_path=frozen/p.REFERENCE_REL;digest=self.put(manifest_path,b'{}')
        wanted=hashlib.sha256(b'exact').hexdigest()
        with patch.object(p,'EXPECTED_INPUTS',{p.REFERENCE_REL:digest}),patch.object(p,'reference_input_pins',return_value=({'required.json':wanted},{'required.json':['declared original pin']})):
            with self.assertRaisesRegex(ValueError,'No exact authenticated reference input'):
                p.authenticated_reference_sources(repo,frozen,overlay,{}, {})
            self.put(repo/'required.json',b'drift')
            with self.assertRaisesRegex(ValueError,'No exact authenticated reference input'):
                p.authenticated_reference_sources(repo,frozen,overlay,{}, {})
            self.put(overlay/'historical/exact.json',b'exact')
            files,authority,chosen=p.authenticated_reference_sources(repo,frozen,overlay,{}, {'historical/exact.json':wanted})
            self.assertEqual(chosen['required.json'],overlay/'historical/exact.json')
            self.assertEqual(files['required.json'],wanted)
    def test_overlay_materialization_preserves_priority_and_pins(self):
        source=self.root/'source';snapshot=self.root/'snapshot'
        files={name:self.put(snapshot/name,blob) for name,blob in
          [('missing.py',b'missing overlay'),('same.py',b'same bytes'),('different.py',b'old overlay')]}
        self.put(source/'same.py',b'same bytes');native=self.put(source/'different.py',b'protected native')
        receipt=p.materialize_overlay_sources(source,snapshot,files)
        self.assertEqual(receipt['counts'],dict(materialized_missing=1,existing_identical=1,existing_different_overlay_precedence=1))
        self.assertEqual(p.sha(source/'different.py'),native)
        self.assertEqual(p.sha(source/'missing.py'),files['missing.py'])
        self.assertEqual(receipt['existing_files_overwritten'],0)
        # Builder's final freeze makes retained inputs read-only as well.
        for path in source.iterdir():path.chmod(0o444)
        p.verify_overlay_materialization(source,files,receipt)
        (source/'missing.py').chmod(0o644);(source/'missing.py').write_bytes(b'drift')
        with self.assertRaisesRegex(ValueError,'physical/source pin differs'):
            p.verify_overlay_materialization(source,files,receipt)
    def test_overlay_materialization_enables_module_discovery_without_import(self):
        source=self.root/'source';tools=source/'code/model/tools';tools.mkdir(parents=True)
        snapshot=self.root/'snapshot';rel='code/model/tools/run_e5f_due_stayer_matched_check.py'
        # If accidentally imported, this fixture deliberately fails. Finding its
        # source path is the specific operation that failed in the real package.
        digest=self.put(snapshot/rel,b'raise RuntimeError("must not import native code")\n')
        self.assertIsNone(importlib.machinery.PathFinder.find_spec('run_e5f_due_stayer_matched_check',[str(tools)]))
        receipt=p.materialize_overlay_sources(source,snapshot,{rel:digest})
        spec=importlib.machinery.PathFinder.find_spec('run_e5f_due_stayer_matched_check',[str(tools)])
        self.assertIsNotNone(spec);self.assertEqual(Path(spec.origin),tools/'run_e5f_due_stayer_matched_check.py')
        self.assertEqual(receipt['counts']['materialized_missing'],1)
        p.verify_overlay_materialization(source,{rel:digest},receipt)
    def test_overlay_materialization_rejects_missing_snapshot(self):
        source=self.root/'source';snapshot=self.root/'snapshot'
        with self.assertRaisesRegex(ValueError,'Authenticated file mismatch'):
            p.materialize_overlay_sources(source,snapshot,{'missing.py':'expected-sha'})
        self.assertFalse(source.exists())
    def test_verify_rejects_uninventoried_control_drift(self):
        stage=self.root/'stage';control='inputs/fit_manifest.json'
        digest=self.put(stage/control,b'{}')
        p.dump(stage/'inventory.json',dict(schema=p.SCHEMA,files={control:digest}))
        self.put(stage/control,b'{"tampered":true}')
        with self.assertRaisesRegex(ValueError,'Authenticated file mismatch'):p.verify(stage)

class ConstructorCompletionTests(unittest.TestCase):
    setUp=PackagingTests.setUp
    tearDown=PackagingTests.tearDown
    put=PackagingTests.put
    def fixture(self):
        stage=self.root/'stage';stage.mkdir();(stage/'jobs').mkdir();(stage/'results').mkdir()
        self.stage=stage;self.dest=self.root/'completion';proof='results/preflight_preflight_v1'
        plan=dict(smoke=True,schema='fixture',kind='two_unanticipated_permanent',identity={'unchanged':'identity'},
            baseline_psi=.1,psi_bound_ratios=[.01,2.],stages=[1,2],rows=[1,2,3,4],weights=[0,1,0,1],
            target_contract={'unchanged':'targets'},legacy_source_overlay={'unchanged':'overlay'},
            gates={},seed={},fit={},endpoint={},path={},budget={'total_seconds':1680,'maximum_policy_calls':400},
            horizons=[6],smoke_seed_endpoint_padding=True,empirical_controls={'total_seconds':21480},
            smoke_protocol='native_two_stage_execution_only_v2')
        pins={}
        for rel in ('prepare_two_shock_torch.py','inputs/fit_manifest.json','inputs/base_fit_plan.json',
            'frozen/source/'+p.DRIVER_REL,'frozen/source/'+p.RUNTIME_REL,'frozen/source/'+p.HELPER_REL,
            'frozen/source_overlay/overlay.json'):
            pins[rel]=self.put(stage/rel,b'authenticated immutable fixture bytes')
        pins['launch_torch.sh']=self.put(stage/'launch_torch.sh',p.launcher_text(stage,stage,self.root/'repo').encode())
        plan['source_pins']={name:dict(path=str(stage/'frozen/source'/rel),sha256=pins['frozen/source/'+rel])
            for name,rel in [('two_shock_driver',p.DRIVER_REL),('two_shock_runtime',p.RUNTIME_REL)]}
        p.dump(stage/'inputs/smoke_manifest.json',plan);pins['inputs/smoke_manifest.json']=p.sha(stage/'inputs/smoke_manifest.json')
        pins.update({f'unscanned/{i}':'not accessed' for i in range(4599-len(pins))})
        inv=dict(schema=p.SCHEMA,local_root=str(stage),remote_root=str(stage),canonical_root=str(self.root/'repo'),files=pins,identity=plan['identity'])
        p.dump(stage/'inventory.json',inv);digest=p.sha(stage/'inventory.json')
        self.patches=[patch.object(p,'COMPLETION_INVENTORY',digest),patch.object(p,'COMPLETION_REMOTE',str(stage))]
        for mock in self.patches:mock.start();self.addCleanup(mock.stop)
        start=dict(mode='preflight',wall_seconds=600,inventory_sha256=digest,cpus=1,memory_gib=24,
            numba_threads=1,blas_threads=1,no_auto_retry=True,start_epoch=1000,deadline_epoch=1600,
            output=str(stage/proof),slurm_job_id=None)
        terminal=dict(mode='preflight',exit_code=124,no_auto_retry=True,start_epoch=1000,deadline_epoch=1600,finished_epoch=1588,slurm_job_id=None)
        package=dict(status='PASS_ZERO_SOLVES_PACKAGE',files=4599,native_calls=0,model_solves=0,scientific_validation=False)
        driver=dict(status='PASS',native_calls=0,scientific_validation=False,production_ready=False,schema=plan['schema'],
            fingerprints=p.completion_fingerprints(plan),horizons=[6],total_seconds=1680,policy_call_stop_cap=400)
        for name,value in zip(p.COMPLETION_PROOFS,(start,terminal,package,package,driver)):p.dump(stage/proof/name,value)
        self.args=argparse.Namespace(stage=stage,destination=self.dest,remote_packet=self.dest,
            proof_relative=proof,name='constructor_completion_v1')
        return plan

    def prepare(self):
        self.fixture();return p.prepare_constructor_completion(self.args)

    def test_completion_packet_exact_body_small_and_immutable(self):
        result=self.prepare();original=p.sha(self.stage/'inventory.json')
        m=p.verify_constructor_completion(self.dest,result['inventory_sha256'])
        self.assertEqual((self.dest/'constructor.py').read_text(),p.constructor_body((self.stage/'launch_torch.sh').read_text()))
        self.assertEqual(m['original_overall_preflight'],'FAILED_EXIT_124')
        self.assertEqual(p.sha(self.stage/'inventory.json'),original)
        self.assertLess(sum(x.stat().st_size for x in self.dest.rglob('*') if x.is_file()),100000)
        subprocess.run(['bash','-n',str(self.dest/'launch_constructor.sh')],check=True)
        with self.assertRaisesRegex(ValueError,'Fresh destination'):p.prepare_constructor_completion(self.args)

    def test_completion_rejects_bad_proof_before_creation(self):
        self.fixture()
        for name,key,value in [('launcher_terminal.json','exit_code',0),('host_verification.json','native_calls',1),
            ('container_verification.json','files',4598),('run/preflight.json','fingerprints',{})]:
            path=self.stage/self.args.proof_relative/name;original=path.read_bytes();data=p.read(path);data[key]=value;p.dump(path,data)
            with self.assertRaises(ValueError):p.prepare_constructor_completion(self.args)
            self.assertFalse(self.dest.exists());path.write_bytes(original)
        path=self.stage/self.args.proof_relative/'host_verification.json';path.unlink()
        with self.assertRaises(FileNotFoundError):p.prepare_constructor_completion(self.args)
        self.assertFalse(self.dest.exists())

    def test_completion_rejects_source_and_receipt_mutation(self):
        result=self.prepare()
        for path in (self.stage/'frozen/source'/p.DRIVER_REL,self.stage/self.args.proof_relative/'host_verification.json',self.dest/'manifest.json'):
            original=path.read_bytes();path.chmod(0o644);path.write_bytes(b'mutated')
            with self.assertRaises(ValueError):p.verify_constructor_completion(self.dest,result['inventory_sha256'])
            path.write_bytes(original)
        with self.assertRaises(ValueError):p.verify_constructor_completion(self.dest,'wrong')

    def test_constructor_unique_markers_and_strict_checks(self):
        text=p.launcher_text(self.root/'local',self.root/'remote',self.root/'repo')
        for changed in (text+"\nPYNATIVE\n",text.replace('assert runtime.rt.total_native_calls==0','pass'),text.replace("<<'PYNATIVE'","<<'OTHER'")):
            with self.assertRaises(ValueError):p.constructor_body(changed)

    def test_completion_harmless_process_cleanup_and_bound(self):
        marker=self.root/'leaked'
        program='import os,time,pathlib,sys\nchild=os.fork()\nif child==0:\n time.sleep(.5)\n pathlib.Path(sys.argv[1]).write_text("leaked")\nelse:\n time.sleep(10)\n'
        with open(os.devnull,'rb') as source,open(os.devnull,'wb') as output:
            status=p.bounded_completion_process([sys.executable,'-c',program,str(marker)],environment=os.environ.copy(),stdin=source,stdout=output,seconds=.1)
            self.assertEqual(status,124)
            with self.assertRaises(ValueError):p.bounded_completion_process([],environment={},stdin=source,stdout=output,seconds=586)
        time.sleep(.6);self.assertFalse(marker.exists())

    def runtime(self,native_calls=0,write_native=True,status=0):
        import resource
        result=self.prepare()
        def fake(command,**kwargs):
            self.assertEqual(command[-3],str(self.stage/'frozen/source'/p.DRIVER_REL))
            self.assertLessEqual(kwargs['seconds'],585)
            self.assertEqual(kwargs['environment']['NUMBA_NUM_THREADS'],'1')
            if write_native:p.dump(self.stage/'results/constructor_completion_v1/native_constructor/native_status.json',
                dict(status='PASS_ZERO_SOLVES_NATIVE_CONSTRUCTOR',native_calls=native_calls,numba_threads=1,
                identity={'unchanged':'identity'},scientific_validation=False))
            return status
        with patch.object(p,'bounded_completion_process',side_effect=fake),patch.object(resource,'setrlimit'):
            args=argparse.Namespace(packet=self.dest,inventory_sha256=result['inventory_sha256'])
            if native_calls or not write_native or status:
                with self.assertRaises((ValueError,FileNotFoundError)):p.run_constructor_completion(args)
            else:p.run_constructor_completion(args)
        return args

    def test_runtime_actual_status_duplicate_guard_and_receipts(self):
        args=self.runtime();out=self.stage/'results/constructor_completion_v1'
        terminal=p.read(out/'launcher_terminal.json')
        self.assertTrue(terminal['composite_ready']);self.assertEqual(terminal['exit_code'],0)
        self.assertEqual(terminal['original_overall_preflight'],'FAILED_EXIT_124')
        with self.assertRaises(FileExistsError):p.run_constructor_completion(args)

    def test_runtime_rejects_native_calls(self):
        self.runtime(native_calls=1)
        terminal=p.read(self.stage/'results/constructor_completion_v1/launcher_terminal.json')
        self.assertFalse(terminal['composite_ready']);self.assertEqual(terminal['exit_code'],1)

    def test_runtime_missing_native_status_fails(self):
        self.runtime(write_native=False)
        self.assertFalse(p.read(self.stage/'results/constructor_completion_v1/launcher_terminal.json')['composite_ready'])

    def test_runtime_timeout_preserves_failed_terminal_and_claim(self):
        self.runtime(write_native=False,status=124)
        terminal=p.read(self.stage/'results/constructor_completion_v1/launcher_terminal.json')
        self.assertEqual(terminal['exit_code'],124);self.assertFalse(terminal['composite_ready'])
        self.assertIsNone(terminal['native_status_pin'])
        self.assertTrue((self.stage/'jobs/constructor_completion.claim').exists())
        self.assertEqual(p.read(self.stage/self.args.proof_relative/'launcher_terminal.json')['exit_code'],124)

    def test_runtime_existing_output_never_overwritten(self):
        result=self.prepare();out=self.stage/'results/constructor_completion_v1';out.mkdir();(out/'keep').write_text('keep')
        with self.assertRaises(FileExistsError):p.run_constructor_completion(argparse.Namespace(packet=self.dest,inventory_sha256=result['inventory_sha256']))
        self.assertEqual((out/'keep').read_text(),'keep');self.assertTrue((self.stage/'jobs/constructor_completion.claim').exists())


if __name__=='__main__':unittest.main()
