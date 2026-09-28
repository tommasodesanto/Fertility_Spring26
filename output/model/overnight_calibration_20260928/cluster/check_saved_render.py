#!/usr/bin/env python3
"""Torch-only real saved-checkpoint render gate. Zero equilibrium solves."""
from __future__ import annotations
import argparse,hashlib,importlib.util,json,os,signal,sys,time,traceback
from pathlib import Path
CONTRACT_SHA='c83aaff1a90b1ba5bb0919151840e745e6816cb0a3a023ce746e750f77226c5d'
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS')

def read(path):return json.loads(Path(path).read_text())
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module);return module

def main():
    p=argparse.ArgumentParser();p.add_argument('--contract',type=Path,required=True);p.add_argument('--smoke',type=Path,required=True);p.add_argument('--smoke-sha256',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm only'
    # Inspect the completed six-case packet before importing any runtime/code.
    smoke=read(a.smoke)
    assert sha(a.smoke)==a.smoke_sha256 and smoke['status']=='exact_loop_smoke_passed'
    assert len(smoke['records'])==6 and all(row['status']=='success' for row in smoke['records'])
    assert sha(a.contract)==smoke['contract_sha256']==CONTRACT_SHA
    declared=read(a.contract);pin=declared['files']['driver'];assert sha(pin['path'])==pin['sha256']
    os.environ['EXPECTED_E5F_NIGHT_SHA256']=CONTRACT_SHA
    for key in THREADS:os.environ[key]='1'
    os.environ['MPLBACKEND']='Agg'
    driver=load('saved_render_verified_driver',pin['path']);c,objs=driver.verify(a.contract)
    verified=driver.verified_smoke(c,objs,a.smoke,a.smoke_sha256,CONTRACT_SHA)
    row=next(r for r in verified['records'] if r['lane']=='primary')
    source=Path(row['case_path']);names=sorted(c['standard_diagnostic_names'])
    assert len(driver.table(source/'target_fit.csv','moment'))==14
    assert len(driver.table(source/'parameters.csv','parameter'))==31
    source_pins={name:sha(source/name) for name in ('receipt.json','initial_state.pkl.gz','target_fit.csv','parameters.csv','stationary_solves.json')}
    source_png={name:sha(source/'standard_diagnostics'/name) for name in names}
    a.output.mkdir(parents=True,exist_ok=False)
    child=a.output/'render_check';now=time.time();deadline=min(now+300,c['budget']['absolute_end_epoch'])
    assert deadline>now
    context=dict(candidate_id=child.name,stage='render',contract_sha256=CONTRACT_SHA,source_sha256=pin['sha256'],target_sha256=c['lanes'][row['lane']]['objective']['sha256'],point_sha256=driver.canon(row['point']))
    request=dict(id=child.name,lane=row['lane'],point=row['point'],context=context,design='real_saved_checkpoint_render_gate',saved_case=row,stage='render',contract_sha256=CONTRACT_SHA,scientific_candidate_id=driver.identity(c,row['lane'],row['point']),normalization_inputs=c['normalization'],controller_pid=os.getpid(),deadline_epoch=deadline,graphs=True)
    request_path=a.output/'request.json';driver.write(request_path,request)
    supervisor_path=c['files']['recovery_search']['path'];sys.path.insert(0,str(Path(supervisor_path).parent));supervisor=load('saved_render_owned_supervision',supervisor_path)
    command=[sys.executable,pin['path'],'--stage','render','--contract',str(a.contract),'--output',str(child),'--request',str(request_path)]
    process=None;started=time.monotonic();last=0.
    def interrupted(signum,frame):raise InterruptedError('Render-check parent interrupted')
    signal.signal(signal.SIGTERM,interrupted)
    try:
        process=supervisor.ManagedProcess(command,a.output/'render.log',deadline,os.environ.copy())
        while True:
            code=process.poll()
            if code is not None:break
            if time.time()-last>=5:
                driver.write(a.output/'heartbeat.json',dict(status='rendering_saved_checkpoint',pid=process.process.pid,elapsed_seconds=time.monotonic()-started,deadline_epoch=deadline,objective_evaluations=0));last=time.time()
            time.sleep(.5)
        assert code==0 and not process.deadline_expired,'Render subprocess failed or exceeded300-second deadline'
        receipt=read(child/'render_receipt.json');startup=read(child/'startup.json')
        assert receipt['status']=='saved_checkpoint_rendered' and receipt['context']==context and receipt['source']==row and receipt['objective_evaluations']==0
        assert startup['context']==context and startup['pid']==process.process.pid and startup['parent_pid']==os.getpid() and startup['normalization_inputs']==c['normalization']
        case=Path(receipt['case_path']);assert case.resolve()==(child/'case').resolve()
        table_check=driver.compare_tables(source,case)
        assert len(driver.table(case/'target_fit.csv','moment'))==14 and len(driver.table(case/'parameters.csv','parameter'))==31
        assert sha(case/'target_fit.csv')==source_pins['target_fit.csv'] and sha(case/'parameters.csv')==source_pins['parameters.csv']
        assert sorted(p.name for p in (case/'standard_diagnostics').glob('*.png'))==names
        rendered_png={name:sha(case/'standard_diagnostics'/name) for name in names}
        from PIL import Image
        sizes={}
        for name in names:
            with Image.open(case/'standard_diagnostics'/name) as im:sizes[name]=list(im.size);im.verify()
        assert all(width>0 and height>0 for width,height in sizes.values())
        assert rendered_png==source_png,'Real regenerated PNGs differ from the authenticated smoke packet'
        assert {name:sha(source/name) for name in source_pins}==source_pins,'Saved scientific source changed'
        driver.verify(a.contract)
        result=dict(status='real_saved_render_passed',contract_sha256=CONTRACT_SHA,smoke_sha256=a.smoke_sha256,source_case=str(source),source_pins=source_pins,table_check=table_check,parameter_rows=31,target_rows=14,png_count=17,png_byte_identical_to_smoke=True,png_sha256=rendered_png,png_sizes=sizes,objective_evaluations=0,child_pid=process.process.pid,elapsed_seconds=time.monotonic()-started,helper_sha256=sha(__file__),command=command)
        driver.write(a.output/'complete.json',result)
        print(json.dumps({key:result[key] for key in ('status','contract_sha256','target_rows','parameter_rows','png_count','png_byte_identical_to_smoke','objective_evaluations','elapsed_seconds')}))
    except BaseException as exc:
        driver.write(a.output/'complete.json',dict(status='real_saved_render_failed',contract_sha256=CONTRACT_SHA,error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc(),objective_evaluations=0,elapsed_seconds=time.monotonic()-started));raise
    finally:
        if process is not None:process.close()

if __name__=='__main__':main()
