"""Three bounded Slurm-array stages; authentication and original GN mathematics."""
import argparse, ast, importlib.util, json, math, os, subprocess, sys, time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
PLAN=json.loads((HERE/'plan.json').read_text())

def require(x,message):
    if not x:raise RuntimeError(message)

def load_runtime(path):
    sys.path.insert(0,str(path));s=importlib.util.spec_from_file_location('old_runner',path/'runner.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def authenticate(m,runtime,source_root,inventory_path,oldrun):
    require(m.inputs.sha(inventory_path)==PLAN['original_inventory_sha256'],'Original inventory drift')
    inventory=json.loads(inventory_path.read_text())
    tree=ast.parse((runtime/'phase_b_pilot.py').read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='validate_parameter_estimates')
    namespace={'_require':require};exec(compile(ast.Module(body=[fn],type_ignores=[]),str(runtime/'phase_b_pilot.py'),'exec'),namespace);validate=namespace['validate_parameter_estimates']
    for rel,h in inventory['files'].items():require(m.inputs.sha(source_root/rel)==h,'Original staged source drift '+rel)
    for name,h in PLAN['runtime_files'].items():require(m.inputs.sha(runtime/name)==h,'Runtime drift '+name)
    require(m.inputs.canonical(json.loads((runtime/'source_pins.json').read_text()))==PLAN['source_fingerprint'],'Source fingerprint drift')
    require(m.inputs.canonical(m.PLAN['target_contract'])==PLAN['target_fingerprint'],'Target fingerprint drift')
    contract=json.loads((oldrun/'input_contract.json').read_text());seed,bounds,_=m.inputs.seed_and_bounds('floor_s2');coordinates=m.inputs.parameters('floor_s2')
    require(contract['lane']=='floor_s2' and contract['seed']==seed and contract['bounds']==bounds,'Old center or bounds drift')
    require(contract['target_contract_sha256']==PLAN['target_fingerprint'] and contract['target_contract']==m.PLAN['target_contract'],'Old target drift')
    require(contract['free_coordinates']==list(coordinates),'Old coordinate order drift')
    rows=json.loads((oldrun/'cases.json').read_text());baseline=rows[0];twin=rows[1]
    require(baseline['label']=='000_baseline' and baseline['parameters']==seed and twin['parameters']==seed,'Baseline center drift')
    steps=m.step_sizes(seed,bounds);reused=[baseline,twin];hashes={}
    for i,k in enumerate(coordinates[:-1]):
        row=rows[i+2];point=dict(seed);h=steps[k]
        if point[k]+h>bounds[k][1]:h=-h
        point[k]+=h
        require(row['label']==f'{i+2:03d}_r0_probe_{k}' and row['parameters']==point,'Old derivative vector drift '+k)
        reused.append(row)
    for row in reused:
        require(row['status']=='passed','Unpassed old case')
        report=oldrun/row['label']/'phase_b_ge/selected_root';row['report']=str(report)
        rr=m.residual(m.readtable(report/'target_fit.csv')).tolist()
        require(rr==row['residual'],'Old residual drift '+row['label'])
        require(abs(sum(x*x for x in rr)-row['loss'])<1e-8,'Old loss drift')
        actual={r['parameter']:float(r['estimate']) for r in m.readtable(report/'parameters.csv')};expected=m.expected_parameters(row['parameters'],(120,9),'floor')
        require(len(actual)==len(expected)==31,'Old parameter table length drift');validate(dict(expected_parameters=expected),m.PLAN['reference_parameter_table'],actual)
        hashes[row['label']]=m.report_hashes(report)
        receipt=json.loads((oldrun/row['label']/'native_selected_repeat.json').read_text())
        require(receipt['status']=='exact_full_ge_repeat_passed' and receipt['target_rows']==14 and receipt['parameter_rows']==31,'Old native repeat absent')
        require(receipt['standard_plot_hashes']=={Path(n).name:h for n,h in hashes[row['label']].items() if n.endswith('.png')},'Old native PNG repeat drift')
    require(baseline['residual']==twin['residual'],'Old baseline residual repeat drift');m.compare_repeated(Path(baseline['report']),Path(twin['report']))
    return reused,seed,bounds,coordinates,steps,float(contract['starting_price']),hashes

def jacobian(m,baseline,probes,coordinates,steps,bounds):
    spans=np.asarray([bounds[k][1]-bounds[k][0] for k in coordinates]);seed=baseline['parameters'];hs=[]
    for k in coordinates:
        h=steps[k]
        if seed[k]+h>bounds[k][1]:h=-h
        hs.append(h)
    J=np.column_stack([(np.asarray(row['residual'])-baseline['residual'])/(h/spans[i]) for i,(row,h) in enumerate(zip(probes,hs))])
    sv=np.linalg.svd(J,compute_uv=False);rank=int(np.linalg.matrix_rank(J))
    identity=dict(status='round_center_local_jacobian' if rank==len(coordinates) else 'underidentified_local_jacobian',rank=rank,center_label=baseline['label'],center_parameters=seed,free_coordinates=list(coordinates),scored_targets=10,valid_probes=len(probes),fresh_at_selected=False,singular_values=sv.tolist(),normalized_coordinate_jacobian=J.tolist())
    points=[]
    if rank==len(coordinates):
        rr=np.asarray(baseline['residual']);ridge=max(float(sv[0]**2)*1e-4,1e-10)
        delta=np.linalg.solve(J.T@J+ridge*np.eye(len(coordinates)),-J.T@rr)
        trust=np.asarray([3*steps[k]/spans[i] for i,k in enumerate(coordinates)]);delta=np.clip(delta,-trust,trust)
        identity.update(ridge=ridge,normalized_trust=trust.tolist(),normalized_step=delta.tolist())
        for damping in (.5,.2,1.):points.append({k:float(np.clip(seed[k]+damping*delta[i]*spans[i],*bounds[k])) for i,k in enumerate(coordinates)})
    return identity,points

def run_stage(m,a,name,tasks,deadline,submitted):
    require(time.time()<deadline,'Global sixty-minute budget expired');folder=a.out/name;folder.mkdir()
    for t in tasks:t.update(deadline_epoch=deadline)
    path=folder/'tasks.json';m.write(path,dict(stage=name,tasks=tasks,deadline_epoch=deadline))
    remaining=max(1,int(deadline-time.time()));minutes=max(1,math.ceil(remaining/60))
    cmd=['sbatch','--parsable',f'--array=0-{len(tasks)-1}',f'--time={minutes}',f'--export=ALL,FAST_STAGE={name},FAST_TASK_FILE={path},FAST_RESULT_ROOT={a.out},FAST_DEADLINE_EPOCH={deadline}',str(a.launcher)]
    job=subprocess.check_output(cmd,text=True).strip().split(';')[0];require(job.isdigit(),'Invalid sbatch result');submitted.append(job)
    m.write(a.out/'jobs.json',dict(jobs=submitted,deadline_epoch=deadline));m.write(a.out/'latest.json',dict(status='array_running',stage=name,job=job,deadline_epoch=deadline))
    while True:
        if time.time()>=deadline:
            subprocess.run(['scancel',job],check=False);raise TimeoutError('Global sixty-minute budget reached; active array cancelled')
        receipts=[folder/str(i)/'case.json' for i in range(len(tasks))]
        if all(p.exists() for p in receipts):break
        active=subprocess.check_output(['squeue','-h','-j',job,'-o','%i'],text=True).strip()
        if not active:
            for _ in range(3):
                if all(p.exists() for p in receipts) or any((folder/str(i)/'failure.json').exists() for i in range(len(tasks))):break
                time.sleep(2)
            require(all(p.exists() for p in receipts),'Array terminal without all complete GE receipts');break
        completed=[json.loads(p.read_text()) for p in receipts if p.exists()]
        if completed:
            m.write(a.out/'latest_completed.json',dict(stage=name,cases=completed,deadline_epoch=deadline))
            prior=json.loads((a.out/'best_so_far.json').read_text())['best'];valid=[r for r in completed if r['status']=='passed']
            best=min([prior]+valid,key=lambda r:r['loss']);m.write(a.out/'best_so_far.json',dict(status='provisional_until_two_new_repeats',best=best))
        time.sleep(min(10,max(.1,deadline-time.time())))
    rows=[json.loads(p.read_text()) for p in receipts]
    require(all(r['status'] in ('passed','budget_exhausted','inadmissible_numerical') for r in rows),'Fatal unchanged GE gate failure')
    if name!='gn':require(all(r['status']=='passed' for r in rows),'Required derivative/repeat GE did not pass')
    for row in rows:
        if row['status']!='passed':continue
        report=folder/str(next(i for i,t in enumerate(tasks) if t['label']==row['label']))/row['label']/'phase_b_ge/selected_root';row['report']=str(report)
        require(m.residual(m.readtable(report/'target_fit.csv')).tolist()==row['residual'],'New residual artifact drift')
        m.report_hashes(report)
    return rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--runtime',type=Path,required=True);p.add_argument('--source-root',type=Path,required=True);p.add_argument('--inventory',type=Path,required=True);p.add_argument('--oldrun',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--launcher',type=Path);p.add_argument('--authenticate-only',action='store_true');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    m=load_runtime(a.runtime);started=time.time();deadline=min(started+3600,float(os.environ.get('FAST_COMMON_DEADLINE_EPOCH','inf')));require(math.isfinite(deadline) and deadline>started,'Original absolute deadline expired');submitted=[];new=[]
    try:
        reused,seed,bounds,coordinates,steps,price,hashes=authenticate(m,a.runtime,a.source_root,a.inventory,a.oldrun)
        m.write(a.out/'authentication.json',dict(status='passed',source_fingerprint=PLAN['source_fingerprint'],target_contract_sha256=PLAN['target_fingerprint'],report_hashes=hashes,reused_full_ge=len(reused),original_center=seed,starting_price=price))
        if a.authenticate_only:return
        require(a.launcher is not None,'Launcher required for native arrays')
        candidates=[r for r in reused if r['kind']!='baseline_repeat'];best=min(candidates,key=lambda r:r['loss']);m.write(a.out/'best_so_far.json',dict(status='provisional_until_two_new_repeats',best=best))
        point=dict(seed);k=coordinates[-1];require(k==PLAN['missing_coordinate'],'Missing derivative order drift');h=steps[k]
        if point[k]+h>bounds[k][1]:h=-h
        point[k]+=h
        new+=run_stage(m,a,'fd',[dict(label='010_r0_probe_tenure_choice_kappa',kind='finite_difference',parameters=point,starting_price=price)],deadline,submitted)
        identification,points=jacobian(m,reused[0],reused[2:]+new,coordinates,steps,bounds);m.write(a.out/'identification.json',identification)
        if points:
            tasks=[dict(label=f'011_r0_gn_{d}',kind='damped_Gauss_Newton',parameters=x,starting_price=price) for d,x in zip((.5,.2,1.),points)]
            new+=run_stage(m,a,'gn',tasks,deadline,submitted)
        candidates+=[r for r in new if r['status']=='passed'];selected=min(candidates,key=lambda r:r['loss']);m.write(a.out/'best_so_far.json',dict(status='provisional_until_two_new_repeats',best=selected));m.write(a.out/'selected_provisional.json',selected)
        repeats=run_stage(m,a,'repeat',[dict(label=f'014_selected_repeat_{i}',kind='selected_repeat',parameters=selected['parameters'],starting_price=price) for i in range(2)],deadline,submitted)
        for i,row in enumerate(repeats):
            require(row['residual']==selected['residual'],'Selected repeat residual differs');m.write(a.out/f'selected_repeat_{i}.json',m.compare_repeated(Path(selected['report']),Path(row['report'])))
        new+=repeats;require(len(new)<=6,'New full-GE cap exceeded')
        result=dict(status='selected_verified',lane='floor_s2',selected=selected,baseline=reused[0],selected_repeats=2,completed_new_full_ge=len(new),maximum_new_full_ge=6,identification=identification,elapsed_seconds=time.time()-started,deadline_epoch=deadline,experimental_not_adopted=True)
        m.write(a.out/'cases.json',reused+new);m.write(a.out/'completed.json',result)
    except BaseException as e:
        m.write(a.out/'failure.json',dict(type=type(e).__name__,message=str(e),elapsed_seconds=time.time()-started,no_auto_retry=True));raise
if __name__=='__main__':main()
