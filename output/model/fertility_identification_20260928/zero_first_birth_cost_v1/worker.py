"""Authenticated native evaluation for one new overnight calibration proposal.

All solving, checkpoint reads, hashing and rendering run on Torch. Original
sources, targets and the reference are immutable. The optional two-birth model
uses the separately tested v2 overlay; no equations are implemented here.
"""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key]='1'
import argparse,copy,csv,gzip,hashlib,importlib.util,json,math,pickle,sys,time,traceback,types
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
PRIOR=BASE/'two_stream_overnight_v1'
ROOT=BASE.parents[2]


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''):digest.update(block)
    return digest.hexdigest()


def canon(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def read(path):return json.loads(Path(path).read_text())


def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n');temporary.replace(path)


def require(condition,message):
    if not condition:raise RuntimeError(message)


def verify_config(path):
    require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm required')
    require(not sys.flags.optimize,'Assertions must remain enabled')
    c=read(path)
    require(sha(path)==os.environ.get('EXPECTED_ZERO_COST_CONFIG_SHA256'),'New configuration pin changed')
    require(c['schema']=='zero_first_birth_cost_v1' and not c.get('synthetic',False),'Production configuration required')
    require(set(c['lanes'])=={'zero_cost'},'Single restricted stream required')
    require(len(c['parameters'])==9 and 'first_birth_fixed_cost' not in c['parameters'] and len(set(c['scored_moments']))==10,'Nine free coordinates and ten scored rows required')
    require(c['fixed_parameters']=={'first_birth_fixed_cost':0.0},'Fixed-cost contract changed')
    require(c['budget']==dict(total_seconds=14400,final_reserve_seconds=4200,case_seconds=2100,
        maximum_stationary_solves=8,max_evaluations=16),'Budgets changed')
    for pin in c['pins'].values():
        require(sha(pin['path'])==pin['sha256'],'Changed source/input: '+pin['path'])
    require(Path(c['pins']['worker']['path']).resolve()==Path(__file__).resolve(),'Wrong worker source')
    return c


def driver_module():
    sys.path.insert(0,str(BASE/'two_births_optimized_v2'))
    name='verified_two_birth_experiment_driver'
    spec=importlib.util.spec_from_file_location(name,BASE/'two_births_optimized_v2/driver.py')
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return module


def validate_request(c,req,config_path):
    require(req['config_sha256']==sha(config_path),'Request configuration mismatch')
    lane=c['lanes'][req['lane']]
    for key in ('source_fingerprint','target_fingerprint'):
        require(req[key]==lane[key],'Request '+key+' mismatch')
    require(req['scored_moments']==c['scored_moments'],'Scored row order mismatch')
    require(set(req['point'])==set(c['parameters']),'Wrong proposal coordinates')
    expected_fixed={'first_birth_fixed_cost':c['lanes']['zero_cost']['anchor_fixed_cost']} if req['role']=='anchor_replay' else c['fixed_parameters']
    require(req['fixed_parameters']==expected_fixed,'Fixed parameter request changed')
    for name,value in req['point'].items():
        lo,hi=c['bounds'][name]
        require(math.isfinite(value) and lo<=value<=hi,'Proposal outside original bounds: '+name)
    require(math.isfinite(req['initial_psi']) and req['initial_psi']>0,'Invalid normalization start')
    require(0<req['case_budget_seconds']<=2100 and req['maximum_stationary_solves']==8,'Request budget changed')
    require(time.time()<req['deadline_epoch']<=c['hard_end_epoch'],'Expired request')
    require(req['deadline_epoch']<=time.time()+req['case_budget_seconds']+2,'Request exceeds case budget')


def same_estimate(name,actual,requested):
    if name=='beta_annual':
        require(abs(actual-requested)<=2e-12,'Annual beta differs beyond inherited tolerance')
    else:
        require(actual==requested,'Estimated parameter differs from request: '+name)


def same_accounting_value(actual,expected):
    """Allow only CSV/operator-association roundoff in target accounting."""
    return math.isclose(actual,expected,rel_tol=1e-12,abs_tol=1e-12)


def evaluate(c,req,out,config_path):
    started=time.monotonic();evaluator=None;phase='authentication'
    context=dict(candidate_id=req['candidate_id'],stage='initial',contract_sha256=req['config_sha256'],
        source_sha256=c['pins']['worker']['sha256'],target_sha256=c['pins']['objective']['sha256'],point_sha256=canon(req['point']))
    identities={key:req[key] for key in ('candidate_id','config_sha256','lane','source_fingerprint','target_fingerprint')}
    out.mkdir(parents=True,exist_ok=False)
    write(out/'startup.json',dict(**identities,request=req,pid=os.getpid(),epoch=time.time()))
    try:
        d=driver_module()
        native,obj,evaluator,reference,reference_point,reference_receipt=d.load_runtime(out)
        require(canon(obj['target_rows'])==req['target_fingerprint'],'Complete target/weight fingerprint changed')
        native_bounds={row['parameter']:[row['lower'],row['upper']] for row in obj['parameter_restrictions']}
        require(c['original_bounds']==native_bounds,'Complete original parameter restrictions changed')
        require(c['bounds']=={name:native_bounds[name] for name in c['parameters']},'Parameter restrictions changed')
        require(native_bounds['first_birth_fixed_cost'][0]==0.0,'Original fixed-cost lower bound changed')
        original_bind=evaluator.bind
        requested_cost=req['fixed_parameters']['first_birth_fixed_cost']
        def bind(self,point):
            full=dict(point,first_birth_fixed_cost=requested_cost)
            require(full['first_birth_fixed_cost']==requested_cost,'Fixed cost changed before binding')
            P=original_bind(full)
            require(P.first_birth_fixed_cost==requested_cost,'Bound model cost differs')
            return P
        evaluator.bind=types.MethodType(bind,evaluator)
        evaluator.c=copy.deepcopy(evaluator.c)
        evaluator.c['normalization']['initial_psi']=req['initial_psi']
        evaluator.c['normalization']['maximum_stationary_solves']=8
        evaluator.c['economic_changes']=(['Authenticated replay of selected original one-birth candidate; no economic change.']
            if req['role']=='anchor_replay' else c['lanes'][req['lane']]['economic_changes'])
        phase='native_evaluation'
        case=out/'case'
        full_point=dict(req['point'],first_birth_fixed_cost=requested_cost)
        receipt=evaluator.evaluate(full_point,case,deadline_epoch=req['deadline_epoch'],graphs=True)
        require(receipt['point']['first_birth_fixed_cost']==requested_cost,'Native receipt changed fixed cost')
        ledger=read(case/'stationary_solves.json')
        require(1<=len(ledger)<=8 and all(row['status']=='completed' for row in ledger),'Stationary solve ledger incomplete')
        require(receipt['objective_stationary_solves']==len(ledger),'Stationary solve count mismatch')
        require(math.isfinite(receipt['normalization']['psi_child']) and receipt['normalization']['psi_child']>0,'Invalid normalized child benefit')
        require(abs(receipt['normalization']['completed_fertility']-2.1)<=5e-4,'Completed-fertility normalization failed')
        require(abs(receipt['adult_entry_gate']['fertility_gap'])<=5e-4,'Demographic renewal failed')
        phase='export_validation'
        with (case/'target_fit.csv').open(newline='') as stream:fit=list(csv.DictReader(stream))
        with (case/'parameters.csv').open(newline='') as stream:parameters=list(csv.DictReader(stream))
        require(len(fit)==14 and len(parameters)==31,'Complete tables missing')
        with Path(c['pins']['reference_parameters.csv']['path']).open(newline='') as stream:
            expected_parameter_names={row['parameter'] for row in csv.DictReader(stream)}
        actual_parameter_names=[row['parameter'] for row in parameters]
        require(len(set(actual_parameter_names))==31 and set(actual_parameter_names)==expected_parameter_names,
                'Parameter export names differ from frozen 31-name reference')
        require(all(math.isfinite(float(row['estimate'])) for row in parameters),'Nonfinite parameter estimate')
        require(set(req['point']).issubset(actual_parameter_names),'Estimated coordinates missing from export')
        fixed_row=next(row for row in parameters if row['parameter']=='first_birth_fixed_cost')
        require(float(fixed_row['estimate'])==requested_cost,'Saved parameter cost differs')
        native_rows={row['restriction_id']:row for row in obj['target_rows']}
        require(len({row['moment'] for row in fit})==14 and set(row['moment'] for row in fit)==set(native_rows),
                'Target-fit names differ from frozen fourteen-row target system')
        for row in fit:
            target=native_rows[row['moment']]
            target_value=model_value=gap=float('nan')
            target_value=float(row['target']);model_value=float(row['model']);gap=float(row['gap'])
            require(all(math.isfinite(x) for x in (target_value,model_value,gap)),'Nonfinite target/model/gap export')
            require(same_accounting_value(gap,model_value-target_value),'Target-fit gap is not model minus target')
            require(target_value==float(target['target']),'Target value changed')
            expected=target['actual_weight']
            require((row['weight']=='' if expected is None else float(row['weight'])==expected),'Weight changed')
            role='normalization' if expected is None else 'validation' if expected==0 else 'scored'
            require(row['role']==role,'Target role changed')
            if expected is None:
                require(row['loss_contribution']=='','Normalization row must not contribute to loss')
            else:
                contribution=float(row['loss_contribution'])
                require(math.isfinite(contribution) and same_accounting_value(contribution,expected*gap*gap),
                        'Loss contribution accounting changed')
        indexed={row['moment']:row for row in fit}
        residuals=[math.sqrt(float(indexed[name]['weight']))*float(indexed[name]['gap']) for name in c['scored_moments']]
        require(all(math.isfinite(value) for value in residuals),'Nonfinite scored residual')
        loss=sum(value*value for value in residuals)
        require(abs(loss-receipt['loss'])<=1e-8*max(1,loss),'Loss accounting changed')
        for row in parameters:
            name=row['parameter']
            if name in req['point']:
                same_estimate(name,float(row['estimate']),req['point'][name])
                require([float(row['lower']),float(row['upper'])]==c['bounds'][name],'Reported bounds changed')
                row['status']='estimated candidate in '+req['lane']+' calibration stream; not adopted'
            elif name=='first_birth_fixed_cost':
                row['status']='fixed at authenticated anchor cost for replay; not adopted' if req['role']=='anchor_replay' else 'fixed exactly zero in experimental restricted diagnostic; not adopted'
        with (case/'parameters.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(parameters[0]),lineterminator='\n');writer.writeheader();writer.writerows(parameters)
        plot_names=sorted(path.name for path in (case/'standard_diagnostics').glob('*.png'))
        require(plot_names==sorted(native['standard_diagnostic_names']),'Standard 17 plots changed')
        require(len(plot_names)==17,'Standard diagnostic packet incomplete')
        receipt.update(**identities,status='verified_zero_cost_diagnostic_not_adopted',free_count=9,
            estimated_coordinates=9,normalized_coordinates=1,fixed_parameters=req['fixed_parameters'],scientific_promotion=False,
            reference_label=d.LABEL,reference_checkpoint_sha256=reference_receipt['case_checkpoint_sha256'],
            worker_elapsed_seconds=time.monotonic()-started,
            economic_changes=evaluator.c['economic_changes'],
            fertility_plot_scope='Original one-birth model')
        write(case/'receipt.json',receipt)
        write(case/'scientific_identity.json',dict(**identities,
            birth_rule='original_one_birth',point_sha256=canon(req['point']),fixed_parameters=req['fixed_parameters'],
            normalized_psi=receipt['normalization']['psi_child'],
            checkpoint_sha256=receipt['case_checkpoint_sha256'],
            effective_source_manifest_sha256=receipt.get('effective_source_manifest_sha256'),
            scientific_promotion=False))
        artifacts=[]
        for path in sorted(case.rglob('*')):
            if path.is_file() and not path.name.endswith('.pkl.gz'):
                artifacts.append(dict(path=str(path.resolve()),sha256=sha(path)))
        success=dict(**identities,status='passed',loss=loss,residuals=residuals,point=req['point'],
            fixed_parameters=req['fixed_parameters'],
            psi=receipt['normalization']['psi_child'],case_path=str(case.resolve()),
            model_evaluations=receipt['objective_stationary_solves'],elapsed_seconds=time.monotonic()-started,
            checkpoint_sha256=receipt['case_checkpoint_sha256'],artifacts=artifacts,
            complete_target_rows=14,complete_parameter_rows=31,standard_plots=17,
            scientific_promotion=False)
        write(out/'SUCCESS.json',success)
    except Exception as exc:
        failure=dict(**identities,status='fatal',classification='unknown_or_integrity_failure',
            phase=phase,error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
        ledger_path=out/'case/stationary_solves.json'
        ledger=read(ledger_path) if ledger_path.exists() else []
        if (phase=='native_evaluation' and type(exc) is RuntimeError and
            str(exc)=='Old-steady-state fertility normalization missed tolerance: 23-solve cap' and
            len(ledger)==8 and all(row['status']=='completed' for row in ledger)):
            failure.update(status='censored',classification='authenticated_eight_solve_resource_cap',authenticated=True,model_evaluations=8)
        elif phase=='native_evaluation' and evaluator is not None:
            try:failure.update(evaluator.classify_failure(exc,context))
            except Exception as error:failure['classifier_error']=str(error)
        write(out/'FAILURE.json',failure)
        raise


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True)
    p.add_argument('--request',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();c=verify_config(a.config);req=read(a.request);validate_request(c,req,a.config)
    evaluate(c,req,a.output,a.config)


if __name__=='__main__':main()
