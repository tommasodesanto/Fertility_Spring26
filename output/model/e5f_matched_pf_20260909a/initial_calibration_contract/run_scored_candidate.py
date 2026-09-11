"""Pinned initial-solve -> original gates -> checkpoint observer -> external score.

Run schema e5f_scored_initial_candidate_run_v1 requires source_root, case_id,
seconds (<=2100), wrapper_sha256; initial_solve_contract, working_objective,
scorer and validator each {path,sha256}. working_objective additionally requires
canonical_sha256. objective_source_files maps every objective source key to
{path,hash_kind}, where hash_kind is bytes or canonical_json. All paths are
explicit or relative to the run-contract folder. Original solve contract stays
unchanged: new_balanced, normalize/observe_early True, eight solves per repetition,
one or two repetitions, <=1800 seconds. Two repetitions are the exact-loop smoke.
Raw outputs are never rewritten. A zero exit certifies this candidate's score,
not calibration optimality, identification or policy readiness.
"""
from __future__ import annotations
import argparse
import copy
import csv
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import signal
import subprocess
import sys
import time

SCHEMA = 'e5f_scored_initial_candidate_run_v1'
DRIVER = 'code/model/tools/run_e5f_initial_revision_probe.py'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path=Path(path); temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
    temporary.replace(path)


def verify(path, pin):
    require(digest(path)==pin, f'Fingerprint mismatch: {path}')


def resolve(base, path):
    return (base/Path(path)).resolve()


def load_module(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def preflight(contract_path, expected_sha256):
    """No model import or subprocess may precede complete independent pin checks."""
    contract_path=Path(contract_path).resolve(); base=contract_path.parent
    verify(contract_path,expected_sha256); c=read(contract_path)
    require(c.get('schema')==SCHEMA,'Wrong scored-run schema')
    require(isinstance(c.get('case_id'),str) and c['case_id'].strip(),'Missing candidate id')
    require(type(c.get('seconds')) is int and 1<=c['seconds']<=2100,'Invalid wrapper budget')
    require(not sys.flags.optimize,'Original assert-based validator requires Python optimization disabled')
    verify(Path(__file__),c['wrapper_sha256'])
    files={key:resolve(base,c[key]['path']) for key in
           ('initial_solve_contract','working_objective','scorer','validator')}
    for key,path in files.items(): verify(path,c[key]['sha256'])
    scorer=load_module('pinned_initial_scorer',files['scorer'])
    objective=read(files['working_objective'])
    require(scorer.fingerprint(objective)==c['working_objective']['canonical_sha256'],
            'Working objective canonical fingerprint mismatch')
    declarations=c['objective_source_files']
    require(set(declarations)==set(objective['source_fingerprints']),'Objective source inventory mismatch')
    for key,entry in declarations.items():
        path=resolve(base,entry['path']); kind=entry['hash_kind']
        require(kind in ('bytes','canonical_json'),'Unsupported source fingerprint kind')
        actual=digest(path) if kind=='bytes' else scorer.fingerprint(read(path))
        require(actual==objective['source_fingerprints'][key],f'Objective source fingerprint mismatch: {key}')
    initial=read(files['initial_solve_contract']); root=resolve(base,c['source_root'])
    require(initial.get('schema')=='e5f_parenthood_initial_probe_v1','Wrong original solve schema')
    require(initial['normalize'] is True and initial['observe_early'] is True,'Complete early normalized solve required')
    require(initial['repetitions'] in (1,2) and initial['maximum_stationary_solves_per_repetition']==8,
            'Preserve original eight-solve repetition gate')
    require(1<=initial['seconds']<=1800 and c['seconds']>=initial['seconds'], 'Original solve/wrapper budget mismatch')
    require(initial['calibrated_smm'] is False and initial['payroll_tax']==.179
            and initial['housing_supply_elasticity']==.63 and initial['fertility_normalization']==2.1,
            'Original scientific contract changed')
    require(set(initial['structural_candidate'])==set(scorer.PARAMETERS),'Complete nine-coordinate candidate required')
    sources=initial['source_sha256']
    inventory={str(p.relative_to(root)) for p in (root/'code/model').rglob('*.py')}
    require(inventory and inventory.issubset(sources),'Original model source inventory incomplete')
    for relative,pin in sources.items():
        path=resolve(root,relative)
        require(path.is_relative_to(root),'Original source escapes snapshot')
        verify(path,pin)
    verify(Path(initial['normalized_checkpoint']),initial['normalized_checkpoint_sha256'])
    require(DRIVER in sources,'Original driver omitted from source pins')
    require('code/model/tools/e5f_recent_parent_flow_observer.py' in sources,'Recent-parent observer omitted from source pins')
    return c,initial,objective,scorer,load_module('pinned_original_validator',files['validator']),root,files


def run_child(command, *, cwd, env, log, seconds):
    require(seconds>0,'Scored-candidate time budget exhausted')
    with Path(log).open('w') as stream:
        child=subprocess.Popen(command,cwd=cwd,env=env,stdout=stream,stderr=subprocess.STDOUT,
                               start_new_session=True)
        try:
            code=child.wait(timeout=seconds)
        except BaseException:
            try: os.killpg(child.pid,signal.SIGKILL)
            except ProcessLookupError: pass
            child.wait()
            raise
    require(code==0,f'Child failed with exit {code}; see {log}')


def csv_rows(path):
    with Path(path).open() as stream: return list(csv.DictReader(stream))


def write_table(path, rows):
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        writer.writerows({k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) else v
                          for k,v in row.items()} for row in rows)


def validate_raw(raw, initial, initial_sha, validator):
    """Use original numeric validators; do not substitute fixed historical case ids."""
    require(not (raw/'failure.json').exists() and not (raw/'timeout.json').exists(),'Original solve failed/timed out')
    require(read(raw/'contract.json')==dict(initial,contract_sha256=initial_sha,case='new_balanced'),
            'Original output contract differs')
    top=read(raw/'summary.json')
    require(top['case']=='new_balanced' and top['repetitions']==initial['repetitions'],'Original repetition identity differs')
    require(top['elapsed_seconds']<=initial['seconds'],'Original time gate failed')
    reps=[];total=0
    for number in range(1,initial['repetitions']+1):
        folder=raw/f'repetition_{number:02d}'; summary=read(folder/'summary.json')
        early=read(folder/'early_measurement.json');parameters=csv_rows(folder/'parameters.csv')
        ges=read(folder/'stationary_solves.json')
        require(len(parameters)==17,'Full seventeen-parameter table required')
        view=dict(top,repetitions=1,stationary_solves=len(ges),final=summary)
        validator.validate_summary(view,initial,early,parameters,ges)
        verify(folder/'initial_state.pkl.gz',summary['checkpoint_sha256'])
        total+=len(ges);reps.append(dict(folder=folder,summary=summary,early=early,parameters=parameters))
    require(total==top['stationary_solves']<=8*initial['repetitions'],'Original aggregate solve budget failed')
    require(top['final']==reps[-1]['summary'],'Original final summary differs')
    validator.validate_market(read(reps[-1]['folder']/'market_quantity_units.json'))
    graphs=sorted((reps[-1]['folder']/'standard_diagnostics').glob('*.png'))
    require(len(graphs)==17,'Original standard seventeen-graph packet incomplete')
    if len(reps)==2:
        require(reps[0]['early']==reps[1]['early'],'Original early observation replay differs')
        require(reps[0]['parameters']==reps[1]['parameters'],'Original parameter replay differs')
        for key in ('price','legacy_stationary_moments'):
            require(reps[0]['summary'][key]==reps[1]['summary'][key],'Original replay differs: '+key)
        norm=[{k:v for k,v in rep['summary']['normalization'].items() if k!='stationary_solve_seconds'} for rep in reps]
        require(norm[0]==norm[1],'Original normalization replay differs')
    return reps,[dict(path=str(p.relative_to(raw)),sha256=digest(p)) for p in graphs]


def observe_checkpoint(checkpoint, pin, root, output, case_id):
    """Internal bounded child, reached only after parent validates original gates."""
    verify(checkpoint,pin)
    sys.path[:0]=[str(root/'code/model/tools'),str(root/'code/model')]
    import run_e5f_open_population_transition as transition
    from intergen_eqscale_seq_optimized import solver as model
    from e5f_recent_parent_flow_observer import observe_recent_parent_flow,SNAPSHOT,AGE_PROJECTION
    transition.calendar.model=model
    require(os.environ.get('NUMBA_DISABLE_JIT')=='0','Compiled observation required')
    with gzip.open(checkpoint,'rb') as stream: packet=pickle.load(stream)
    require(packet['parameters'].use_numba_scatter is True,'Checkpoint must request compiled scatter')
    result=observe_recent_parent_flow(packet['evaluation'],packet['parameters'],diagnostic_enabled=True,
        snapshot=SNAPSHOT,age_projection=AGE_PROJECTION,diagnostic_allow_residence_proxy=True,
        input_provenance=dict(case_id=case_id,checkpoint_sha256=pin))
    write(output,result)


def run(contract, contract_sha256, output):
    started=time.monotonic()
    c,initial,objective,scorer,validator,root,files=preflight(contract,contract_sha256)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,NUMBA_DISABLE_JIT='0',PYTHONOPTIMIZE='0',
             PYTHONPATH=f'{root}/code/model/tools:{root}/code/model',
             OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
    def remaining(): return c['seconds']-(time.monotonic()-started)
    try:
        write(out/'preflight.json',dict(status='verified',run_contract_sha256=contract_sha256,
            objective_canonical_sha256=c['working_objective']['canonical_sha256'],
            objective_file_sha256=c['working_objective']['sha256'],
            initial_solve_contract_sha256=c['initial_solve_contract']['sha256']))
        run_child([sys.executable,str(root/DRIVER),'--contract',str(files['initial_solve_contract']),
            '--contract-sha256',c['initial_solve_contract']['sha256'],'--case','new_balanced',
            '--output',str(out/'raw')],cwd=root,env=env,log=out/'initial_solve.log',
            seconds=min(initial['seconds'],remaining()))
        reps,graphs=validate_raw(out/'raw',initial,c['initial_solve_contract']['sha256'],validator)
        scores=[];recent_packets=[];receipts=[]
        for number,rep in enumerate(reps,1):
            dest=out/f'scored_repetition_{number:02d}';dest.mkdir()
            checkpoint=rep['folder']/'initial_state.pkl.gz';pin=rep['summary']['checkpoint_sha256']
            run_child([sys.executable,str(Path(__file__).resolve()),'observe','--checkpoint',str(checkpoint),
                '--checkpoint-sha256',pin,'--source-root',str(root),'--output',str(dest/'recent_parent.json'),
                '--case-id',c['case_id']],cwd=root,env=env,log=dest/'observer.log',seconds=remaining())
            recent=read(dest/'recent_parent.json')
            inputs=dict(early_measurement=rep['early'],recent_parent_observation=recent,
                        normalization=rep['summary']['normalization'],parameters=rep['parameters'])
            # Certification denotes the explicitly maintained approximation;
            # the raw observer's diagnostic/exact-ACS flags are not rewritten.
            receipt=dict(schema='e5f_initial_score_receipt_v1',status='verified',
                source_fingerprints=objective['source_fingerprints'],checkpoint_sha256=pin,
                numerical_gates_verified=True,recent_parent_certified=True,
                recent_parent_approximation_id=objective['recent_parent_approximation']['approximation_id'],
                input_sha256={k:scorer.fingerprint(v) for k,v in inputs.items()})
            scored=scorer.score_initial(objective,expected_contract_sha256=c['working_objective']['canonical_sha256'],
                                       evaluation_receipt=receipt,**inputs)
            require(len(scored['target_fit'])==13 and len(scored['parameters'])==17,'Incomplete score readout')
            write(dest/'verified_evaluation_receipt.json',receipt)
            write(dest/'score.json',scored);write_table(dest/'target_fit.csv',scored['target_fit'])
            write_table(dest/'parameters.csv',scored['parameters'])
            scores.append(scored);recent_packets.append(recent);receipts.append(scorer.fingerprint(receipt))
        if len(scores)==2:
            require(scores[0]['loss']==scores[1]['loss'],'Exact repeated loss differs')
            projections=[]
            for recent in recent_packets:
                packet=copy.deepcopy(recent)
                packet['metadata']['policy_input_provenance'].pop('checkpoint_sha256')
                projections.append(packet)
            require(projections[0]==projections[1],'Recent-parent numerical replay differs')
        result=dict(status='verified_scored_candidate',case_id=c['case_id'],repetitions=len(scores),
            loss=scores[0]['loss'],exact_loss_equality=(len(scores)==2),
            objective_canonical_sha256=c['working_objective']['canonical_sha256'],
            objective_file_sha256=c['working_objective']['sha256'],initial_solve_contract_sha256=c['initial_solve_contract']['sha256'],
            run_contract_sha256=contract_sha256,evaluation_receipt_sha256=receipts,original_graphs=graphs,
            calibrated_smm=False,benchmark_certified=False,elapsed_seconds=time.monotonic()-started)
        require(result['elapsed_seconds']<=c['seconds'],'Wrapper deadline exceeded')
        write(out/'summary.json',result);return result
    except BaseException as exc:
        write(out/'failure.json',dict(status='failed_scored_candidate',error_type=type(exc).__name__,
            error=str(exc),elapsed_seconds=time.monotonic()-started))
        raise


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='mode',required=True)
    runner=sub.add_parser('run');runner.add_argument('--contract',type=Path,required=True)
    runner.add_argument('--contract-sha256',required=True);runner.add_argument('--output',type=Path,required=True)
    observer=sub.add_parser('observe');observer.add_argument('--checkpoint',type=Path,required=True)
    observer.add_argument('--checkpoint-sha256',required=True);observer.add_argument('--source-root',type=Path,required=True)
    observer.add_argument('--output',type=Path,required=True);observer.add_argument('--case-id',required=True)
    args=parser.parse_args()
    if args.mode=='run': print(json.dumps(run(args.contract,args.contract_sha256,args.output)))
    else: observe_checkpoint(args.checkpoint,args.checkpoint_sha256,args.source_root,args.output,args.case_id)


if __name__=='__main__':main()
