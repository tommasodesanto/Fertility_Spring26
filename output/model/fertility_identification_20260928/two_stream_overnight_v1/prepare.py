"""Pin the new two-stream experiment on Torch; never dispatch or resume a search."""
from __future__ import annotations
import csv,datetime,json,os,sys,time
from pathlib import Path
from worker import HERE,BASE,ROOT,sha,canon,read,write,require,driver_module


def same_estimate(name,actual,receipt_value):
    if name=='beta_annual':
        require(abs(actual-receipt_value)<=2e-12,'Annual beta differs beyond inherited tolerance')
    else:
        require(actual==receipt_value,'Anchor table differs from receipt point: '+name)


def main():
    require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm required')
    require(not (HERE/'config.json').exists(),'Preparation is exclusive-create')
    d=driver_module()
    native,obj,evaluator,reference,reference_point,reference_receipt=d.load_runtime(HERE/'preparation')
    pins={}
    def pin(name,path):
        path=Path(path).resolve();pins[name]=dict(path=str(path),sha256=sha(path));return pins[name]
    for name in ('worker','search','tests_search','prepare'):
        pin(name,HERE/(name+'.py'))
    pin('run',HERE/'run.sh')
    pin('integration_smoke',HERE/'integration_smoke.py')
    pin('original_contract',BASE/'contract_v1/contract.json')
    pin('objective',BASE/'contract_v1/primary_objective.json')
    pin('original_source_manifest',native['source_manifest']['path'])
    for name in ('receipt.json','target_fit.csv','parameters.csv'):
        pin('reference_'+name,BASE/'resume_v1/selected_export/primary'/name)
    two=BASE/'two_births_optimized_v2'
    smoke=read(two/'smoke_v1/worker/receipt.json')
    require(smoke['status']=='passed' and smoke['unit_tests']==16,'Two-birth regression smoke failed')
    for name,digest in smoke['source_pins'].items():
        require(sha(two/name)==digest,'Verified two-birth source changed')
        pin('two_birth_'+name,two/name)
    require(read(two/'run_v1/worker/receipt.json')['status']=='passed','Two-birth anchor failed')
    pin('two_birth_smoke',two/'smoke_v1/worker/receipt.json')
    pin('two_birth_passed_run',two/'run_v1/worker/receipt.json')
    pair=BASE/'numerical_pair_v1/run_v1'
    require(read(pair/'paired_comparison.json')['screens']['passed'],'One-birth paired verification failed')
    pin('one_birth_pair',pair/'paired_comparison.json')
    parameters=list(reference_point)
    bounds={row['parameter']:[row['lower'],row['upper']] for row in obj['parameter_restrictions']}
    scored=[row['restriction_id'] for row in obj['target_rows'] if row['actual_weight'] is not None and row['actual_weight']>0]
    require(len(parameters)==len(scored)==10,'Identification count changed')
    fingerprint=canon(obj['target_rows'])
    require(fingerprint==native['lanes']['primary']['target_weight_fingerprint'],'Target fingerprint mismatch')
    lanes={}
    for lane,case,seed in [('one_birth',pair/'retained_start/case',2609281),
                           ('two_birth',two/'run_v1/worker/case',2609282)]:
        receipt=read(case/'receipt.json')
        require(sha(case/'initial_state.pkl.gz')==receipt['case_checkpoint_sha256'],'Anchor checkpoint changed')
        with (case/'parameters.csv').open(newline='') as stream:rows={r['parameter']:r for r in csv.DictReader(stream)}
        point={name:float(receipt['point'][name]) for name in parameters}
        require(set(receipt['point'])==set(parameters),'Anchor receipt point coordinate set changed')
        for name in parameters:
            same_estimate(name,float(rows[name]['estimate']),point[name])
        require(abs(receipt['normalization']['completed_fertility']-2.1)<=5e-4,'Anchor misses normalization')
        require(abs(receipt['adult_entry_gate']['fertility_gap'])<=5e-4,'Anchor misses renewal')
        require(receipt['target_system_sha256']==pins['objective']['sha256'],'Anchor objective differs')
        for name in ('receipt.json','target_fit.csv','parameters.csv'):
            pin(lane+'_anchor_'+name,case/name)
        identity=dict(original_source_manifest=native['source_manifest']['sha256'],birth_rule=lane,
            two_birth_sources=smoke['source_pins'] if lane=='two_birth' else {},
            target_weight_fingerprint=fingerprint,bounds=bounds)
        changes=[
            'Experimental recalibration: search the same ten coordinates within original bounds.',
            'Child-benefit parameter adjusted to completed-fertility target2.1 and unchanged demographic renewal gate.',
            'Earnings, initial wealth/income distributions, transfers/floors, housing/credit primitives, targets and weights unchanged.'
        ]
        if lane=='two_birth':changes += [
            'Experimental: optional second birth attempt after a successful birth, at most two within each four-year cell.',
            'Experimental: existing later-birth Gumbel scale/inclusive value and independent second conception draw.',
            'Retained common-event linear within-period age projection; separately dated births not modeled.'
        ]
        lanes[lane]=dict(initial_point=point,initial_psi=float(rows['psi_child']['estimate']),
            source_fingerprint=canon(identity),source_identity=identity,target_fingerprint=fingerprint,
            seed=seed,anchor_case=str(case),anchor_checkpoint_sha256=receipt['case_checkpoint_sha256'],
            anchor_loss=receipt['loss'],economic_changes=changes,scientific_promotion=False)
    end=datetime.datetime(2026,9,29,11,30,tzinfo=datetime.timezone.utc).timestamp()
    require(end-time.time()>7*3600,'Insufficient overnight window; prepare a new explicit budget')
    config=dict(schema='two_stream_search_v1',parameters=parameters,scored_moments=scored,bounds=bounds,
        lanes=lanes,worker_command=[sys.executable,str(HERE/'worker.py'),'--config',str(HERE/'config.json')],
        hard_end_epoch=end,budget=dict(total_seconds=25200,final_reserve_seconds=4200,case_seconds=2100,
            maximum_stationary_solves=8,max_evaluations=36),pins=pins,
        reference_label=d.LABEL,reference_checkpoint_sha256=reference_receipt['case_checkpoint_sha256'],
        identification='Ten scored rows for ten searched coordinates plus one separately normalized coordinate; rank not yet certified.',
        authorization='September28 author request: split overnight calibration between original and experimental two-birth versions.',
        prepared_epoch=time.time(),scientific_promotion=False,transition_run=False)
    write(HERE/'config.json',config)
    write(HERE/'preparation_receipt.json',dict(status='prepared_not_launched',config_sha256=sha(HERE/'config.json'),
        source_pins=pins,parameters=parameters,scored_moments=scored,
        worker_count=2,max_objectives_total=72,max_stationary_solves_total=576,
        actual_elapsed_time_takes_precedence_over_case_ceilings=True,
        hard_end_epoch=end,model_solves=0,scientific_promotion=False))
    print(json.dumps(dict(status='prepared_not_launched',config_sha256=sha(HERE/'config.json'))))


if __name__=='__main__':main()
