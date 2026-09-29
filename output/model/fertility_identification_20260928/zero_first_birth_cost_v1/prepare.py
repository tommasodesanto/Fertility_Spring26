"""Prepare an isolated, source-pinned restricted diagnostic on Torch; no model solve."""
from __future__ import annotations
import csv, datetime, json, os, sys, time
from pathlib import Path
from worker import HERE, BASE, PRIOR, sha, canon, read, write, require

REMOTE=Path('/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project')
SELECTED=REMOTE/'output/model/fertility_identification_20260928/two_stream_overnight_v1/run_v1/one_birth/one_birth_024_gn1_0/case'

def main():
    require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm required')
    require(not (HERE/'config.json').exists(),'Preparation is exclusive-create')
    old=read(PRIOR/'config.json')
    require(sha(PRIOR/'config.json')=='419b46d7cffdf51d95760f2999aec67a66b10b08bb9198164ac3d6f7d96323d3','Original prepared configuration changed')
    require(old['schema']=='two_stream_search_v1','Original search contract changed')
    receipt=read(SELECTED/'receipt.json')
    require(receipt['status']=='verified_overnight_candidate_not_adopted','Selected receipt status changed')
    require(receipt['candidate_id']=='one_birth_024_gn1_0' and receipt['lane']=='one_birth','Wrong selected case')
    require(abs(receipt['loss']-7.826226594410982)<1e-12,'Selected loss changed')
    require(receipt['point']['first_birth_fixed_cost']==0.35270914196085973,'Anchor first-birth cost changed')
    require(receipt['normalization']['psi_child']==0.12184551693359474,'Anchor normalized benefit changed')
    require(receipt['target_system_sha256']==old['pins']['objective']['sha256'],'Selected target contract changed')
    require(receipt['case_checkpoint_sha256']==read(SELECTED/'scientific_identity.json')['checkpoint_sha256'],'Selected checkpoint identity disagrees')
    require(receipt['case_checkpoint_sha256']=='308245b919c8a50d25972b29e55b5d4fca1ae57d01155ff6050724dc61610451','Selected checkpoint identity changed')
    with (SELECTED/'parameters.csv').open(newline='') as f:rows={r['parameter']:r for r in csv.DictReader(f)}
    require(len(rows)==31 and float(rows['first_birth_fixed_cost']['estimate'])==receipt['point']['first_birth_fixed_cost'],'Anchor parameter table differs')
    names=[n for n in old['parameters'] if n!='first_birth_fixed_cost']
    require(len(names)==9 and len(old['scored_moments'])==10,'Original scored/free count changed')
    point={n:receipt['point'][n] for n in names}
    for name in names:
        require(abs(float(rows[name]['estimate'])-point[name])<=(2e-12 if name=='beta_annual' else 0),'Anchor point/table disagreement: '+name)
    pins={}
    for name in ('worker','search','tests_search','prepare'):
        path=HERE/(name+'.py');pins[name]=dict(path=str(path),sha256=sha(path))
    path=HERE/'run.sh';pins['run']=dict(path=str(path),sha256=sha(path))
    tests=read(HERE/'TESTS.json')
    require(tests['status']=='passed' and tests['test_count']==5 and tests['synthetic'] is True and
        tests['model_imports']==tests['model_solves']==0,'Exact-loop tests did not pass')
    tested_sources={name+'.py':pins[name]['sha256'] for name in ('worker','search','tests_search','prepare')}
    tested_sources['run.sh']=pins['run']['sha256']
    require(tests['source_hashes']==tested_sources,'Tested sources changed')
    pins['tests_receipt']=dict(path=str(HERE/'TESTS.json'),sha256=sha(HERE/'TESTS.json'))
    for name in ('objective','original_contract','original_source_manifest','reference_parameters.csv','reference_receipt.json','reference_target_fit.csv','two_birth_driver.py'):
        pins[name]=old['pins'][name]
    for name in ('receipt.json','parameters.csv','target_fit.csv','scientific_identity.json'):
        path=SELECTED/name;pins['selected_'+name]=dict(path=str(path),sha256=sha(path))
    for name,pin in pins.items():
        require(sha(pin['path'])==pin['sha256'],'Pinned small source/input changed: '+name)
    source_identity=dict(original_source_manifest=old['pins']['original_source_manifest']['sha256'],
        evaluator_source=old['pins']['two_birth_driver.py']['sha256'],birth_rule='original_one_birth',
        selected_checkpoint_sha256=receipt['case_checkpoint_sha256'],fixed_first_birth_cost=0.0,
        target_weight_fingerprint=old['lanes']['one_birth']['target_fingerprint'])
    lane=dict(anchor_point=point,initial_point=point,anchor_fixed_cost=receipt['point']['first_birth_fixed_cost'],
        initial_psi=receipt['normalization']['psi_child'],anchor_case=str(SELECTED),anchor_loss=receipt['loss'],
        anchor_checkpoint_sha256=receipt['case_checkpoint_sha256'],source_identity=source_identity,
        source_fingerprint=canon(source_identity),target_fingerprint=old['lanes']['one_birth']['target_fingerprint'],
        economic_changes=['Experimental: first_birth_fixed_cost fixed exactly zero; original one-birth opportunity retained.',
            'Nine other original search coordinates retain original bounds; psi_child renormalized to completed fertility 2.1.',
            'Original ten scored targets, three validation rows, demographic renewal gate, earnings, entry distributions, transfers/floors, preferences except fixed cost, and housing/credit primitives retained.'],
        scientific_promotion=False)
    end=time.time()+14400
    config=dict(schema='zero_first_birth_cost_v1',parameters=names,scored_moments=old['scored_moments'],
        bounds={name:old['bounds'][name] for name in names},original_bounds=old['bounds'],
        fixed_parameters={'first_birth_fixed_cost':0.0},lanes={'zero_cost':lane},
        worker_command=[sys.executable,str(HERE/'worker.py'),'--config',str(HERE/'config.json')],
        hard_end_epoch=end,budget=dict(total_seconds=14400,final_reserve_seconds=4200,case_seconds=2100,
            maximum_stationary_solves=8,max_evaluations=16),pins=pins,
        reference_label=old['reference_label'],reference_checkpoint_sha256=old['reference_checkpoint_sha256'],
        identification='Ten scored rows for nine free coordinates; psi separately normalized. Rank unverified.',
        authorization='Author-authorized isolated diagnostic; launch only after lead review.',
        prepared_epoch=time.time(),scientific_promotion=False,transition_run=False)
    write(HERE/'config.json',config)
    write(HERE/'preparation_receipt.json',dict(status='prepared_not_launched',config_sha256=sha(HERE/'config.json'),
        max_evaluations=16,max_stationary_solves=128,anchor_checkpoint_sha256=receipt['case_checkpoint_sha256'],
        model_solves=0,scientific_promotion=False))

if __name__=='__main__':main()
