"""Prepare two isolated exact-zero lanes on Torch; no model solve."""
from __future__ import annotations
import csv, os, sys, time
from pathlib import Path
from worker import HERE, BASE, PRIOR, sha, canon, read, write, require

REMOTE=Path('/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project')
SELECTED=REMOTE/'output/model/fertility_identification_20260928/two_stream_overnight_v1/run_v1/one_birth/one_birth_024_gn1_0/case'
FIXED={'first_scale_zero':'kappa_fert','later_scale_zero':'kappa_fert_continuation'}

def main():
    require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm required')
    require(not (HERE/'config.json').exists(),'Exclusive preparation required')
    old=read(PRIOR/'config.json')
    require(sha(PRIOR/'config.json')=='419b46d7cffdf51d95760f2999aec67a66b10b08bb9198164ac3d6f7d96323d3','Original configuration changed')
    receipt=read(SELECTED/'receipt.json')
    require(receipt['candidate_id']=='one_birth_024_gn1_0' and receipt['lane']=='one_birth','Wrong E01 case')
    require(abs(receipt['loss']-7.826226594410982)<1e-12,'E01 loss changed')
    require(receipt['point']['first_birth_fixed_cost']==0.35270914196085973,'E01 positive cost changed')
    require(receipt['normalization']['psi_child']==0.12184551693359474,'E01 psi changed')
    require(receipt['case_checkpoint_sha256']==read(SELECTED/'scientific_identity.json')['checkpoint_sha256'],'Checkpoint identity mismatch')
    require(len(old['scored_moments'])==10 and len(old['parameters'])==10,'Target/parameter count changed')
    with (SELECTED/'parameters.csv').open(newline='') as f: rows={r['parameter']:r for r in csv.DictReader(f)}
    require(len(rows)==31,'Full E01 parameters unavailable')
    point={name:receipt['point'][name] for name in old['parameters']}
    for name,value in point.items():
        require(abs(float(rows[name]['estimate'])-value)<=(2e-12 if name=='beta_annual' else 0), 'E01 parameter mismatch: '+name)
    tests=read(HERE/'TESTS.json')
    require(tests['status']=='passed' and tests['synthetic'] is True and tests['model_solves']==0,'Torch tests required')
    pins={}
    for name in ('worker.py','search.py','tests_search.py','tests_solver.py','prepare.py','run.sh','patch_solver.py','overlay_runtime.py'):
        path=HERE/name;pins[name]=dict(path=str(path),sha256=sha(path))
        require(tests['source_hashes'][name]==pins[name]['sha256'],'Tested source changed: '+name)
    path=HERE/'TESTS.json';pins['tests_receipt']=dict(path=str(path),sha256=sha(path))
    for name in ('objective','original_contract','original_source_manifest','reference_parameters.csv','reference_receipt.json','reference_target_fit.csv','two_birth_driver.py'):
        pins[name]=old['pins'][name]
    for name in ('receipt.json','parameters.csv','target_fit.csv','scientific_identity.json'):
        path=SELECTED/name;pins['selected_'+name]=dict(path=str(path),sha256=sha(path))
    solver=BASE.parents[2]/'code/model/intergen_eqscale_seq_optimized/solver.py'
    pins['solver']=dict(path=str(solver),sha256=sha(solver))
    for name,pin in pins.items():require(sha(pin['path'])==pin['sha256'],'Pinned source/input changed: '+name)
    lanes={}
    for lane,fixed_name in FIXED.items():
        anchor={name:value for name,value in point.items() if name!=fixed_name}
        identity=dict(original_source_manifest=pins['original_source_manifest']['sha256'],
            evaluator_source=pins['two_birth_driver.py']['sha256'],solver=pins['solver']['sha256'],
            overlay=pins['patch_solver.py']['sha256'],birth_rule='original_one_birth',
            selected_checkpoint_sha256=receipt['case_checkpoint_sha256'],externally_fixed_zero=fixed_name,
            target_weight_fingerprint=old['lanes']['one_birth']['target_fingerprint'])
        lanes[lane]=dict(anchor_point=anchor,initial_point=anchor,anchor_fixed_scale=point[fixed_name],
            initial_psi=receipt['normalization']['psi_child'],anchor_case=str(SELECTED),anchor_loss=receipt['loss'],
            anchor_checkpoint_sha256=receipt['case_checkpoint_sha256'],source_identity=identity,
            source_fingerprint=canon(identity),target_fingerprint=old['lanes']['one_birth']['target_fingerprint'],
            economic_changes=[f'Experimental: {fixed_name} externally fixed exactly zero, outside original [0.02, 50] bound; original one-birth opportunity and conception risk retained.',
                'Other nine original coordinates, including positive-cost initial center and freely estimated cost [0,8] in refit, retain original bounds; psi_child renormalized to completed fertility 2.1.',
                'Ten scored targets, three validation rows, renewal gate, earnings, entry distributions, timing, transfers/floors, housing and credit primitives retained.'],
            scientific_promotion=False)
    config=dict(schema='zero_fertility_taste_v1',parameters=old['parameters'],scored_moments=old['scored_moments'],
        bounds=old['bounds'],original_bounds=old['bounds'],fixed_parameters=FIXED,lanes=lanes,
        worker_command=[sys.executable,str(HERE/'worker.py'),'--config',str(HERE/'config.json')],
        hard_end_epoch=time.time()+14400,budget=dict(total_seconds=14400,final_reserve_seconds=4200,
            case_seconds=2100,maximum_stationary_solves=8,max_evaluations=16),pins=pins,
        reference_label=old['reference_label'],reference_checkpoint_sha256=old['reference_checkpoint_sha256'],
        identification='Ten scored rows for nine free coordinates per lane; psi separately normalized. Rank unverified; exact-zero choices may be nonsmooth.',
        authorization='Author-authorized preparation; launch only after lead review.',prepared_epoch=time.time(),
        scientific_promotion=False,transition_run=False)
    write(HERE/'config.json',config)
    write(HERE/'preparation_receipt.json',dict(status='prepared_not_launched',config_sha256=sha(HERE/'config.json'),
        max_evaluations_per_lane=16,max_stationary_solves_per_lane=128,anchor_checkpoint_sha256=receipt['case_checkpoint_sha256'],
        model_solves=0,scientific_promotion=False))

if __name__=='__main__':main()
