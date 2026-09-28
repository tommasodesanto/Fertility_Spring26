#!/usr/bin/env python3
"""Prepare a separate overnight contract on Torch. No model imports or solves."""
import argparse,copy,hashlib,importlib.util,json,os,sys
from pathlib import Path
END=1790596800

def read(p):return json.loads(Path(p).read_text())
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for x in iter(lambda:f.read(1<<20),b''):h.update(x)
    return h.hexdigest()
def pin(p):p=Path(p).resolve(strict=True);return dict(path=str(p),sha256=sha(p))
def write(p,x):
    with Path(p).open('x') as f:f.write(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def load(path):
    s=importlib.util.spec_from_file_location('night_contract_helpers',path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def build(a):
    assert sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm preparation only'
    root=a.source_root.resolve(strict=True);out=a.output.resolve();assert not out.exists()
    old=read(a.evening_contract);assert old['schema']=='e5f_evening_v1'
    for record in list(old['files'].values())+[old['source_manifest']]:assert sha(record['path'])==record['sha256'],record['path']
    assert 0<END-a.start_epoch<=86400
    # New source tree may contain new controller files; every retained scientific
    # and reporting source must remain byte-identical to the evening inventory.
    inventory=read(old['source_manifest']['path'])['files']
    for path,digest in inventory.items():assert sha(root/path)==digest,'Evening source changed: '+path
    helper=load(root/'code/model/tools/run_e5f_night_calibration.py')
    anchor=a.anchor_case.resolve(strict=True);repeats=[a.anchor_repeat_one.resolve(strict=True),a.anchor_repeat_two.resolve(strict=True)]
    for repeated in repeats:helper.compare_tables(anchor,repeated)
    receipt=read(anchor/'receipt.json');assert receipt['point']==read(repeats[0]/'receipt.json')['point']==read(repeats[1]/'receipt.json')['point']
    history_path=a.timeout_history or root/'output/model/evening_calibration_20260927/cluster/final_review/checkpoint.json'
    history=read(history_path);design=[]
    for name in ('initial_0044_block','initial_0186_primary'):
        matches=[row for row in history['records'] if row['case']==name]
        assert len(matches)==1 and matches[0]['status']=='censored_timeout'
        design.append(dict(label='timeout_replay_'+name,point=matches[0]['point'],scope='Exact retained point in each lane under1800-second cap; six search slots total'))
    assert receipt['normalization_inputs']==old['normalization'],'Preserve exact canonical normalization start/step'
    assert sha(anchor/'initial_state.pkl.gz')==receipt['case_checkpoint_sha256']
    names=sorted(p.name for p in (repeats[0]/'standard_diagnostics').glob('*.png'));assert names==sorted(old['standard_diagnostic_names'])
    c=copy.deepcopy(old);c.update(schema='e5f_night_v1',evening_contract=pin(a.evening_contract),source_root=str(root),initial_point=receipt['point'],initial_design=design,timeout_replay_history=pin(history_path),omitted_seed_profiles=[],seed=20260928,
        anchor=dict(case_path=str(anchor),receipt=pin(anchor/'receipt.json'),target_fit=pin(anchor/'target_fit.csv'),parameters=pin(anchor/'parameters.csv'),repeats=[dict(case_path=str(p),receipt=pin(p/'receipt.json'),target_fit=pin(p/'target_fit.csv'),parameters=pin(p/'parameters.csv')) for p in repeats]),
        proposal_widths={key:value*.5 for key,value in old['proposal_widths'].items()},
        budget=dict(absolute_start_epoch=a.start_epoch,absolute_end_epoch=END,total_seconds=END-a.start_epoch,workers=24,max_objective_cases=746,max_search_cases=720,maximum_diagnostic_objectives=12,objective_cap_seconds=1800,search_reserve_seconds=3600,export_reserve_seconds=600),
        search_design='Local/moderate multivariate proposals around separate lane incumbents starting from the repeated evening block0347 point. Eight-mode cycle: one moderate joint, three local joint, housing subspace, fertility subspace, two coordinate. No full-box proposals tonight. Full scientific bounds remain unchanged; this computational restriction is not a reachability claim.',
        search_design_evidence='Evening broad coverage:77 of84 rejected by unchanged housing gates,seven timed out,zero successes. All61 timeout ledgers contained4-8 certified GEs,none a first-GE timeout. Night doubles only objective wall-time to1800seconds;23-solve cap and scientific gates unchanged. Two exact timeout points replay in allthree lanes within720 search slots.',
        budget_accounting=dict(search_at_most=720,controller_smokes=6,final_repeats_at_most=8,other_diagnostic_objectives_at_most=12,total_at_most=746,hourly_plot_exports='zero-solve saved-checkpoint renderer; consumes worker/time capacity, not objective slots'),
        authorization='Prepared only. Six exact-loop lane smokes must match all14 anchor physical rows and all31 parameters; search additionally requires a separately pinned approval receipt.',
        diagnostic_status='No additional diagnostic objectives dispatched by this version; twelve slots remain reserved, not automatically used.')
    c['files'].update(driver=pin(root/'code/model/tools/run_e5f_night_calibration.py'),builder=pin(root/'code/model/tools/prepare_e5f_night_contract.py'),controller_tests=pin(root/'code/model/tools/test_e5f_night_calibration.py'),evening_driver=pin(root/'code/model/tools/run_e5f_evening_calibration.py'))
    assert c['files']['runtime']==old['files']['runtime'] and c['normalization']==old['normalization'] and c['fixed']==old['fixed']
    out.mkdir(parents=True)
    new_inventory={str(p.relative_to(root)):sha(p) for p in sorted((root/'code/model').rglob('*.py')) if '__pycache__' not in p.parts}
    write(out/'source_manifest.json',dict(files=new_inventory,file_count=len(new_inventory),retained_evening_sources_verified=True));c['source_manifest']=pin(out/'source_manifest.json')
    write(out/'contract.json',c)
    write(out/'preparation_receipt.json',dict(status='prepared_unverified_numerically',contract=pin(out/'contract.json'),anchor=c['anchor'],global_end_epoch=END,search_cutoff_epoch=END-3600,repeat_cutoff_epoch=END-600,original_normalization_preserved=True))
    print(json.dumps(dict(contract=pin(out/'contract.json'),status='prepared_no_model_solve',search_cutoff_epoch=END-3600,global_end_epoch=END)))

def main():
    p=argparse.ArgumentParser()
    for name in ('source-root','evening-contract','anchor-case','anchor-repeat-one','anchor-repeat-two','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--timeout-history',type=Path)
    p.add_argument('--start-epoch',type=int,required=True);build(p.parse_args())
if __name__=='__main__':main()
