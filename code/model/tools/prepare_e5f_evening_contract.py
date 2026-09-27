#!/usr/bin/env python3
"""Torch-only, standard-library contract builder; never imports/runs the model.

Run inside the staged project's exact-path container after source review. All
fourteen targets and nine inherited bounds come from the byte-pinned ancestry.
This builder prepares an immutable contract; it cannot authorize search.
"""
from __future__ import annotations
import argparse,copy,hashlib,json,math,os,sys
from pathlib import Path
START=1790539920
END=1790561520
FREE=('H0','beta_annual','chi','first_birth_fixed_cost','kappa_fert','kappa_fert_continuation','theta0','delta_alpha_jump','child_benefit_curvature','tenure_choice_kappa')
VALIDATION=('nchs_share30','family_rooms','old_dispersion')
BLOCKS={'fertility':('cps_childlessness','cps_exactly_one','nchs_mean_age','early_fertility'),'housing':('mean_rooms','ownership_30_55','first_birth_rooms','recent_parent_ownership'),'wealth':('wealth_earnings','bequest_wealth')}
SCALE_FLOORS={'cps_childlessness':.1,'cps_exactly_one':.1,'nchs_mean_age':1.,'early_fertility':.1,'mean_rooms':1.,'ownership_30_55':.1,'first_birth_rooms':1.,'recent_parent_ownership':.1,'wealth_earnings':1.,'bequest_wealth':.01}

def read(p):return json.loads(Path(p).read_text())
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for part in iter(lambda:f.read(1<<20),b''):h.update(part)
    return h.hexdigest()
def canon(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def pin(p):p=Path(p).resolve(strict=True);return dict(path=str(p),sha256=sha(p))
def write(p,x):
    with Path(p).open('x') as f:f.write(json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+'\n')

def lane_objectives(original,tenure_lower=.001,tenure_upper=.1):
    rows={r['restriction_id']:r for r in original['target_rows']}
    scored=set().union(*map(set,BLOCKS.values()))
    assert set(rows)==scored|set(VALIDATION)|{'initial_normalization'} and len(rows)==14
    assert rows['initial_normalization']['actual_weight'] is None
    assert {r['parameter'] for r in original['parameter_restrictions']}==set(FREE[:-1])
    results={}
    for lane in ('primary','identity','block'):
        obj=copy.deepcopy(original)
        obj['parameter_restrictions'].append(dict(parameter='tenure_choice_kappa',lower=tenure_lower,upper=tenure_upper,status='author-authorized exploratory positive log interval; not externally estimated'))
        obj['evening_weight_contract']=dict(lane=lane,definition='inherited active weights' if lane=='primary' else 'identity after fixed target-relative row scaling' if lane=='identity' else 'one third times each block mean of inherited standardized errors',blocks=BLOCKS,scale_floors=SCALE_FLOORS)
        for r in obj['target_rows']:
            name=r['restriction_id'];inherited=r['actual_weight'];r['inherited_weight']=inherited
            if name=='initial_normalization':r['role']='normalization';continue
            if name in VALIDATION:r['actual_weight']=0.;r['role']='validation';continue
            assert math.isfinite(inherited) and inherited>0
            block=next(k for k,values in BLOCKS.items() if name in values);scale=max(abs(float(r['target'])),SCALE_FLOORS[name])
            r.update(role='scored',weight_block=block,relative_scale=scale)
            if lane=='identity':r['actual_weight']=1/(scale*scale)
            elif lane=='block':r['actual_weight']=inherited/(3*len(BLOCKS[block]))
        results[lane]=obj
    return results

def initial_design(point,bounds,widths):
    rows=[];omitted=[]
    for value,label in ((.05,'tenure_anchor_0.05'),(.004,'tenure_profile_0.004'),(.006,'tenure_profile_0.006')):
        p=dict(point,tenure_choice_kappa=value);assert bounds['tenure_choice_kappa'][0]<=value<=bounds['tenure_choice_kappa'][1]
        rows.append(dict(label=label,point=p,scope='fixed seed, benefit renormalized'))
    for k in ('chi','delta_alpha_jump','child_benefit_curvature'):
        lo,hi=bounds[k];step=min(.2*widths[k],.01*(hi-lo))
        for sign,label in ((-1,'minus'),(1,'plus')):
            value=point[k]+sign*step
            if value<lo or value>hi:omitted.append(dict(parameter=k,direction=label,reason='outside inherited bounds'));continue
            p=dict(point);p[k]=value;rows.append(dict(label=f'{k}_seed_profile_{label}',point=p,scope='fixed seed, benefit renormalized'))
    assert len({canon(r['point']) for r in rows})==len(rows) and all(r['point']!=point for r in rows)
    return rows,omitted

def build(a):
    if sys.platform!='linux' or not os.environ.get('SLURM_JOB_ID','').isdigit():raise RuntimeError('Torch Slurm preparation only')
    root=a.source_root.resolve(strict=True);out=a.output.resolve()
    if out.exists():raise FileExistsError('Use a fresh contract directory')
    ancestry=read(a.ancestry_contract);objective_pin=pin(a.inherited_objective)
    assert objective_pin==ancestry['objective'],'Inherited objective must be exact ancestry file'
    original=read(a.inherited_objective);seed=read(a.seed_receipt)
    assert seed['target_system_sha256']==objective_pin['sha256'],'Seed uses another target system'
    point={k:float(seed['point'][k]) for k in FREE[:-1]};point['tenure_choice_kappa']=.005
    objectives=lane_objectives(original,a.tenure_lower,a.tenure_upper)
    bounds={r['parameter']:(r['lower'],r['upper']) for r in objectives['primary']['parameter_restrictions']}
    assert all(math.isfinite(point[k]) and lo<=point[k]<=hi for k,(lo,hi) in bounds.items())
    widths={k:float(ancestry['proposal_widths'][k]) for k in FREE[:-1]};widths['tenure_choice_kappa']=math.log(2)
    design,omitted=initial_design(point,bounds,widths)
    names=sorted(p.name for p in a.diagnostics_directory.glob('*.png'));assert len(names)==len(set(names))==17
    tools=root/'code/model/tools';files={k:pin(tools/name) for k,name in {'driver':'run_e5f_evening_calibration.py','runtime':'e5f_evening_calibration_runtime.py','builder':'prepare_e5f_evening_contract.py','controller_tests':'test_e5f_evening_calibration.py','runtime_tests':'test_e5f_evening_calibration_runtime.py'}.items()}
    files['recovery_search']=pin(a.recovery_search)
    files['weight_review']=pin(root/'docs/model/e5f_evening_weight_review_20260927.md')
    files['measurement_review']=pin(root/'docs/model/e5f_target_measurement_review_20260927.md')
    files['first_birth_builder']=pin(root/'code/data/psid_followup_mar2026/sa_rooms_first_birth_v2.do')
    files['bequest_builder']=pin(root/'code/data/scf/build_bequest_flow_2007.py')
    inventory={str(p.relative_to(root)):sha(p) for p in sorted((root/'code/model').rglob('*.py')) if '__pycache__' not in p.parts}
    assert inventory and all((root/path).is_file() for path in inventory)
    out.mkdir(parents=True)
    write(out/'source_manifest.json',dict(files=inventory,scope='All active code/model Python sources; inherited tool/data inventories separately authenticated by native ancestry runtime',file_count=len(inventory)))
    lanes={}
    for lane,obj in objectives.items():
        path=out/f'objective_{lane}.json';write(path,obj);lanes[lane]=dict(objective=pin(path),canonical_sha256=canon(obj),target_weight_fingerprint=canon(obj['target_rows']))
    fixed=copy.deepcopy(ancestry['fixed']);fixed.pop('tenure_choice_kappa',None);fixed.update(due=True,delta_alpha=0.,sigma=2.)
    contract=dict(schema='e5f_evening_v1',source_root=str(root),source_manifest=pin(out/'source_manifest.json'),files=files,lanes=lanes,native_ancestry_contract=pin(a.ancestry_contract),reference_case=str(a.reference_case.resolve(strict=True)),seed_receipt=pin(a.seed_receipt),inherited_objective=objective_pin,initial_point=point,initial_design=design,omitted_seed_profiles=omitted,normalization=copy.deepcopy(ancestry['normalization']),fixed=fixed,validation_rows=list(VALIDATION),standard_diagnostic_names=names,log_parameters=['tenure_choice_kappa'],proposal_widths=widths,seed=20260927,budget=dict(absolute_start_epoch=START,absolute_end_epoch=END,total_seconds=21600,workers=24,max_objective_cases=384,max_search_cases=366,maximum_native_prerequisite_cases=2,objective_cap_seconds=900,search_reserve_seconds=2700,export_reserve_seconds=300),economic_changes=copy.deepcopy(ancestry['economic_changes'])+[dict(status='author-adopted tonight',change='Native DUE stayer credit; tenure choice scale searched; three specified validation moments receive zero objective weight; three separate weighting lanes')],pending_observer_mismatches=copy.deepcopy(ancestry['pending_observer_mismatches']),search_design='Six lane smokes; fixed seed tenure anchors and small coordinate profiles consume search slots; then deterministic broad/joint/coordinate proposals; each lane selected separately; two fresh repeats per selected lane; no optimizer convergence claim',authorization='Preparation only. Search requires separately pinned lead approval of the exact six-smoke receipt.',budget_accounting=dict(controller_smokes=6,native_prerequisites_at_most=2,search_at_most=366,selected_repeats=6,total_at_most=380))
    assert contract['normalization']['maximum_stationary_solves']==23
    contract['tenure_seed_status']='0.005 is the inherited fixed tenure scale at the selected nine-coordinate seed; the evening contract searches this coordinate in its declared positive interval.'
    contract['provenance_corrections']={
        'scope':'New-contract metadata correction only; targets, samples, observers and immutable ancestry are unchanged.',
        'source':files['measurement_review'],
        'first_birth_rooms':dict(builder=files['first_birth_builder'],arm='A2h',receipt_estimate=1.465293280235685,adopted_target=1.465,
            event_clock='K = year - f_c_y',comparison='calendar years +3/+4 since first birth relative to calendar years -3/-2, not interview numbers',
            model_limit='Matched-state destination housing contrast after one four-year model period remains a proxy; it does not reproduce empirical prepath, cohort/event weighting or baseline head/spouse selection.'),
        'bequest_wealth':dict(builder=files['bequest_builder'],authority='Adopted 2007 SCF receipt remains authority for selected row; builder existence corrects frozen metadata, and its historical source fingerprint was not established by the measurement review.',
            model_limit='Empirical annual child-directed mortality-weighted flow differs from model positive post-saving estates over all child-history states; recipient allocation remains an acknowledged approximation.')}
    write(out/'contract.json',contract)
    write(out/'preparation_receipt.json',dict(status='prepared_no_model_import_or_solve',contract=pin(out/'contract.json'),source_manifest=contract['source_manifest'],lanes=lanes,weighted_rows=10,displayed_rows=14,searched_coordinates=10,fitted_quantities_including_normalization=11,absolute_end_epoch=END))
    print(json.dumps(dict(contract=str(out/'contract.json'),sha256=sha(out/'contract.json'),total_cases_at_most=380,absolute_end_epoch=END)))

def main():
    p=argparse.ArgumentParser()
    for name in ('source-root','ancestry-contract','reference-case','seed-receipt','inherited-objective','recovery-search','diagnostics-directory','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--tenure-lower',type=float,default=.001);p.add_argument('--tenure-upper',type=float,default=.1);build(p.parse_args())
if __name__=='__main__':main()
