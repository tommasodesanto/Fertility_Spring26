"""Build immutable overnight contracts on Torch from the verified capped seed."""
from pathlib import Path
import hashlib,json,sys

def digest(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(Path(p).read_text())
def main():
    root=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a')
    batch=root/'batches/night_surprises_20260912'
    seed=Path('/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/capped_beta_099_20260911/results/cases/selected_exact_repetitions')
    c=read(batch/'terminal_template.json');c['schema']='e5f_candidate_terminal_v1'
    for name,path in {'initial_contract':seed/'initial_contract.json','initial_checkpoint':seed/'evaluation/raw/repetition_02/initial_state.pkl.gz','initial_summary':seed/'evaluation/raw/repetition_02/summary.json'}.items():
        c[name]=dict(path=str(path),sha256=digest(path))
    c['originating_source_commit']=read(seed/'initial_contract.json')['source_commit']
    c['source_sha256']=dict(read(seed/'initial_contract.json')['source_sha256'])
    c['source_sha256'].update({str(p.relative_to(root)):digest(p) for p in (root/'code/model').rglob('*.py')})
    c['psi_change_from_initial']=0.
    sys.path[:0]=[str(root/'code/model/tools'),str(root/'code/model')]
    import run_e5f_candidate_terminal as terminal
    terminal.validate_contract(c);terminal.verify_sources(c);terminal.validate_initial_sources(c,read(seed/'initial_contract.json'))
    score=seed/'evaluation/scored_repetition_01/score.json'
    rc=read(batch/'history_template.json')['root_controls'];rc['max_evaluations']=8
    pins=[*batch.glob('*.py'),batch/'empirical_blocks.csv',score]
    p=dict(source_root=str(root),output_root=str(batch/'results'),total_seconds=35400,policy_reserve_seconds=7200,
        terminal_template=c,initial_score_path=str(score),target_fingerprint='c0e266d3a0d430343c469d780d1aedb45fa87f8763c9c938889e0c37daa31de2',
        empirical_blocks=str(batch/'empirical_blocks.csv'),maximum_trials_per_window=4,seed_steps=[-.005,-.015,-.03],
        warm_price_2007=.7071194654472099,warm_pension_2007=3.512600520619363,history_root_controls=rc,
        fertility_fit_tolerance=.005,file_sha256={str(f):digest(f) for f in pins},
        outside_entry_status='diagnostic_outstanding_not_estimated',outside_origin_entry_share=.169,
        author_authorization='Overnight fitting with current rooms target, beta estimated capped at .99, and conditional policy experiments; September 12 conversation',
        maximum_terminal_solves_per_arm=16,maximum_forecast_root_mappings_per_trial=72,
        production_eligible=False)
    assert read(score)['contract_sha256']==p['target_fingerprint']
    (batch/'plan.json').write_text(json.dumps(p,indent=2)+'\n')
    print(json.dumps(dict(plan=str(batch/'plan.json'),sha256=digest(batch/'plan.json'),source_files=len(c['source_sha256']))))
if __name__=='__main__':main()
