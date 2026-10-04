"""Derive 15 provisional recovery starts from pinned failed-parent checkpoints; zero solves."""
import hashlib,json,math,sys
from pathlib import Path
PARENT=Path('/scratch/td2248/projects/estate_birth_calibration_20261003_v3')
STAGE=Path('/scratch/td2248/projects/estate_birth_recovery_20261004_v1')
PREFIX='output/model/experiments/birth_count_choice/estate_a_recovery_20261004_v1'
PARENT_SHA='d14a39bcbb55067060ca492943f1e0067000c10733a48c3a3aff3b06e61e2afe'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def req(v,m):
    if not v:raise RuntimeError(m)
def main():
    req(sha(PARENT/'inventory.json')==PARENT_SHA,'Parent inventory drift')
    pin=read(PARENT/'inventory.json');original_path=PARENT/'source/output/model/experiments/birth_count_choice/estate_a_calibration_v1/start_plan.json'
    req(sha(original_path)==pin['start_plan_sha256'],'Parent start plan drift')
    original=read(original_path);req(original['target_fingerprint']==pin['target_fingerprint'] and original['weight_fingerprint']==pin['weight_fingerprint'],'Parent objective drift')
    rows=[]
    for task in range(10):
        arm='binary' if task<5 else 'count3';chain=task%5;folder=PARENT/'results'/f'production_{arm}_chain_{chain}'
        contract=read(folder/'run/start_contract.json');best=read(folder/'run/best_so_far.json')['best'];cases=read(folder/'run/cases.json');term=read(folder/'launcher_terminal.json');launch=read(folder/'launcher_start.json')
        req(term['exit_code']!=0 and term['slurm_job_id']==launch['slurm_job_id'] and launch['stage_inventory_sha256']==PARENT_SHA and launch['mode']=='production' and launch['arm']==arm and launch['chain']==chain,'Parent failure/source/task drift')
        req(contract['arm']==arm and contract['chain']==chain and contract['starts_file_sha256']==pin['start_plan_sha256'] and contract['target_fingerprint']==pin['target_fingerprint'] and contract['weight_fingerprint']==pin['weight_fingerprint'] and contract['bounds']==original['bounds'],'Parent contract drift')
        req(best['status']=='passed' and best in cases and best==min((c for c in cases if c['status']=='passed'),key=lambda c:c['loss']),'Parent best/cases drift')
        rows.append(dict(task=task,arm=arm,chain=chain,parameters=best['parameters'],provisional_loss=best['loss'],best_sha256=sha(folder/'run/best_so_far.json'),cases_sha256=sha(folder/'run/cases.json'),contract_sha256=sha(folder/'run/start_contract.json'),terminal_sha256=sha(folder/'launcher_terminal.json'),launcher_start_sha256=sha(folder/'launcher_start.json')))
    bounds=original['bounds'];b=min(rows[:5],key=lambda r:(r['provisional_loss'],r['task']))['parameters'];c=min(rows[5:],key=lambda r:(r['provisional_loss'],r['task']))['parameters']
    def clip(row):return {k:min(bounds[k][1],max(bounds[k][0],float(v))) for k,v in row.items()}
    one=dict(b);two={k:(b[k]+c[k])/2 for k in c}
    three=dict(c);three['beta_annual']=.94;three['first_birth_fixed_cost']=c['first_birth_fixed_cost']*.5
    four=dict(c);four['kappa_fert']=c['kappa_fert']*3;four['kappa_fert_continuation']=c['kappa_fert_continuation']*3
    five=dict(c);five['theta0']=c['theta0']*2;five['psi_child']=c['psi_child']*1.5;five['h_P']=c['h_P']-.25
    for arm,starts in [('binary',[r['parameters'] for r in rows[:5]]),('count3',[r['parameters'] for r in rows[5:]]+[clip(s) for s in (one,two,three,four,five)])]:
        req(len({json.dumps(s,sort_keys=True) for s in starts})==len(starts),'Duplicate starts')
        req(all(set(s)==set(bounds) and all(math.isfinite(v) and bounds[k][0]<=v<=bounds[k][1] for k,v in s.items()) for s in starts),'Out of bounds')
        plan=dict(original);plan['per_chain']=dict(original['per_chain'],wall_seconds=43200,absolute_stop_epoch=1791122400);plan['solve_budget']=dict(original['solve_budget'],actual_chain_wall_cap_hours=12,total_chains=15,search_max_GE_calls=7500,final_fresh_GE_calls=15)
        plan.update(starts=starts,recovery_arm=arm,parent_job_id='19127370',parent_inventory_sha256=PARENT_SHA,parent_receipt_root='/work/parent',parent_rows=rows,original_start_plan_sha256=sha(original_path),provisional_seed_source='/work/parent/inventory.json',provisional_seed_source_sha256=PARENT_SHA,provisional_seed=dict(source='failed_parent_saved_best',verified=False),deterministic_seed=20261004,no_auto_retry=True,no_scientific_adoption=True)
        path=STAGE/'source'/PREFIX/f'plan_{arm}.json';path.parent.mkdir(parents=True,exist_ok=True);req(not path.exists(),'Refusing existing recovery plan');path.write_text(json.dumps(plan,sort_keys=True,indent=2)+'\n')
    receipt=dict(status='provisional_parent_saved_best_only',parent_inventory_sha256=PARENT_SHA,plan_binary_sha256=sha(STAGE/'source'/PREFIX/'plan_binary.json'),plan_count3_sha256=sha(STAGE/'source'/PREFIX/'plan_count3.json'),parent_rows=rows,no_native_postcheck_yet=True)
    (STAGE/'control').mkdir(exist_ok=True);(STAGE/'control/starts_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='parent_rows'}))
if __name__=='__main__':main()
