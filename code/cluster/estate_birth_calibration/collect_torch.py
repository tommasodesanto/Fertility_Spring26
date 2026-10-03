"""Read-only matched collector: reject mixed source/start/target fingerprints."""
import argparse,hashlib,json
from pathlib import Path

def read(path):return json.loads(path.read_text())
def collect(stage,mode):
    inv=read(stage/'inventory.json');inventory_sha=hashlib.sha256((stage/'inventory.json').read_bytes()).hexdigest();rows=[]
    for arm,cap in [('binary',1),('count3',3)]:
        for chain in ([0] if mode=='smoke' else range(5)):
            run=stage/'results'/f'{mode}_{arm}_chain_{chain}'/'run'
            if not (run/'completed.json').exists():
                rows.append(dict(arm=arm,chain=chain,status='pending_or_incomplete',heartbeat=read(run/'heartbeat.json') if (run/'heartbeat.json').exists() else None));continue
            contract=read(run/'start_contract.json');done=read(run/'completed.json');launcher=read(run.parent/'launcher_start.json')
            assert launcher['stage_inventory_sha256']==inventory_sha
            assert launcher['mode']==mode and launcher['arm']==arm and launcher['chain']==chain
            flags=dict(birth_count_choice_enabled=True,birth_count_choice_cap=cap,bequest_net_of_selling_cost=True,estate_flow_net_of_selling_cost=True)
            assert contract['experiment_flags']==flags
            assert contract['birth_cap']==cap and contract['selected_source_sha256']==inv['selected_source_sha256']
            assert contract['starts_file_sha256']==done['starts_file_sha256']==inv['start_plan_sha256']
            assert all(contract[k]==done[k]==inv[k] for k in ('target_fingerprint','weight_fingerprint'))
            if done['status']=='selected_numerically_verified':
                assert done['selected_postcheck']['experiment_flags']==flags
                assert done['selected_postcheck']['observer_contract']=='estate_a_postsaving_net_selling_cost_v1'
                assert len(done['target_fit'])==14 and len(done['parameters'])==31 and len(done['repeat']['standard_plot_hashes'])==17
            rows.append(dict(arm=arm,chain=chain,**done))
    best={arm:min([r for r in rows if r['arm']==arm and r['status']=='selected_numerically_verified'],key=lambda r:r['native_loss'],default=None) for arm in ('binary','count3')}
    return dict(mode=mode,target_fingerprint=inv['target_fingerprint'],weight_fingerprint=inv['weight_fingerprint'],chains=rows,best_by_arm=best,no_adoption=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',type=Path,required=True);p.add_argument('--mode',choices=('smoke','production'),required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not a.out.exists(),'Refusing existing collection';data=collect(a.stage,a.mode);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(chains=len(data['chains']),verified=sum(r['status']=='selected_numerically_verified' for r in data['chains']))))
