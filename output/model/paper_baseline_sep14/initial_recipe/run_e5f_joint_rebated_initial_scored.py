"""Score the joint initial closure through the unchanged full target observer."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import sys


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--helper',type=Path,required=True)
    parser.add_argument('--helper-sha256',required=True)
    parser.add_argument('--joint',type=Path,required=True)
    parser.add_argument('--joint-sha256',required=True)
    parser.add_argument('--template',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--proposal',type=Path)
    args=parser.parse_args()
    spec=importlib.util.spec_from_file_location('frozen_rebate_helper',args.helper)
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    if helper.sha(args.helper)!=args.helper_sha256 or helper.sha(args.joint)!=args.joint_sha256:
        raise ValueError('Initial numerical adapter source pin changed')
    packet=helper.saved_packet(args.template)
    restrictions=helper.validate_scientific_contract(packet)
    original_load=helper.load_module
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    def load(name,path):
        module=original_load(name,path)
        if name=='saved_rebated_scored_wrapper':
            original_child=module.run_child
            def child(command,**kwargs):
                if len(command)>2 and command[2]=='raw' and Path(command[1]).resolve()==args.helper.resolve():
                    root=command[command.index('--source-root')+1]
                    driver=command[command.index('--driver')+1]
                    tail=command[command.index('--')+1:]
                    command=[sys.executable,'-B',str(args.joint),'--source-root',root,
                        '--driver',driver,'--progress',str(out/'joint_progress.json'),'--',*tail]
                return original_child(command,**kwargs)
            module.run_child=child
        return module
    helper.load_module=load
    packet['initial']=copy.deepcopy(packet['initial'])
    packet['initial'].update(initial_root_method='joint_price_psi_equal_rebate',
        joint_adapter_sha256=args.joint_sha256,maximum_fixed_price_evaluations_per_root=21)
    proposal=(json.loads(args.proposal.read_text()) if args.proposal else copy.deepcopy(packet['plan']['resume_proposal']))
    proposal.setdefault('repetitions',1)
    if not args.proposal:proposal['case_id']='joint_rebated_initial_smoke'
    helper.write(out/'launch_contract.json',dict(helper_sha256=args.helper_sha256,
        joint_sha256=args.joint_sha256,proposal=proposal,objective_sha256=packet['run']['working_objective']['canonical_sha256']))
    result=helper.candidate(packet,restrictions,proposal,out/'case')
    helper.write(out/'candidate_result.json',result)
    if result['status']!='verified':
        helper.write(out/'summary.json',dict(status='failed_joint_rebated_initial',result=result));raise SystemExit(2)
    helper.write_selected_tables(out,result)
    summary=dict(status='verified_rebated_initial_smoke',method='joint_price_fertility_rebate',
        source_root=str(packet['source_root']),loss=result['loss'],checkpoint=result['accounting'][-1])
    for name in ('latest_completed.json','best_so_far.json','summary.json'):helper.write(out/name,summary)
    print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
