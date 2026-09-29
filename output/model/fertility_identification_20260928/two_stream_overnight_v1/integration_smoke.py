"""Two real evaluations per stream through the exact owned dispatcher, on Torch.

The synthetic tests cover the full search loop. This verifies model wiring,
normalization, actual exported artifacts and repeated dispatch for both models.
"""
from __future__ import annotations
import argparse,csv,math,os,time
from pathlib import Path
import search
from worker import HERE,read,write,sha,require,verify_config


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--lane',choices=('one_birth','two_birth'),required=True)
    args=parser.parse_args();config_path=HERE/'config.json';c=verify_config(config_path)
    require(sha(HERE/'search.py')==c['pins']['search']['sha256'],'Controller source changed')
    out=HERE/'integration_smoke_v1'/args.lane
    stream=search.Stream(c,args.lane,out,sha(config_path))
    # Additional hard cap on this two-call smoke; the production budget is not modified.
    stream.end=min(stream.end,time.time()+4200)
    stream.cutoff=stream.end
    lane=c['lanes'][args.lane]
    with (Path(lane['anchor_case'])/'target_fit.csv').open(newline='') as f:
        fit={r['moment']:r for r in csv.DictReader(f)}
    expected=[math.sqrt(float(fit[k]['weight']))*float(fit[k]['gap']) for k in c['scored_moments']]
    results=[]
    for index in range(2):
        result=stream.evaluate(lane['initial_point'],lane['initial_psi'],f'integration_repeat{index}',final=index>0)
        require(result is not None,'Real integration evaluation failed')
        residual=max(abs(a-b) for a,b in zip(expected,result['residuals']))
        loss=abs(result['loss']-lane['anchor_loss'])
        require(residual<=.01 and loss<=.05,'Anchor numerical replay differs beyond retained screens')
        results.append(dict(success_path=result['success_path'],success_sha256=result['success_sha256'],
            checkpoint_sha256=result['checkpoint_sha256'],maximum_weighted_residual_difference=residual,
            loss_difference=loss,model_evaluations=result['model_evaluations'],elapsed_seconds=result['elapsed_seconds']))
    write(out/'SMOKE.json',dict(status='passed',lane=args.lane,config_sha256=sha(config_path),
        evaluations=results,real_model_evaluations=2,source_fingerprint=lane['source_fingerprint'],
        target_fingerprint=lane['target_fingerprint'],full_tables_and_17_plots_per_evaluation=True,
        elapsed_seconds=time.time()-stream.start,scientific_promotion=False))


if __name__=='__main__':main()
