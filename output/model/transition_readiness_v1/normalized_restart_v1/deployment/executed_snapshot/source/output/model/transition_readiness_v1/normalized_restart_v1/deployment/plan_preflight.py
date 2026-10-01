import argparse,json,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--plan',required=True);p.add_argument('--output',required=True);a=p.parse_args()
root=Path(__file__).resolve().parents[5];sys.path.insert(0,str(root/'code/model/experiments/transition_readiness'))
import one_shock_floor as c
plan=json.loads(Path(a.plan).read_text());result=c.preflight(plan);c.write(Path(a.output)/'fit_preflight.json',result);print(json.dumps(result))
