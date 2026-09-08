#!/bin/bash
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=00:35:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_rebated_tax_20260908a/slurm-%j.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_rebated_tax_20260908a/mpl
cd /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a
python - <<'PY'
from pathlib import Path
import json,hashlib,csv
r=Path("/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_rebated_tax_20260908a/run");rows=[]
for name in ("tax1-equal-rebate","tax2-equal-rebate"):
 p=r/name;s=json.loads((p/"receipt.json").read_text())
 assert s["status"]=="complete" and s["contract_sha256"]=="931def1c35d3a6b4575bcb9e8cd65e345f71804d1b25d74f3671a8719b680a7d"
 for path,sha in s["artifact_sha256"].items():
  h=hashlib.sha256()
  with (p/path).open("rb") as f:
   for chunk in iter(lambda:f.read(2**20),b""):h.update(chunk)
  assert h.hexdigest()==sha,path
 rows.append(dict(policy=name,**s["endpoint"]["quantities"],**s["endpoint"]["ledger"]))
a,b=rows
result=dict(births_per_household_pct=100*(b["births_per_household"]/a["births_per_household"]-1),rooms_per_household_pct=100*(b["rooms_per_household"]/a["rooms_per_household"]-1),ownership_pp=100*(b["ownership"]-a["ownership"]),young_ownership_whole_nodes_pp=100*(b["young_ownership_whole_nodes"]-a["young_ownership_whole_nodes"]),young_mean_rooms_pct=100*(b["young_mean_rooms_whole_nodes"]/a["young_mean_rooms_whole_nodes"]-1),asset_price_pct=100*(b["asset_price"]/a["asset_price"]-1),rebate_level_change=b["equal_transfer_period_units"]-a["equal_transfer_period_units"])
with (r/"tax_endpoint_quantities.csv").open("w") as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(r/"comparison_receipt.json").write_text(json.dumps(dict(status="complete",contract_sha256="931def1c35d3a6b4575bcb9e8cd65e345f71804d1b25d74f3671a8719b680a7d",comparison="rebated tax2 versus rebated tax1",scope="2023 impact only; no causal component decomposition",effects=result),indent=2)+"\n")
print(json.dumps(result))
PY
