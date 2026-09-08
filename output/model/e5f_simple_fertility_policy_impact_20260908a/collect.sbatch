#!/bin/bash
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:35:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_policy_impact_20260908a/slurm-%A-%a.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_policy_impact_20260908a/mpl
cd /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a
python - <<'PY'
from pathlib import Path
import csv,json,hashlib
r=Path("/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_policy_impact_20260908a/run")
rows=[]
for name in ("baseline","dependent-child-ltv95","supply-plus-20"):
 p=r/name;receipt=json.loads((p/"receipt.json").read_text())
 assert receipt["status"]=="complete" and receipt["contract_sha256"]=="1f7beed6c5ededada14cce5b7b480aee7fcc4619a3963fb09b0b6b7620ccd3ca"
 for path,sha in receipt["artifact_sha256"].items():
  h=hashlib.sha256()
  with (p/path).open("rb") as f:
   for block in iter(lambda:f.read(2**20),b""):h.update(block)
  assert h.hexdigest()==sha,path
 rows.append(dict(policy=name,**receipt["evaluations"][0]["quantities"]))
base=rows[0];comp=[]
for row in rows[1:]:
 comp.append(dict(policy=row["policy"],births_pct=100*(row["births_per_household"]/base["births_per_household"]-1),rooms_pct=100*(row["rooms_per_household"]/base["rooms_per_household"]-1),ownership_pp=100*(row["ownership"]-base["ownership"]),asset_price_pct=100*(row["asset_price"]/base["asset_price"]-1)))
for name,data in (("case_quantities.csv",rows),("impact_comparison.csv",comp)):
 with (r/name).open("w") as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
(r/"comparison_receipt.json").write_text(json.dumps(dict(status="complete",contract_sha256="1f7beed6c5ededada14cce5b7b480aee7fcc4619a3963fb09b0b6b7620ccd3ca",scope="2023 impacts; no future transition",comparison=comp),indent=2)+"\n")
print(json.dumps(comp))
PY
