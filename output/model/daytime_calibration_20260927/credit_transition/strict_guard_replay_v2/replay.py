"""Replay the strict gate on authenticated saved states; zero model solves."""
import copy,gzip,hashlib,json,pickle,signal,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[5]
sys.path[:0]=[str(ROOT/'code/model/tools'),str(ROOT/'code/model')]
import numpy as np
import run_dynamic_population_transition as calendar
OUT=Path(__file__).resolve().parent
BASE=ROOT/'tmp/e5f_overnight_local_20260927/portable/night_launch_v4/primary_continuation/search/de_0093/case/initial_state.pkl.gz'
DATED=OUT.parent/'root_smoke_v1/run/candidate_packets/evaluation_0008/date_0000.pkl.gz'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def load(p,h):
 assert sha(p)==h
 with gzip.open(p,'rb') as f:return pickle.load(f)
signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('90-second saved-state inspection cap')))
signal.alarm(90)
a=load(BASE,'090c9ebda662bf7837c4f4cf1d816159bc9d203a9babe7c70d00d0c4be1e575e')
b=load(DATED,'3210426153bc8a3241b2c74c355459c7b5a997bee4403ed5f924e7ca670d9306')
g=np.asarray(a['stationary_g_pre']);before=g.tobytes();result={}
for name,packet in [('baseline',a),('dated',b)]:
 P=copy.deepcopy(packet['parameters']);P.native_exact_inherited_distribution=True
 P.native_inherited_distribution_evidence_dir=str(OUT/name)
 policy=packet['evaluation'].policy
 try:
  out,moved=calendar.gate_pre_fertility_distribution(g,policy,P,a['b_grid'],packet['shared'])
  assert out.tobytes()==before and moved==0
  result[name]={'status':'preserved','projection':moved}
 except calendar.InheritedDistributionInfeasible as e:
  result[name]={'status':'rejected_without_modification','dead_mass':e.dead_mass,'evidence_path':e.evidence_path}
 assert g.tobytes()==before
result.update(no_model_solves=True,initial_array_unchanged=True,calendar_source_sha256=sha(Path(calendar.__file__)),limitation='Guard replay only; no changed-price Bellman solve or new equilibrium. Baseline positive infeasible tail is not silently dropped.')
(OUT/'complete.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
