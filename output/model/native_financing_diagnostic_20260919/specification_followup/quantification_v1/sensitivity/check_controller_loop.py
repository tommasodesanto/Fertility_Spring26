from pathlib import Path
import sys,tempfile,json,types,threading
sys.path.insert(0,str(Path(__file__).resolve().parents[6]/'code/model/tools'))
import run_e5f_income_fit_sensitivity as m
root=Path.cwd(); prior=json.loads((root/'output/model/native_financing_diagnostic_20260919/overnight/final_search/cases.json').read_text()); source=next(x for x in prior if x['case']==60)['receipt'];seed=m.structural_parameters(source);calls=[]
def evaluator(plan_path,out,parameters,**kw):
 calls.append(str(out));r=json.loads(json.dumps(source));r['loss']=m.ANCHOR_LOSS+sum((parameters[k]-seed[k])**2 for k in seed);r['_evaluation']=str(out/'evaluation')
 for p in r['parameters']:
  if p['parameter'] in parameters:p['estimate']=parameters[p['parameter']]
 return r
old=types.SimpleNamespace(adapter_evaluator=evaluator,validate_score_contract=lambda r,p:None,STOP=threading.Event())
m.load_old_controller=lambda p:old
m.load_old_adapter=lambda p:types.SimpleNamespace(validate_plan=lambda p,**kw:None)
m.artifact_signature=lambda *a,**kw:{'fixture':'same'}
anchor={'score':source,'parameters':seed,'psi':m.ANCHOR_PSI,'evaluation':'fixture','score_path':'fixture/score.json','numeric_signature':m.numeric_signature(source),'array_signature':{'fixture':'same'}}
m.resolve_anchor=lambda *a,**kw:anchor
plan=json.loads((root/'output/model/native_financing_diagnostic_20260919/overnight/plan.remote.json').read_text())
with tempfile.TemporaryDirectory() as t:
 p=Path(t);(p/'plan.json').write_text(json.dumps(plan))
 args=types.SimpleNamespace(mode='smoke',plan=p/'plan.json',selected_summary=p/'selected.json',output=p/'smoke',case_timeout=900,workers=8,global_timeout=300,finish_utc='2026-09-21T11:30:00Z',smoke_summary=None,expected_plan_sha256=None)
 a=m.run_panel(args);assert a['evaluations']==3 and len(calls)==3
 args.mode='production';args.output=p/'production';args.smoke_summary=p/'smoke/smoke_summary.json'
 b=m.run_panel(args);assert b['objective_evaluations']==28 and len(calls)==28,(b.get('objective_evaluations'),len(calls))
 assert b['selected']['objective']==m.ANCHOR_LOSS
 assert len(json.loads((p/'production/finite_difference_jacobian.json').read_text())['columns'])==24
 print(json.dumps({'status':'passed','mock_objective_calls':len(calls),'smoke_calls':3,'production_new_calls':25,'probe_columns':24,'selected_anchor_repeated':True}))
