import importlib.util,json,time
from pathlib import Path
from experiments.ces_normalized_shares import adapter
from production import equilibrium,native_price,native_phase_b
root=Path('/work/deployment'); packet=json.loads((root/'followup_tools/historical_start_point.json').read_text())
spec=importlib.util.spec_from_file_location('diagnostic_packet',root/'followup_tools/historical_packet.py');diagnostic=importlib.util.module_from_spec(spec);spec.loader.exec_module(diagnostic)
out=Path('/work/results/run');out.mkdir(parents=True,exist_ok=False)
price=6.318956013969669*1.6;end=time.time()+900
P,grid=adapter.load_inputs(packet['parameters'])
with adapter.install():
 context=equilibrium.build_context(P,grid,out,price_start=price,deadline=end,max_lifecycle=2,closure='population_one')
 budget=equilibrium.Budget(out,end,2)
 label='selected_diagnostic_upper'
 live=native_price.solve_fixed_price(context,float(P.unsecured_credit_limit),price,budget,label,out/'stage')
 observed=native_phase_b._observe_with_deadline(context,live,label,final=False)
 plots=diagnostic.standard_packet(context,live,out/'standard_diagnostics')
 result=dict(status='fixed_price_diagnostic_only',price=price,observed=observed,parameters=packet['parameters'],source=packet['source'],source_parameter_csv_sha256=packet['parameter_csv_sha256'],plot_names=plots,no_GE_certification=True,no_search=True,effective_input_fingerprint=adapter.canonical.effective_input_fingerprint(P,grid))
 (out/'completed.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(dict(status=result['status'],price=price,renewal_residual=observed['renewal_residual'])))
