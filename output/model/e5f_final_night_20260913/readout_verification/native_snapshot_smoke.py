from pathlib import Path
import json,importlib.util,sys,time,hashlib
b=Path(sys.argv[1]);label=sys.argv[2];helper=sys.argv[3];case=sys.argv[4]
r=b/'readout_source/collect_e5f_final_history_readout.py';profile=b/'readout_source/collect_e5f_patch_readout.py'
spec=importlib.util.spec_from_file_location('readout_smoke',r);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
manifest_path=b/f'history_manifest_{label}.json';sm=json.loads(manifest_path.read_text());plan=json.loads(Path(sm['prior_plan']).read_text());source=Path(plan['source_root']);pins={Path(p).resolve():h for p,h in sm['file_sha256'].items()};pins.update({Path(p).resolve():h for p,h in plan['file_sha256'].items()});pins[profile]=m.file_sha256(profile);pins[r]=m.file_sha256(r);m.verify_pins(pins)
folder=b/f'histories_{label}_6'/case/'window_2007/trial_00';pack_path=folder/'first_period_diagnostics.pkl.gz';accepted_path=folder/'accepted_forecast.pkl.gz';root_path=folder/'root_receipt.json'
for path in [pack_path,accepted_path,root_path]:pins[path]=m.file_sha256(path)
payload={'runtime_paths':[str(b/helper),str(profile.parent),str(source/'code/model/tools'),str(source/'code/model')],'profile_kernel':str(profile)}
started=time.monotonic();rt=m.configure_runtime(source_root=source,pins=pins,kernels=[p for p in pins if p.suffix=='.py'],manifest=payload,manifest_path=manifest_path)
pack=m.validate_pack(m.load_gzip_pickle(pack_path),'2007 saved native snapshot');accepted=m.load_gzip_pickle(accepted_path);result=m.validate_accepted(accepted,m.read_json(root_path));row=m.row_for_year(result.path.rows,2007)
check=m.snapshot_aggregate_check(pack,row,rt['profile'],rt['fertility'].period_fertility_diagnostics)
e,P,grid,shared=[pack[k] for k in ['evaluation','parameters','b_grid','shared']]
values={}
values['fertility_stock_timing']=rt['fertility_stock'].observe_initial_fertility(e,P,age_projection='uniform_birth_time')
values['housing_wealth']=rt['housing_wealth'].observe_initial_housing_wealth(e,P,grid,shared,diagnostic_enabled=True,age_projection='uniform_within_age_cell',diagnostic_allow_family_proxies=True,include_wealth=True,include_birth_response=False)
recent=rt['recent_parent'];values['recent_parent']=recent.observe_recent_parent_flow(e,P,diagnostic_enabled=True,snapshot=recent.SNAPSHOT,age_projection=recent.AGE_PROJECTION,diagnostic_allow_residence_proxy=True)
m.verify_pins(pins)
out=b/'readout_native_smokes';out.mkdir(exist_ok=True)
m.write_json(out/f'{label}.json',dict(status='PASS',calendar_year=2007,source_label=label,seconds=time.monotonic()-started,aggregate_check=check,observations=values,model_solve_performed=False,scope='Single saved2007snapshot and accepted forecast validation. The final2019to2023dated birth-room observer remains pending.',source_files={str(p):h for p,h in pins.items()}))
print(json.dumps(dict(status='PASS',calendar_year=2007,label=label,maximum_abs=check['maximum_abs'],seconds=time.monotonic()-started)),flush=True)
