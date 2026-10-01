import hashlib,json,pathlib,shlex,subprocess
HERE=pathlib.Path(__file__).resolve().parent
REMOTE='/scratch/td2248/projects/transition_readiness_v1/current_floor/results/actualfit_v2_candidate2_report'
sha=lambda path:hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()
code='import hashlib,json,pathlib;r=pathlib.Path('+repr(REMOTE)+');print(json.dumps({n:hashlib.sha256((r/n).read_bytes()).hexdigest() for n in ("report/report.json","launcher_start.json","launcher_terminal.json")}))'
result=subprocess.run(['ssh','-4','-o','BatchMode=yes','torch','python -c '+shlex.quote(code)],capture_output=True,text=True,check=True)
manifest=json.loads(result.stdout)
def copy(rel,digest):
 path=HERE/rel;path.parent.mkdir(parents=True,exist_ok=True)
 subprocess.run(['scp','-4','-o','BatchMode=yes','-q','torch:'+REMOTE+'/'+rel,str(path)],check=True)
 if sha(path)!=digest:raise RuntimeError('Copied hash mismatch: '+rel)
 return path
for rel,digest in manifest.items():copy(rel,digest)
report=json.loads((HERE/'report/report.json').read_text())
assert report['native_model_calls']==0 and report['diagnostic_only'] and not report['fit_certified'] and not report['production_certified']
assert [r['period'] for r in report['retained_dates']]==[0,8,15]
terminal=json.loads((HERE/'launcher_terminal.json').read_text());assert terminal['exit_code']==0
names=None;selected=[]
for date in report['retained_dates']:
 assert len(date['plots'])==17
 if names is None:names=set(date['plots'])
 assert set(date['plots'])==names
 for name,pin in date['plots'].items():
  rel='report/diagnostics/date_%03d/standard_diagnostics/%s'%(date['period'],name)
  manifest[rel]=pin['sha256'];copy(rel,pin['sha256'])
for period,stem in ((0,'housing_market'),(8,'market_clearing_residuals'),(15,'fertility_by_age')):
 matches=[name for name in names if stem in name]
 assert len(matches)==1,(stem,matches)
 selected.append(str(HERE/('report/diagnostics/date_%03d/standard_diagnostics/%s'%(period,matches[0]))))
receipt=dict(status='saved_native_diagnostics_collected',plots=51,retained_periods=[0,8,15],native_model_calls=0,
 diagnostic_only=True,fit_certified=False,production_certified=False,visual_review_pending=True,
 report_sha256=manifest['report/report.json'],launcher_exit_code=terminal['exit_code'],
 sha256_by_relative_path=manifest,selected_png_paths=selected)
(HERE/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ('plots','report_sha256','selected_png_paths')},indent=2))
