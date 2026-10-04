"""Collect only complete, internally exact CES chains; never adopt a result."""
import argparse,csv,hashlib,json,sys
from pathlib import Path
def read(p):return json.loads(Path(p).read_text())
def rows(p,expected):
 with Path(p).open(newline='') as f:r=list(csv.DictReader(f))
 if len(r)!=expected:raise RuntimeError(f"{p} has {len(r)}, expected {expected}")
 return r
def report_path(raw,launch):
 # Container /work/results maps to the host chain launch directory; receipts include /run/.
 prefix='/work/results'
 if not isinstance(raw,str) or not (raw==prefix or raw.startswith(prefix+'/')):raise RuntimeError('report path is outside container /work/results root')
 rel=Path(raw[len(prefix):].lstrip('/'))
 if rel.is_absolute() or '..' in rel.parts:raise RuntimeError('report path contains traversal')
 launch=Path(launch).resolve(); report=(launch/rel).resolve()
 if not report.is_relative_to(launch):raise RuntimeError('report path escapes chain launch directory')
 return report
def verified_reports(child,launch):
 report=report_path(child['selected_root']['report'],launch)
 repeat=report_path(child['selected_repeat_final']['report'],launch)
 fits=rows(report/'target_fit_experimental.csv',14); params=rows(report/'parameters.csv',31)
 hashes=child.get('standard_plot_hashes',{})
 plot_dir=report/'standard_diagnostics'
 if len(hashes)!=17:raise RuntimeError('postcheck must identify 17 standard plots')
 for name,digest in hashes.items():
  if Path(name).name!=name or not isinstance(digest,str):raise RuntimeError('invalid standard plot receipt entry')
  plot=plot_dir/name
  if not plot.is_file():raise RuntimeError('missing standard plot: '+str(plot))
  actual=hashlib.sha256(plot.read_bytes()).hexdigest()
  if actual!=digest:raise RuntimeError('standard plot hash mismatch: '+name)
 return report,repeat,fits,params
def markdown(chain,child,launch):
 report,repeat,fits,params=verified_reports(child,launch)
 out=[f"## Chain {chain}","",f"Postcheck: `{child['status']}`. Root report: `{report}`. Internal repeat: `{repeat}`.","","### Full target fit","","| Moment | Target | Model | Gap | Weight | Loss contribution |","|---|---:|---:|---:|---:|---:|"]
 out += [f"| {r['moment']} | {r['target']} | {r['model']} | {r['gap']} | {r['weight']} | {r['loss_contribution']} |" for r in fits]
 out += ["","### Parameters","","| Parameter | Estimate | Lower | Upper | Status | Near bound |","|---|---:|---:|---:|---|---|"]
 out += [f"| {r['parameter']} | {r['estimate']} | {r['lower']} | {r['upper']} | {r.get('status','')} | {r.get('near_bound','')} |" for r in params]
 return '\n'.join(out)
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--root",type=Path,required=True);ap.add_argument("--mode",choices=("smoke","production"),required=True);a=ap.parse_args();inv=read(a.root/"inventory.json");plan=read(a.root/'source/output/model/experiments/ces_normalized_shares/overnight_v1/start_plan.json');accepted=[];reject=[]
 for i in range(4 if a.mode=="production" else 1):
  run=a.root/"results"/f"{a.mode}_chain_{i}"/"run"
  needed=(run/'completed.json',run/'start_contract.json',run/'native_postcheck/completed.json')
  if not all(p.is_file() for p in needed):reject.append(dict(chain=i,reason='pending_or_incomplete'));continue
  c,s,child=map(read,needed)
  keys=('target_fingerprint','weight_fingerprint','starts_file_sha256','selected_source_sha256','source_checkpoint_sha256')
  if c.get('status')!='selected_numerically_verified' or child.get('status')!='full_native_postcheck_passed' or any(s.get(k)!=(inv['start_plan_sha256'] if k=='starts_file_sha256' else plan.get(k)) or child.get(k)!=(inv['start_plan_sha256'] if k=='starts_file_sha256' else plan.get(k)) for k in keys) or s.get('target_fingerprint')!=inv['target_fingerprint'] or s.get('weight_fingerprint')!=inv['weight_fingerprint'] or s.get('starts_file_sha256')!=inv['start_plan_sha256']:
   reject.append(dict(chain=i,reason='receipt_identity_or_completion_mismatch'));continue
  try:
   if child.get('target_rows')!=14 or child.get('parameter_rows')!=31 or len(child.get('standard_plot_hashes',{}))!=17:raise RuntimeError('postcheck shape')
   markdown(i,child,run.parent)
  except Exception as e:reject.append(dict(chain=i,reason='missing_report_artifact',detail=str(e)));continue
  accepted.append(dict(chain=i,status=c['status'],completed=c,postcheck=child))
 status='rejected_incomplete_or_mixed_chain' if reject else 'collected_no_adoption'
 receipt=dict(mode=a.mode,accepted=accepted,rejected=reject,status=status,no_adoption=True)
 (a.root/f"{a.mode}_collection.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
 text=['# CES normalized-share collection','',f'Status: **{status}**. This collector neither adopts nor publishes a calibration.']
 if reject:text += ['', '## Rejected or incomplete chains','',* [f"- Chain {x['chain']}: {x['reason']}" for x in reject]]
 for x in accepted:text += ['',markdown(x['chain'],x['postcheck'],a.root/'results'/f"{a.mode}_chain_{x['chain']}")]
 (a.root/'RESULTS.md').write_text('\n'.join(text)+'\n')
 print(json.dumps(dict(status=status,accepted=len(accepted),rejected=len(reject))))
 if reject:raise SystemExit(1)
if __name__=="__main__":main()
