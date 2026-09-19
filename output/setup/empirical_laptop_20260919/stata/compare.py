import csv, hashlib, json, math
from pathlib import Path
root=Path(__file__).resolve().parents[4]
here=Path(__file__).resolve().parent
baseline=root/'code/data/psid_followup_mar2026/output/sa_rooms_first_birth_household_aligned_v1'
produced=here/'sa_rooms_first_birth_household_aligned_v1'
reports={}
for name in ('target_receipt.csv','event_study_estimates.csv'):
    old=list(csv.DictReader((baseline/name).open()))
    new=list(csv.DictReader((produced/name).open()))
    assert len(old)==len(new) and list(old[0])==list(new[0]), name
    diffs=[]; numeric=0; largest=0.
    for i,(a,b) in enumerate(zip(old,new)):
        for k in a:
            if k=='runtime_seconds':continue
            try:x,y=float(a[k]),float(b[k])
            except ValueError:
                if a[k]!=b[k]:diffs.append(dict(row=i,column=k,old=a[k],new=b[k]))
                continue
            numeric+=1
            if math.isnan(x) and math.isnan(y):continue
            gap=abs(x-y);largest=max(largest,gap)
            if not math.isclose(x,y,rel_tol=0,abs_tol=1e-9):diffs.append(dict(row=i,column=k,old=x,new=y,gap=gap))
    reports[name]=dict(rows=len(new),numeric_cells=numeric,max_abs_gap=largest,differences=diffs,passed=not diffs)
pins=json.loads((here/'baseline_hashes.json').read_text())
unchanged=all(hashlib.sha256((root/p).read_bytes()).hexdigest()==h for p,h in pins.items())
report=dict(files=reports,baseline_unchanged=unchanged,passed=unchanged and all(x['passed'] for x in reports.values()),excluded_fields=['runtime_seconds'],absolute_tolerance=1e-9)
(here/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
assert report['passed']
