"""Candidate pre-2007 fertility age profiles; no objective or uncertainty adoption."""
from pathlib import Path
import csv, hashlib, json

RAW=Path('/Users/tommasodesanto/Desktop/Projects/Datasets/CPS/extract3/cps_00003.dat'); W=261
receipt=Path(__file__).parents[3]/'parameter_target_audit'/'fertility'/'fertility_availability.json'
old=json.loads(receipt.read_text()); out=Path(__file__).parent; out.mkdir(parents=True,exist_ok=True)
assert hashlib.sha256(Path(old['schema_source']).read_bytes()).hexdigest()==old['schema_sha256']
bands=[(20,24),(25,29),(30,34),(35,39),(40,44)]; rows=[]; allrows=[]; excluded=[]
with RAW.open('rb') as f:
 for yr in (2004,2006):
  p=old['partitions'][str(yr)]; f.seek(p['byte_start']); raw=f.read(p['byte_end_exclusive']-p['byte_start'])
  assert hashlib.sha256(raw).hexdigest()==p['partition_sha256']
  assert len(raw)==W*p['all_sample_rows']
  for i in range(0,len(raw),W):
   r=raw[i:i+W]; age=int(r[146:148]); sex=int(r[148:149]); c=int(r[238:241]); wt=float(r[250:260])/10000
   assert len(r)==W and r.endswith(b'\n') and (int(r[:4]),int(r[9:11]))==(yr,6)
   if sex==2 and 20<=age<=44:
    assert c in range(21) or c==999
    if c==999 and wt>0: excluded.append((yr,age))
   if sex==2 and 20<=age<=44 and c in range(21) and wt>0: allrows.append((yr,age,c,wt))
for yr in (2004,2006):
 for lo,hi in bands:
  x=[z for z in allrows if z[0]==yr and lo<=z[1]<=hi]; rows.append((yr,lo,hi,x))
for lo,hi in bands: rows.append(('pooled',lo,hi,[z for z in allrows if lo<=z[1]<=hi]))
def m(k,lo,hi,x):
 s=sum(z[3] for z in x); return {'window':k,'age_lower':lo,'age_upper':hi,'n':len(x),'total_weight':s,'mean_CEB_uncapped':sum(z[2]*z[3] for z in x)/s,'mean_CEB_capped5':sum(min(z[2],5)*z[3] for z in x)/s,'mean_model_coded_CEB':sum((0 if z[2]==0 else 1 if z[2]==1 else 2 if z[2]==2 else 3.602359422009)*z[3] for z in x)/s,'share_0':sum(z[3] for z in x if z[2]==0)/s,'share_1':sum(z[3] for z in x if z[2]==1)/s,'share_2':sum(z[3] for z in x if z[2]==2)/s,'share_3plus':sum(z[3] for z in x if z[2]>=3)/s,'exactly_one_given_mothers':sum(z[3] for z in x if z[2]==1)/sum(z[3] for z in x if z[2]>0),'frever999_excluded_n':0}
result=[m(k,lo,hi,x) for k,lo,hi,x in rows]
checks=[]
for row in result:
 del row['frever999_excluded_n']
 row['frever999_excluded_positive_weight_n']=sum(1 for y,a in excluded if row['age_lower']<=a<=row['age_upper'] and (row['window']=='pooled' or y==row['window']))
 assert abs(sum(row['share_'+n] for n in ('0','1','2','3plus'))-1)<1e-12
 assert abs(row['mean_model_coded_CEB']-(row['share_1']+2*row['share_2']+3.602359422009*row['share_3plus']))<1e-12
 if row['age_lower']==40:
  label=str(row['window']) if row['window']!='pooled' else '2004+2006 pooled weighted records'
  anchor=next(r for r in old['cps_moments'] if r['window']==label)
  gaps={key:row[key]-anchor[oldkey] for key,oldkey in {'mean_CEB_uncapped':'mean_children_ever_born_uncapped','mean_CEB_capped5':'mean_children_ever_born_capped5','share_0':'childless_share','exactly_one_given_mothers':'exactly_one_given_mother'}.items()}
  assert max(abs(v) for v in gaps.values())<1e-12
  checks.append(dict(window=row['window'],gaps=gaps))
with (out/'age_profile_candidates.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=result[0],lineterminator='\n');writer.writeheader();writer.writerows(result)
metadata={'status':'candidate only; no target activation or weights','source':str(RAW),'schema_source':str(old['schema_source']),'schema_sha256':old['schema_sha256'],'partition_receipt':{y:old['partitions'][y] for y in ('2004','2006')},'prior_receipt_sha256':hashlib.sha256(receipt.read_bytes()).hexdigest(),'builder_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sample':'June2004/2006,femaleSEX2,statedageband,FREVER0–20,FRSUPPWT>0','weights':'FRSUPPWT/10000; pooled observed supplement weights','uncertainty':'Not computed; no empirical SE or synthetic scale adopted','fixed_effects':None,'clustering':'Point estimates only; uncertainty not computed','top_bin_representative':3.602359422009,'top_bin_source':'code/data/cps_fertility/output/cps_fertility_targets.csv:capped_top_bin_mean','measurement':'Raw/capped5 empirical CEB and fixed-representative model-coded CEB are distinct observations','existing_40_44_anchor_checks':checks,'rows':result}
(out/'age_profile_candidates.json').write_text(json.dumps(metadata,indent=2)+'\n')
print(json.dumps(dict(status='PASS',rows=len(result),existing_anchor_max_error=max(abs(v) for x in checks for v in x['gaps'].values()),positive_weight_NIU_exclusions=len(excluded))))
