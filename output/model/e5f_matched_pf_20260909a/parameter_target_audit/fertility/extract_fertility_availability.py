"""Bounded read-only empirical inventory/extraction; no target activation."""
from pathlib import Path
import csv,hashlib,json,collections,datetime
OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[4]
RAW=Path('/Users/tommasodesanto/Desktop/Projects/Datasets/CPS/extract3/cps_00003.dat')
LOADER=RAW.with_name('loader.do')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
width=261;n=RAW.stat().st_size//width
assert RAW.stat().st_size%width==0
records={};parts={}
with RAW.open('rb') as f:
 def key(i):
  f.seek(i*width);r=f.read(11);return int(r[:4]),int(r[9:11])
 def lower(k):
  a,b=0,n
  while a<b:
   m=(a+b)//2
   if key(m)<k:a=m+1
   else:b=m
  return a
 samplekeys=[key((n-1)*j//100) for j in range(101)]
 assert samplekeys==sorted(samplekeys)
 for yr in [2004,2006,2008]:
  a,b=lower((yr,6)),lower((yr,7))
  assert a<b and key(a)==(yr,6) and key(b-1)==(yr,6)
  assert a==0 or key(a-1)<(yr,6)
  assert b==n or key(b)>(yr,6)
  f.seek(a*width);rows=[];h=hashlib.sha256();codes=collections.Counter()
  for j in range(b-a):
   r=f.read(width);h.update(r)
   assert len(r)==width and r.endswith(b'\n')
   assert (int(r[:4]),int(r[9:11]))==(yr,6)
   if int(r[148:149])==2 and 40<=int(r[146:148])<=44:
    c=int(r[238:241]);w=float(r[250:260])/10000.;codes[c]+=1
    assert c in range(21) or c==999
    if c!=999 and w>0:rows.append((c,w))
  records[str(yr)]=rows
  parts[str(yr)]={'byte_start':a*width,'byte_end_exclusive':b*width,'all_sample_rows':b-a,'partition_sha256':h.hexdigest(),'female40_44_frever_codes':dict(codes)}
def moments(label,rows):
 sw=sum(w for c,w in rows);parents=sum(w for c,w in rows if c>0);ones=sum(w for c,w in rows if c==1)
 return {'window':label,'n':len(rows),'sum_supplement_weights':sw,'mean_children_ever_born_uncapped':sum(c*w for c,w in rows)/sw,'mean_children_ever_born_capped5':sum(min(c,5)*w for c,w in rows)/sw,'childless_share':sum(w for c,w in rows if c==0)/sw,'exactly_one_given_mother':ones/parents,'weighted_mothers_denominator':parents,'weighted_exactly_one_numerator':ones,'standard_errors':'NOT COMPUTED; not production target or weight'}
rows=[moments(y,v) for y,v in records.items()]
rows.append(moments('2004+2006 pooled weighted records',records['2004']+records['2006']))
cache=ROOT/'code/data/nchs_natality_timing/first_birth_counts_year_age.csv'
d=list(csv.DictReader(cache.open()));timing=[]
for lo,hi in [(2003,2006),(2011,2014),(2020,2023),(2003,2007),(2011,2015),(2019,2023)]:
 rs=[r for r in d if lo<=int(r['year'])<=hi];assert sorted(set(int(r['year']) for r in rs))==list(range(lo,hi+1))
 total=sum(int(r['n_first_births']) for r in rs)
 def mid(a):return 20 if a<22 else 44 if a>=42 else 20+4*((a-18)//4)
 timing.append({'start_year':lo,'end_year':hi,'first_birth_count':total,'mean_age_boundary_collapsed_midpoint':sum(int(r['n_first_births'])*mid(int(r['age'])) for r in rs)/total,'share_first_births_age30plus':sum(int(r['n_first_births']) for r in rs if int(r['age'])>=30)/total,'operator':'pooled period first-birth counts, ages12-49, no female exposure denominator','uncertainty':'not computed; existing cohort-window SE cannot be reused automatically'})
result={'status':'availability and preliminary extraction; NOT adopted targets','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cps_source':str(RAW),'raw_size_bytes':RAW.stat().st_size,'raw_mtime_ns':RAW.stat().st_mtime_ns,'raw_hash_scope':'exact selected June partitions hashed; whole6GB source not hashed in bounded audit','schema_source':str(LOADER),'schema_sha256':sha(LOADER),'schema_1based_inclusive':{'year':[1,4],'month':[10,11],'age':[147,148],'sex':[149,149],'frever':[239,241],'frsuppwt':[251,260]},'sample':'June, femaleSEX2,age40-44,FREVER0-20,FRSUPPWT>0;excludeFREVER999NIU','weights':'IPUMS fertility supplement FRSUPPWT divided10000 as supplied loader; official Census supplement weight crosswalk not independently certified','topcode':'cap FREVER at5 to match2024 public-file target convention; uncapped values retained separately','partitions':parts,'cps_moments':rows,'nchs_source':str(cache),'nchs_sha256':sha(cache),'nchs_period_timing':timing,'cautions':['Pre2007 completed-cohort fertility differs from periodTFR/replacement2.1; no automatic simultaneous equality in steady-state mapping','Exactly-one conditional on motherhood is new parity shape variation, not a continuation-birth hazard','2020-23 includes COVID period','No bootstrap or target/weight revision in this audit']}
(OUT/'fertility_availability.json').write_text(json.dumps(result,indent=2)+'\n')
for name,table in [('cps_initial_window_candidates.csv',rows),('nchs_period_timing_candidates.csv',timing)]:
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
print(json.dumps({'output':str(OUT),'cps':rows,'nchs':timing},indent=2))
