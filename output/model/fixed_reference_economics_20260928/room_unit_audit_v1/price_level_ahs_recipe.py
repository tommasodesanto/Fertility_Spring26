import csv, io, json, math, zipfile, hashlib
from pathlib import Path
P=Path('/tmp/fertility_ahs2007_supply_check')
def num(v):
 try:return float(v.strip("' "))
 except:return None
fields=['STATUS','HHAGE','ROOMS','WGT90GEO','TENURE','RENT','FRENT','VALUE','PROJ','SUBRNT','VCHER','RCNTRL','ZINC','ZINC2']
groups={}
def add(name,w,rooms,cost,annual=True):
 t=groups.setdefault(name,{'n':0,'weight':0.,'rooms':0.,'cost':0.,'cost_per_room_sum':0.})
 t['n']+=1;t['weight']+=w;t['rooms']+=w*rooms;t['cost']+=w*cost;t['cost_per_room_sum']+=w*cost/rooms
sN=sH=0.;n=0;excluded={}
with zipfile.ZipFile(P/'ahs2007_flat.zip') as z:
 with z.open('ahs2007n.csv') as f:
  reader=csv.reader(io.TextIOWrapper(f,encoding='utf-8-sig'));h=next(reader);ix={k:h.index(k) for k in fields}
  for row in reader:
   x={k:row[i] for k,i in ix.items()}
   if x['STATUS'].strip("'")!='1':continue
   age=num(x['HHAGE']);rooms=num(x['ROOMS']);w=num(x['WGT90GEO'])
   if age is None or not 18<=age<=85 or rooms is None or rooms<=0 or w is None or w<=0:continue
   sN+=w;sH+=w*rooms;n+=1
   if x['TENURE'].strip("'")=='2':
    rent=num(x['RENT']);freq=num(x['FRENT'])
    if rent is None or rent<=1 or freq is None or not 0<freq<=52:
     reason='income_dependent' if rent==1 else 'invalid_or_topcoded_frequency'
     excluded[reason]=excluded.get(reason,0)+1;continue
    add('renter_numeric_contract',w,rooms,rent*freq)
    if not any(x[k].strip("'")=='1' for k in ['PROJ','SUBRNT','VCHER','RCNTRL']): add('renter_no_reported_subsidy_or_control',w,rooms,rent*freq)
    if freq==12:add('renter_monthly_only',w,rooms,rent*12)
   elif x['TENURE'].strip("'")=='1':
    value=num(x['VALUE'])
    if value is not None and value>0:add('owner_positive_value',w,rooms,value,False)
for t in groups.values():
 t['mean_rooms']=t['rooms']/t['weight'];t['mean_cost']=t['cost']/t['weight'];t['room_weighted_price']=t['cost']/t['rooms'];t['household_weighted_price']=t['cost_per_room_sum']/t['weight']
out={'sample':'STATUS=1,HHAGE18--85,ROOMS>0,WGT90GEO>0','n':n,'weighted_households':sN,'weighted_rooms':sH,'mean_rooms':sH/sN,'groups':groups,'excluded_rent_records':excluded,'time':'rent*FRENT annually, frequency topcode53 excluded; RENT=1 income-dependent excluded','quantity':'occupied housing only, literal public-use rooms topcode21','raw_sha256':hashlib.sha256((P/'ahs2007_flat.zip').read_bytes()).hexdigest()}
(P/'ahs_price_quantities.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
