"""Saved CPS June 2004/2006 teen stock diagnostic; Torch only, zero model solves."""
import hashlib
import json
import math
import os
from pathlib import Path

SOURCE = Path('/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/early_fertility_target_20260926')
OUT = Path('/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project/output/model/fertility_identification_20260928/measurement_audit_v1/age_tail_v1')

def main():
    assert os.environ.get('SLURM_JOB_ID','').isdigit() and os.uname().sysname == 'Linux'
    receipt=json.loads((SOURCE/'fertility_availability.json').read_text())
    result={}
    hashes={}
    for year in (2004,2006):
        path=SOURCE/'input'/f'june_{year}.dat'
        data=path.read_bytes()
        digest=hashlib.sha256(data).hexdigest()
        assert digest==receipt['partitions'][str(year)]['partition_sha256']
        hashes[str(path)]=digest
        assert len(data)%261==0
        for off in range(0,len(data),261):
            row=data[off:off+261]
            assert int(row[:4])==year and int(row[9:11])==6 and row.endswith(b'\n')
            if row[148:149]!=b'2':continue
            age=int(row[146:148])
            if age not in (17,18,19,25):continue
            n=int(row[238:241]); w=int(row[250:260])/10000
            x=result.setdefault(str(age),{'n':0,'invalid_frever':0,'nonpositive_weight':0,'weight':0.,'weighted_capped3':0.,'weighted_any':0.})
            if n==999:x['invalid_frever']+=1;continue
            assert 0<=n<=20
            if w<=0:x['nonpositive_weight']+=1;continue
            x['n']+=1;x['weight']+=w;x['weighted_capped3']+=w*min(n,3);x['weighted_any']+=w*(n>0)
    for x in result.values():
        x['children_capped3']=x.pop('weighted_capped3')/x['weight']
        x['mother_share']=x.pop('weighted_any')/x['weight']
    assert result['25']['n']==1774 and abs(result['25']['children_capped3']-.8095276384290021)<1e-12
    payload={'status':'PASS','job_id':os.environ['SLURM_JOB_ID'],'model_solves':0,'model_imports':0,'checkpoint_reads':0,'sample':'June 2004/2006 CPS women, positive FRSUPPWT, FREVER 0-20; exact completed interview ages. Age25 control reproduces approved builder. Teen ages are diagnostic only, not longitudinally linked to age25 sample.','ages':result,'input_sha256':hashes,'limits':['Age17 represents [17,18), not stock exactly at 18th birthday.','The cross-sectional age17 stock cannot be subtracted directly from age25 mean to identify a cohort contribution.']}
    (OUT/'teen_entry_stock.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({'ages':result,'status':'PASS'}))

if __name__=='__main__':main()
