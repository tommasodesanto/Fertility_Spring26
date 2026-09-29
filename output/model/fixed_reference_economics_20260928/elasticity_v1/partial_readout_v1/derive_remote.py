"""Read-only Torch derivation from authenticated saved receipts and fit tables."""
import csv,hashlib,json,math,sys
from pathlib import Path
mode=sys.argv[1]
root=Path('/scratch/td2248/projects/fixed_reference_elasticity_v2_20260929/results_v2')
old=root/'solve_v2'; continuation=root/'solve_v4'
label='2007 stationary reference — block0506, September 28 verified export'
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def receipt(case,low):
    path=(continuation if low else old)/case/'receipt.json'
    return json.loads(path.read_text()),sha(path)
if mode=='slopes':
    fields=('reference_label','regime','scope','outcome','unit','factor_low','factor_high',
            'value_low','value_high','one_sided_log_elasticity','low_case','high_case',
            'low_receipt_sha256','high_receipt_sha256','plan_sha256')
    writer=csv.DictWriter(sys.stdout,fieldnames=fields);writer.writeheader()
    for regime,center in (('reference','grid_control'),('credit','credit')):
        lo,lo_pin=receipt(regime+'_990',True);hi,hi_pin=receipt(center,False)
        assert lo['plan_sha256']==hi['plan_sha256']=='6354bafd069fbfc4eda21271fc91747a22e3aeb8c90a57293f7dfde01badf28f'
        for scope in ('impact','cohort'):
            for name in ('births_per_household','first_births','second_births')+(
                    ('completed_fertility',) if scope=='cohort' else ()):
                def value(r):
                    if name=='completed_fertility':return float(r[name])
                    summary=r[scope+'_summary'];v=float(summary[name])
                    return v/float(summary['household_mass']) if name in ('first_births','second_births') else v
                a,b=value(lo),value(hi)
                assert a>0 and b>0
                writer.writerow(dict(reference_label=label,regime=regime,scope=scope,outcome=name,
                    unit='children' if name=='completed_fertility' else 'births per household',
                    factor_low=.99,factor_high=1.,value_low=a,value_high=b,
                    one_sided_log_elasticity=(math.log(b)-math.log(a))/(math.log(1.)-math.log(.99)),
                    low_case=regime+'_990',high_case=center,low_receipt_sha256=lo_pin,
                    high_receipt_sha256=hi_pin,plan_sha256=lo['plan_sha256']))
elif mode in ('refinement','parameters'):
    filename='target_fit.csv' if mode=='refinement' else 'parameters.csv'
    old_path=Path('/scratch/td2248/projects/fixed_reference_credit_20260929/results/solve_v1/credit')/filename
    union_path=old/'credit'/filename
    previous=list(csv.DictReader(old_path.open()));current=list(csv.DictReader(union_path.open()))
    assert len(previous)==len(current)==(14 if mode=='refinement' else 31)
    if mode=='refinement':
        fields=('reference_label','moment','target','model_262','model_602','difference_602_minus_262',
                'old_table_sha256','union_table_sha256')
        writer=csv.DictWriter(sys.stdout,fieldnames=fields);writer.writeheader()
        for a,b in zip(previous,current):
            assert a['moment']==b['moment']
            writer.writerow(dict(reference_label=label,moment=a['moment'],target=a['target'],
                model_262=a['model'],model_602=b['model'],
                difference_602_minus_262=float(b['model'])-float(a['model']),
                old_table_sha256=sha(old_path),union_table_sha256=sha(union_path)))
    else:
        fields=('reference_label','parameter','estimate_262','estimate_602','difference_602_minus_262',
                'old_table_sha256','union_table_sha256')
        writer=csv.DictWriter(sys.stdout,fieldnames=fields);writer.writeheader()
        for a,b in zip(previous,current):
            assert a['parameter']==b['parameter']
            writer.writerow(dict(reference_label=label,parameter=a['parameter'],
                estimate_262=a['estimate'],estimate_602=b['estimate'],
                difference_602_minus_262=float(b['estimate'])-float(a['estimate']),
                old_table_sha256=sha(old_path),union_table_sha256=sha(union_path)))
else:raise SystemExit('Unknown extraction mode')
