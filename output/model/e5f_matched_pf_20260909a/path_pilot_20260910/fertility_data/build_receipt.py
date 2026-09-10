"""Reproduce the bounded official published-series extraction; no model imports."""
from pathlib import Path
import csv, hashlib, json, re, sys
from decimal import Decimal
from pypdf import PdfReader
ROOT = Path(__file__).resolve().parent
SRC=ROOT/'sources'
# Optional local scratch install; otherwise use xlrd from the active environment.
sys.path.insert(0,'/tmp/e5f_fertility_xlrd_20260910')
import xlrd
sources={
 'nchs_2015_final':{'file':'nvsr66-1_2015_final.pdf','url':'https://www.cdc.gov/nchs/data/nvsr/nvsr66/NVSR66_01.pdf','publication_date':'2017-01-05','tfr_table':4,'tfr_pdf_page':20,'births_table':1,'births_pdf_page':16},
 'nchs_2023_final':{'file':'nvsr74-1_2023_final.pdf','url':'https://www.cdc.gov/nchs/data/nvsr/nvsr74/nvsr74-1.pdf','publication_date':'2025-03-18','tfr_table':2,'tfr_pdf_page':14,'births_table':1,'births_pdf_page':12},
 'census_hh3':{'file':'census_hh3_download_20260910.xls','url':'https://www2.census.gov/programs-surveys/demo/tables/families/time-series/households/hh3.xls','publication_date':'2025-12','table':'HH-3','download_date':'2026-09-10'},
}
def save_csv(name,rows):
 with (ROOT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
parsed={}
for key in ['nchs_2015_final','nchs_2023_final']:
 s=sources[key];r=PdfReader(SRC/s['file']);parsed[key]={'tfr':{},'births':{}}
 for kind,page in [('tfr',s['tfr_pdf_page']),('births',s['births_pdf_page'])]:
  txt=r.pages[page-1].extract_text()
  for line in txt.splitlines():
   m=re.match(r'^(20\d\d)[\s.]+([\d,]+\.\d|[\d,]+)(?:\s|$)(.*)',line)
   if not m:continue
   y=int(m[1])
   if y in parsed[key][kind]:continue # all-race group comes first
   n=Decimal(m[2].replace(',',''))
   parsed[key][kind][y]={'value':n,'source_line':line}
   if kind=='tfr':
    rest=[Decimal(x) for x in m[3].split()]
    assert len(rest)==10,(y,rest)
    asfr=[rest[i] for i in [0,1,4,5,6,7,8,9]]
    assert sum(asfr)*5==n,(key,y,n,asfr)
    parsed[key][kind][y]['asfr']=asfr
checks=[]
for y in range(2010,2016):
 for k in ['tfr','births']:
  a,b=[parsed[s][k][y]['value'] for s in ['nchs_2015_final','nchs_2023_final']]
  assert a==b;checks.append({'year':y,'series':k,'earlier_report':str(a),'later_report':str(b),'exact_match':True})
rows=[];age_rows=[]
for y in range(2007,2024):
 key='nchs_2015_final' if y<2010 else 'nchs_2023_final';s=sources[key];p=parsed[key]
 rows.append({'year':y,'period_tfr_births_per_woman':str(p['tfr'][y]['value']/1000),'published_tfr_per_1000_women':str(p['tfr'][y]['value']),'live_births':int(p['births'][y]['value']),'source_id':key,'tfr_table':s['tfr_table'],'tfr_pdf_page':s['tfr_pdf_page'],'births_table':1,'births_pdf_page':s['births_pdf_page'],'status':'verified_published_final','source_url':s['url']})
 for group,val in zip(['10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49'],p['tfr'][y]['asfr']):
  age_rows.append({'year':y,'female_denominator_age_group':group,'births_per_1000_females':str(val),'source_id':key,'female_exposure_count':'','exposure_status':'not_extracted_do_not_reverse_engineer_rounded_rates'})
save_csv('annual_fertility_2007_2023.csv',rows);save_csv('annual_published_age_specific_rates.csv',age_rows)
sh=xlrd.open_workbook(str(SRC/sources['census_hh3']['file'])).sheet_by_index(0)
hh={};alternatives=[]
for i in range(10,31):
 label=str(sh.cell_value(i,0));m=re.match(r'^(20\d\d)',label)
 if not m:continue
 y=int(m[1]);n=sh.cell_value(i,1)
 if not 2007<=y<=2023:continue
 rec={'year':y,'published_row_label':label,'households_thousands':n,'households_stock_proxy':int(n*1000),'xlsx_row_1_based':i+1,'source_id':'census_hh3','selection':'selected_revised_if_available' if y not in hh else 'unrevised_alternative_not_selected'}
 alternatives.append(rec)
 if y not in hh:hh[y]=rec
assert set(hh)==set(range(2007,2024))
save_csv('annual_household_stocks_2007_2023.csv',[hh[y] for y in sorted(hh)])
save_csv('household_source_rows_including_alternatives.csv',alternatives)
byyear={r['year']:r for r in rows};blocks=[]
for t in [2007,2011,2015,2019]:
 ys=list(range(t+1,t+5));assert len(ys)==4 and len(set(ys))==4 and all(y in byyear for y in ys);b=sum(byyear[y]['live_births'] for y in ys);tfr=sum(Decimal(byyear[y]['period_tfr_births_per_woman']) for y in ys)/4;h=sum(hh[y]['households_stock_proxy'] for y in ys)
 blocks.append({'decision_year':t,'birth_year_start':t+1,'birth_year_end':t+4,'live_births_total':b,'annual_observations_count':len(ys),'unique_annual_observations_count':len(set(ys)),'period_tfr_arithmetic_mean':str(tfr),'live_births_index_2008_2011':None,'household_stock_proxy_sum':h,'live_births_per_household_stock_proxy_year':b/h,'household_rate_status':'diagnostic_March_stock_proxy_all_head_ages_not_literal_exposure','tfr_aggregation':'equal_weight_mean_of_four_annual_published_TFRs_not_exposure_pooled'})
for r in blocks:r['live_births_index_2008_2011']=r['live_births_total']/blocks[0]['live_births_total']
save_csv('empirical_blocks.csv',blocks)
for key,s in sources.items():
 f=SRC/s['file'];s['sha256']=hashlib.sha256(f.read_bytes()).hexdigest();s['bytes']=f.stat().st_size
(ROOT/'source_manifest.json').write_text(json.dumps(sources,indent=2)+'\n')
(ROOT/'verification.json').write_text(json.dumps({'annual_years_complete':list(range(2007,2024)),'overlap_verification':checks,'every_tfr_equals_five_times_sum_eight_published_age_rates':True,'household_selection_note':'Source first occurrence is revised for2011 and2021; original alternatives retained. All five existing model HH-3 dates agree.','empirical_2007_tfr':2.12,'author_initial_benchmark_distinct_assumption':2.1,'no_model_comparability_certified':True},indent=2)+'\n')
print(json.dumps(blocks,indent=2))
