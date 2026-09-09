from pathlib import Path
import json,hashlib,zipfile
from PIL import Image
R=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
O=R/'output/model/e5f_simple_fertility_tax_transition_20260908a';G=O/'graphs/v2/full'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
packets=[];nfiles=0
for case in ['tax1-equal-rebate','tax2-equal-rebate']:
 numerical=O/'results/full'/case;rh=digest(numerical/'receipt.json')
 for year in range(2023,2064,4):
  p=G/case/f'date_{year}';m=json.loads((p/'graph_manifest.json').read_text())
  assert m['status']=='complete' and m['calendar_year']==year and m['case']==case
  assert m['standard_png_count']==17 and m['numerical_solves']==0 and not m['kernel_changes']
  assert m['input_provenance']['stage_receipt_sha256']==rh
  endpoint=json.loads((numerical/f'date_{year}/endpoint_receipt.json').read_text())
  assert m['input_checkpoint_sha256']==endpoint['checkpoint_sha256']
  assert m['reporting_corrections']['input_checkpoint_sha256_unchanged']==endpoint['checkpoint_sha256']
  assert m['reporting_corrections']['maximum_native_owner_operator_gap']<2e-11
  assert max(m['reporting_corrections']['validation'].values())<2e-10
  assert all(x['formatting_preserved_line_arrays'] for x in m['reporting_corrections']['figure_checks'])
  for name,h in m['artifact_sha256'].items():assert digest(p/name)==h,(p,name)
  for name,h in m['graph_sha256'].items():
   image=p/'standard_diagnostics'/name;assert digest(image)==h
   with Image.open(image) as im:
    assert im.width>=1000 and im.height>=600;im.verify()
   nfiles+=1
  packets.append({'case':case,'year':year,'graph_manifest_sha256':digest(p/'graph_manifest.json'),'checkpoint_sha256':m['input_checkpoint_sha256']})
assert nfiles==374 and len(packets)==22
receipt={'status':'all_graph_artifacts_verified','packets':packets,'png_count':nfiles,'numerical_solves':0,'visual_review':'Smoke inspected directly; endpoint appendix and supplemental plots require PDF visual QA.'}
(O/'morning_packet/graph_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
readme='''# Complete standard graph set\n\nAnnual property tax1% versus2%, equal rebates in both paths. All11dates2023–2063;17standard graphs per date and case (374PNGs). See the morning PDF for summary and complete calibration tables.\n\nEach case/date folder contains its standard_diagnostics figures, numerical-source hash provenance, reporting-correction receipt and underlying diagnostic tables. Original numerical checkpoints remain on Torch. No model was solved by this exporter.\n\nInterpretation: rooms are occupied housing services; prices and user costs are distinct. Policy consumption/housing panels condition on renting with no birth. Fertilitypolicy curves show attempt probabilities; fertility_by_age shows first births divided by prechoice population at risk. Income indices retain permanent-group ordering. 'Young' uses model age nodes26/30/34; exact annual-age ACS alignment is outstanding.\n\nReporting corrections sum all owned products and conception outcomes and use prechoice birth-risk weights. Every correction was checked against the native operator. The maintained future household-entry closure is diagnostic; this is not a resident-population forecast.\n'''
(G/'README.md').write_text(readme)
archive=R/'output/pdf/rebated_tax_complete_graphs_20260909.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
 for p in sorted(G.rglob('*')):
  if p.is_file():z.write(p,p.relative_to(G))
 z.write(O/'morning_packet/graph_verification.json','graph_verification.json')
print(json.dumps({'verified_pngs':nfiles,'archive':str(archive),'bytes':archive.stat().st_size}))
