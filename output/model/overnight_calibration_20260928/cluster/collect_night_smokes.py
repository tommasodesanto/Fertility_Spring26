"""Collect compact smoke evidence remotely; never approve or copy checkpoints."""
import hashlib,json,pathlib,shutil,sys
from PIL import Image,ImageOps,ImageDraw
stage=pathlib.Path(sys.argv[1]);dest=pathlib.Path(sys.argv[2]);dest.mkdir(exist_ok=False)
original=pathlib.Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
smoke=stage/'project/output/model/overnight_calibration_20260928/gated_v1/smoke'
complete=json.loads((smoke/'complete.json').read_text())
for p in smoke.glob('*.json'):shutil.copy2(p,dest/p.name)
for p in smoke.glob('*.log'):
    if p.stat().st_size<=2*1024*1024:shutil.copy2(p,dest/p.name)
    else:
        with p.open('rb') as stream:stream.seek(-65536,2);tail=stream.read()
        (dest/(p.name+'.last64KiB')).write_bytes(tail)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
review={'status':complete['status'],'approval':False,'records':[],'plot_groups':{},'checkpoint_files_copied':False}
for request in sorted(smoke.glob('*.request.json')):
    req=json.loads(request.read_text());case=smoke/req['id']/'case';target=dest/req['id'];target.mkdir()
    for suffix in ('*.json','*.csv'):
        for p in case.glob(suffix):shutil.copy2(p,target/p.name)
    for p in (case.parent/'failure.json',case.parent/'checkpoint.json'):
        if p.is_file():shutil.copy2(p,target/('child_'+p.name))
    plots=sorted((case/'standard_diagnostics').glob('*.png'))
    hashes={p.name:sha(p) for p in plots};group=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()
    record=dict(case=req['id'],lane=req['lane'],plot_count=len(plots),plot_hashes=hashes,plot_group=group,
                source_case=str(case),request_sha256=sha(request),files={p.name:sha(p) for p in target.iterdir() if p.is_file()})
    review['records'].append(record)
    if len(plots)!=17:
        record['plot_status']='incomplete';continue
    summary=case/'standard_diagnostics/summary.json'
    if summary.is_file():shutil.copy2(summary,target/'standard_diagnostics_summary.json')
    if group not in review['plot_groups']:
        review['plot_groups'][group]=dict(representative=req['id'],cases=[])
        saved=target/'standard_diagnostics';saved.mkdir()
        for p in plots:shutil.copy2(p,saved/p.name)
        for page in range(3):
            canvas=Image.new('RGB',(1500,1410),'white');draw=ImageDraw.Draw(canvas)
            for i,p in enumerate(plots[page*6:(page+1)*6]):
                tile=ImageOps.contain(Image.open(p).convert('RGB'),(740,425));x=(i%2)*750;y=(i//2)*470
                draw.text((x+8,y+7),p.name,fill='black');canvas.paste(tile,(x+(750-tile.width)//2,y+35))
            canvas.save(target/f'standard_contact_{page+1}.png')
    review['plot_groups'][group]['cases'].append(req['id'])
review['all_six_plot_packets_byte_identical']=len(review['records'])==6 and len(review['plot_groups'])==1 and all(r['plot_count']==17 for r in review['records'])
(dest/'collection_review.json').write_text(json.dumps(review,indent=2)+'\n')
print(json.dumps({k:v for k,v in review.items() if k not in ('records','plot_groups')},indent=2))

import csv
pairs={}
for lane in ('primary','identity','block'):
    ids=[r['case'] for r in review['records'] if r['lane']==lane]
    checks={}
    for name,n in [('target_fit.csv',14),('parameters.csv',31)]:
        paths=[dest/i/name for i in ids]
        if len(paths)==2 and all(p.exists() for p in paths):
            rows=[list(csv.DictReader(p.open())) for p in paths]
            checks[name]={'expected_rows':n,'row_counts':[len(x) for x in rows],'all_cells_equal':rows[0]==rows[1],'byte_identical':paths[0].read_bytes()==paths[1].read_bytes()}
        else:checks[name]={'status':'missing'}
    pairs[lane]=checks
(dest/'pair_comparison.json').write_text(json.dumps(pairs,indent=2)+'\n')
