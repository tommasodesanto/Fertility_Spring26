from pathlib import Path
import subprocess, os, time, json, signal, datetime
root=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
out=root/'output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/fable_reply'
now=lambda: datetime.datetime.now(datetime.timezone.utc).isoformat()
env=os.environ.copy()
for key in ['ANTHROPIC_API_KEY','ANTHROPIC_AUTH_TOKEN','ANTHROPIC_BASE_URL']:
    env.pop(key,None)
auth=subprocess.run(['/Users/tommasodesanto/.local/bin/claude','auth','status'],env=env,capture_output=True,text=True,check=True)
a=json.loads(auth.stdout)
if not (a.get('loggedIn') and a.get('authMethod')=='claude.ai' and a.get('subscriptionType')=='max' and a.get('apiProvider')=='firstParty'):
    raise SystemExit('Claude Max first-party authentication required; no alternate provider used.')
cmd=['/Users/tommasodesanto/.local/bin/claude','--model','fable','--effort','max','--print','--output-format','stream-json','--verbose','--max-turns','30','--tools','Read,Glob,Grep,WebSearch,WebFetch','--allowedTools','Read,Glob,Grep,WebSearch,WebFetch','--permission-mode','dontAsk','--strict-mcp-config']
r={'status':'running','started_utc':now(),'time_limit_seconds':900,'max_turns':30,'auth':'claude.ai Max firstParty','model_requested':'fable','effort':'max','prompt':'prompt.md'}
(out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n')
t=time.monotonic()
with (out/'stream.jsonl').open('w') as log,(out/'stderr.txt').open('w') as err:
    p=subprocess.Popen(cmd,cwd=root,env=env,stdin=subprocess.PIPE,stdout=log,stderr=err,text=True,start_new_session=True)
    r['pid']=p.pid
    (out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n')
    try:
        p.communicate((out/'prompt.md').read_text(),timeout=900)
        r['status']='completed' if p.returncode==0 else 'failed'
    except subprocess.TimeoutExpired:
        os.killpg(p.pid,signal.SIGTERM)
        try:p.wait(timeout=10)
        except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
        r['status']='timed_out'
r.update(exit_code=p.returncode,finished_utc=now(),elapsed_seconds=round(time.monotonic()-t,1))
final=None
for line in (out/'stream.jsonl').read_text().splitlines():
    try: item=json.loads(line)
    except json.JSONDecodeError: continue
    if item.get('type')=='system' and item.get('subtype')=='init':
        r['session_id']=item.get('session_id');r['model_actual']=item.get('model')
    if item.get('type')=='result':
        final=item.get('result');r['result_subtype']=item.get('subtype');r['is_error']=item.get('is_error')
if final:(out/'final.md').write_text(final+'\n')
(out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
