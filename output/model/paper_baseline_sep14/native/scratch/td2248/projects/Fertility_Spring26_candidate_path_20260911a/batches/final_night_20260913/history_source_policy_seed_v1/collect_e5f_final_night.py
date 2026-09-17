#!/usr/bin/env python3
"""Collect bounded E5F final-night stage receipts and progress snapshots."""
import argparse, hashlib, json
from pathlib import Path
def main():
 p=argparse.ArgumentParser();p.add_argument("runroot",type=Path);p.add_argument("--output",type=Path,required=True);a=p.parse_args(); out=[]
 for d in sorted(x for x in a.runroot.iterdir() if x.is_dir() and not x.name.startswith("logs")):
  row={"stage":d.name}
  for n in ("latest_completed.json","best_so_far.json","final_summary.json","summary.json","failure.json","heartbeat.json"):
   f=d/n
   if f.is_file():
    try: row[n]=json.loads(f.read_text())
    except (OSError,json.JSONDecodeError): row[n]={"malformed":True}
  row["source_hashes"]={str(f.relative_to(d)):hashlib.sha256(f.read_bytes()).hexdigest() for f in d.rglob("*") if f.is_file() and "source" in f.name}
  out.append(row)
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps({"stages":out},indent=2)+"\n")
if __name__=="__main__":main()
