#!/usr/bin/env python3
"""Run the unchanged history driver with independently verified exact reuse."""
import argparse
import hashlib
import json
from pathlib import Path
import time

def cache_size_gib(count, requested=None):
    default = 24 if count == 100 else 6
    if requested is None:
        return default
    if type(requested) is not int or not 1 <= requested <= default:
        raise ValueError('A cache override must reduce the existing memory cap')
    return requested

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--cache-proof',type=Path,required=True)
    parser.add_argument('--cache-sha256',required=True)
    parser.add_argument('--cache-gib',type=int,help='Optional smaller exact-policy cache; no model change')
    parser.add_argument('driver_args',nargs=argparse.REMAINDER)
    args=parser.parse_args()
    import e5f_exact_policy_cache as cache
    proof=json.loads(args.cache_proof.read_text())
    digest=hashlib.sha256(Path(cache.__file__).read_bytes()).hexdigest()
    if (digest!=args.cache_sha256 or proof.get('cache_sha256')!=digest
        or proof.get('status')!='verified' or proof.get('exact_mapping_equal') is not True
        or proof.get('household_and_accounting_mapping_valid') is not True):
        raise ValueError('Exact policy cache lacks its matching native forecast certificate')
    import e5f_rebated_surprises as rebated
    import run_e5f_final_rebated_history as driver
    _,joined,*_=rebated._runtime()
    argv=args.driver_args[1:] if args.driver_args[:1]==['--'] else args.driver_args
    output=Path(argv[argv.index('--output')+1]);started=time.monotonic()
    count=int(argv[argv.index('--count')+1])
    cache_gib=cache_size_gib(count,args.cache_gib)
    with cache.policy_cache(joined.pf,max_bytes=cache_gib*1024**3) as stats:
        try:driver.main(argv)
        finally:
            driver.save(output/'policy_cache_receipt.json',dict(cache_sha256=digest,
                cache_gib=cache_gib,
                proof=str(args.cache_proof),proof_sha256=hashlib.sha256(args.cache_proof.read_bytes()).hexdigest(),
                elapsed_seconds=time.monotonic()-started,**stats.snapshot()))

if __name__=='__main__':main()
