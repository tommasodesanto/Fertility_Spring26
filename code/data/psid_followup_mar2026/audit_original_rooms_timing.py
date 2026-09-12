#!/usr/bin/env python3
"""Prepare a frozen, private sample from the author-recognized Stata source.

The original preparation block is copied, not rewritten. Only selective loading
and carrying a second rooms column are added. No source data or target changes.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
import shutil

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT.parent / 'Codes/code_per tommi_addingcontrolsandfixingthings.do'
SOURCE_SHA256 = '0262dadfb07b5a998bcbe3c32e58593556f4fb62d399ace97cb4e6189474c011'
STATA = Path('/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(work):
    assert sha(SOURCE) == SOURCE_SHA256, 'Author-recognized reference changed; review before rerunning'
    work.mkdir(parents=True, exist_ok=False)
    work.chmod(0o700)
    text = SOURCE.read_text()
    start = text.index('use "')
    end = text.index('\ncompress', start) + len('\ncompress')
    original = text[start:end]
    loading = original.splitlines()[0]
    columns = 'ID year AGEREP EDUYEAR SEX RELCHI1BYEAR MOVEDFREF_ DEATHYEAR WHYMOVED1_ ACTUALROOMS_ HOMEOWN EARNINDR IW'
    patched = original.replace('use "', 'use ' + columns + ' using "', 1)
    new_loading = patched.splitlines()[0]
    patched = patched.replace(new_loading, new_loading + '''
sort ID year
by ID: gen double rooms_aligned = ACTUALROOMS_[_n-1] if _n > 1
by ID: gen int rooms_source_year = year[_n-1] if _n > 1
replace rooms_aligned = . if !inlist(year-rooms_source_year,1,2)
''', 1)
    assert patched.count('keep ID year') == 1
    patched = patched.replace('keep ID year', 'keep rooms_aligned rooms_source_year ID year', 1)
    preamble = f'''clear all
set more off
set processors 2
version 17.0
global weight IW
global outcomes log_e
log using "{work}/prepare.log", replace text
timer clear 1
timer on 1
'''
    ending = f'''
rename first_child_year f_c_y
gen double K = year-f_c_y
quietly summarize f_c_y
gen byte lastcohort = f_c_y == r(max)
gen byte common_rooms = !missing(rooms,rooms_aligned)
gen byte original_complete = !missing(rooms,AGEREP,EDUYEAR)
gen byte aligned_complete = !missing(rooms_aligned,AGEREP,EDUYEAR)
gen byte row_has_reference = inrange(K,-2,-1) & aligned_complete
bysort f_c_y: egen byte cohort_has_reference = max(row_has_reference)
gen byte row_has_m2 = K==-2 & aligned_complete
bysort f_c_y: egen byte cohort_has_m2 = max(row_has_m2)
bysort ID: egen double birth_min = min(f_c_y)
bysort ID: egen double birth_max = max(f_c_y)
count if birth_min != birth_max
local varying_birth_rows = r(N)
drop birth_min birth_max row_has_reference row_has_m2
preserve
    gen long rows = 1
    collapse (sum) rows original_complete aligned_complete common_rooms, by(year)
    export delimited using "{work}/year_sample_counts.csv", replace
restore
preserve
    gen long rows = 1
    collapse (sum) rows original_complete aligned_complete common_rooms (max) cohort_has_reference cohort_has_m2, by(f_c_y K lastcohort)
    export delimited using "{work}/cohort_support.csv", replace
restore
keep ID year AGEREP EDUYEAR f_c_y K lastcohort rooms rooms_aligned common_rooms original_complete aligned_complete cohort_has_reference cohort_has_m2
sort ID year
isid ID year
compress
save "{work}/analysis_sample.dta", replace
timer off 1
timer list 1
di "ORIGINAL_PREPARATION_PASS varying_birth_rows=`varying_birth_rows'"
log close _all
'''
    script = work / 'prepare.do'
    script.write_text(preamble + patched + ending)
    (work / 'source_preparation_exact.txt').write_text(original)
    receipt = {'source': str(SOURCE), 'source_sha256': sha(SOURCE),
               'source_block_sha256': hashlib.sha256(original.encode()).hexdigest(),
               'generated_do_sha256': sha(script), 'status': 'running',
               'changes': ['selectively load only variables used in original preparation',
                           'before sample restrictions add previous-interview rooms with one/two-year gap',
                           'carry additional rooms column through original keep command'],
               'room_codes': 'unchanged in both columns; timing-only diagnostic',
               'timeout_seconds': 180}
    started = time.monotonic()
    process = subprocess.Popen([str(STATA), '-bq', 'do', str(script)], cwd=work)
    try:
        while process.poll() is None:
            elapsed = time.monotonic() - started
            (work / 'heartbeat.json').write_text(json.dumps({'elapsed_seconds': elapsed, 'pid': process.pid}))
            if elapsed > 180:
                process.terminate()
                process.wait(timeout=10)
                raise TimeoutError('Preparation exceeded 180-second cap')
            time.sleep(2)
        log = (work / 'prepare.log').read_text()
        assert process.returncode == 0 and 'ORIGINAL_PREPARATION_PASS varying_birth_rows=' in log
        assert sha(SOURCE) == receipt['source_sha256']
        receipt.update(status='pass', elapsed_seconds=time.monotonic()-started,
                       analysis_sample_sha256=sha(work/'analysis_sample.dta'))
        (work/'analysis_sample.dta').chmod(0o600)
    except Exception as exc:
        receipt.update(status='failed', error=str(exc))
        raise
    finally:
        (work/'preparation_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


def run_local(work, arm, timeout, processors=8, run_label=''):
    """Run one bounded local fit; retain aggregate evidence, never upload data."""
    prep = json.loads((work/'preparation_receipt.json').read_text())
    assert prep['status'] == 'pass'
    assert sha(work/'analysis_sample.dta') == prep['analysis_sample_sha256']
    assert sha(SOURCE) == SOURCE_SHA256
    script = ROOT/'code/data/psid_followup_mar2026/audit_original_rooms_timing.do'
    assert processors in (1,2,4,8)
    assert not run_label or run_label.replace('_','').isalnum()
    out = work/'local'/run_label/arm
    out.mkdir(parents=True, exist_ok=False)
    out.chmod(0o700)
    entry = out/'entry.do'
    entry.write_text(f'''clear all
set more off
set processors {processors}
version 17.0
log using "{out}/estimation.log", replace text
use "{work}/analysis_sample.dta", clear
do "{script}" {arm} "{out}"
log close _all
''')
    receipt = {'arm': arm, 'route': 'local Stata', 'status': 'running',
               'timeout_seconds': timeout, 'data_sha256': prep['analysis_sample_sha256'],
               'estimator_do_sha256': sha(script), 'source_sha256': SOURCE_SHA256,
               'microdata_uploaded': False, 'processors': processors, 'run_label': run_label}
    started = time.monotonic()
    process = subprocess.Popen([str(STATA), '-bq', 'do', str(entry)], cwd=out)
    try:
        while process.poll() is None:
            elapsed = time.monotonic()-started
            log = out/'estimation.log'
            heartbeat = {'arm': arm, 'elapsed_seconds': elapsed, 'pid': process.pid,
                         'log_bytes': log.stat().st_size if log.exists() else 0}
            (out/'heartbeat.json').write_text(json.dumps(heartbeat)+'\n')
            if elapsed > timeout:
                process.terminate()
                process.wait(timeout=15)
                raise TimeoutError(f'{arm} exceeded {timeout}-second cap')
            time.sleep(2)
        log_text = (out/'estimation.log').read_text()
        assert process.returncode == 0 and f'\nORIGINAL_TIMING_ARM_PASS {arm}\n' in log_text
        assert sha(script) == receipt['estimator_do_sha256']
        receipt.update(status='pass', elapsed_seconds=time.monotonic()-started,
                       sample_keys_sha256=sha(out/'private_sample_keys.csv'))
        (out/'private_sample_keys.csv').unlink()
    except Exception as exc:
        receipt.update(status='failed', error=str(exc), elapsed_seconds=time.monotonic()-started)
        raise
    finally:
        (out/'run_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
        evidence = ROOT/'code/data/psid_followup_mar2026/output/first_birth_correction_review/timing_local'/run_label/arm
        evidence.mkdir(parents=True, exist_ok=True)
        for name in ['coefficients.csv','covariance.csv','fitted_support.csv','fit_receipt.csv','run_receipt.json','estimation.log']:
            if (out/name).exists():
                shutil.copy2(out/name,evidence/name)
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('--run', choices=['original_native','original_common','aligned_common'])
    parser.add_argument('--timeout', type=int, default=600)
    parser.add_argument('--processors', type=int, default=8)
    parser.add_argument('--run-label', default='')
    args = parser.parse_args()
    if args.run:
        run_local(args.work.resolve(), args.run, args.timeout, args.processors, args.run_label)
    else:
        prepare(args.work.resolve())
