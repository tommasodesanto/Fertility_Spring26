#!/usr/bin/env python3
"""Prepare the extended frozen sample for the one-change-at-a-time sequence.

Reuses the author-recognized preparation block exactly as audit_original_rooms_timing.py
does, adding only extra carried columns (sex, relationship, household identifiers,
weight, full biological birth history, reported child count) constructed before
the block's restrictions. The base rows and both room columns are asserted
identical to the reference prepared sample.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from audit_original_rooms_timing import SOURCE, SOURCE_SHA256, STATA, sha

REFERENCE_SAMPLE = Path('/tmp/psid_rooms_window_20260923/analysis_sample.dta')
REFERENCE_SHA256 = '3466e1b734d43f3991d96b3360e261fad6c0df770e00aec864f1dcace93b340b'


def prepare(work):
    assert sha(SOURCE) == SOURCE_SHA256
    assert sha(REFERENCE_SAMPLE) == REFERENCE_SHA256
    work.mkdir(parents=True, exist_ok=False)
    work.chmod(0o700)
    text = SOURCE.read_text()
    start = text.index('use "')
    end = text.index('\ncompress', start) + len('\ncompress')
    original = text[start:end]
    loading = original.splitlines()[0]
    relchi = ' '.join(f'RELCHI{c}TYPE RELCHI{c}BYEAR' for c in range(1, 21))
    columns = ('ID year AGEREP EDUYEAR SEX REL CURRENT HHID FID RELCHIREP RELCHI1BYEAR MOVEDFREF_ '
               'DEATHYEAR WHYMOVED1_ ACTUALROOMS_ HOMEOWN EARNINDR IW ' + relchi)
    patched = original.replace('use "', 'use ' + columns + ' using "', 1)
    new_loading = patched.splitlines()[0]
    construction = '''
sort ID year
by ID: gen double rooms_aligned = ACTUALROOMS_[_n-1] if _n > 1
by ID: gen int rooms_source_year = year[_n-1] if _n > 1
replace rooms_aligned = . if !inlist(year-rooms_source_year,1,2)
gen double bio_candidate = .
forvalues c = 1/20 {
    replace bio_candidate = RELCHI`c'BYEAR if RELCHI`c'TYPE == 1 & !missing(RELCHI`c'BYEAR) & (missing(bio_candidate) | RELCHI`c'BYEAR < bio_candidate)
}
bysort ID: egen double bio_first_year = min(bio_candidate)
drop bio_candidate
bysort ID: egen double relchirep_max = max(RELCHIREP)
egen byte fid_tag = tag(HHID year FID) if CURRENT == 1 & !missing(HHID) & HHID > 0 & !missing(FID) & FID > 0
bysort HHID year: egen int n_current_fids = total(fid_tag)
drop fid_tag
gen byte woman = SEX == 2
gen byte current = CURRENT == 1
gen double rel = REL
gen double hhid = HHID
gen double fid = FID
gen double iw = IW
sort ID year
'''
    patched = patched.replace(new_loading, new_loading + construction, 1)
    assert patched.count('keep ID year') == 1
    extra = 'rooms_aligned rooms_source_year bio_first_year relchirep_max n_current_fids woman current rel hhid fid iw'
    patched = patched.replace('keep ID year', f'keep {extra} ID year', 1)
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
keep ID year AGEREP EDUYEAR f_c_y K lastcohort rooms rooms_aligned common_rooms {extra}
sort ID year
isid ID year
preserve
    keep ID year rooms rooms_aligned common_rooms
    rename (rooms rooms_aligned common_rooms) (rooms_new rooms_aligned_new common_rooms_new)
    merge 1:1 ID year using "{REFERENCE_SAMPLE}", assert(match) keepusing(rooms rooms_aligned common_rooms) nogenerate
    assert common_rooms == common_rooms_new
    assert (rooms == rooms_new) | (missing(rooms) & missing(rooms_new))
    assert (rooms_aligned == rooms_aligned_new) | (missing(rooms_aligned) & missing(rooms_aligned_new))
    di "SEQUENCE_BASE_IDENTITY_PASS"
restore
quietly count
local rows = r(N)
quietly count if missing(f_c_y) & relchirep_max == 0
local confirmed_childless_rows = r(N)
quietly count if missing(f_c_y) & !(relchirep_max == 0)
local unknown_or_untimed_rows = r(N)
quietly count if !missing(f_c_y) & (bio_first_year != f_c_y)
local bio_differs_rows = r(N)
compress
save "{work}/analysis_sample.dta", replace
timer off 1
timer list 1
di "SEQUENCE_PREPARATION_PASS rows=`rows' confirmed_childless_rows=`confirmed_childless_rows' unknown_or_untimed_rows=`unknown_or_untimed_rows' bio_differs_rows=`bio_differs_rows'"
log close _all
'''
    script = work / 'prepare.do'
    script.write_text(preamble + patched + ending)
    receipt = {'source': str(SOURCE), 'source_sha256': sha(SOURCE), 'reference_sample_sha256': REFERENCE_SHA256,
               'generated_do_sha256': sha(script), 'status': 'running',
               'changes': ['load extra shelf columns', 'construct rooms_aligned, biological first birth, reported-child maximum, family-unit count, sex/relationship/current/household/weight copies before restrictions',
                           'carry extra columns through original keep command', 'assert base rows and both room columns identical to the reference sample']}
    started = time.monotonic()
    process = subprocess.Popen([str(STATA), '-bq', 'do', str(script)], cwd=work)
    try:
        while process.poll() is None:
            if time.monotonic() - started > 600:
                process.terminate(); process.wait(timeout=10)
                raise TimeoutError('Preparation exceeded 600-second cap')
            time.sleep(2)
        log = (work / 'prepare.log').read_text()
        assert process.returncode == 0 and 'SEQUENCE_BASE_IDENTITY_PASS' in log and 'SEQUENCE_PREPARATION_PASS rows=' in log
        summary = [l for l in log.splitlines() if l.startswith('SEQUENCE_PREPARATION_PASS')][-1]
        receipt.update(status='pass', elapsed_seconds=time.monotonic()-started, summary=summary,
                       analysis_sample_sha256=sha(work/'analysis_sample.dta'))
        (work/'analysis_sample.dta').chmod(0o600)
    except Exception as exc:
        receipt.update(status='failed', error=str(exc)); raise
    finally:
        (work/'preparation_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    prepare(parser.parse_args().work.resolve())
