#!/usr/bin/env python3
"""Replay the pinned first-birth rooms estimator with extra support receipts.

The production .do file, target, sample and regression are unchanged. Generated
Stata code runs in a separate output folder. No microdata are exported.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BUILDER = HERE / 'sa_rooms_first_birth_household_aligned_v1.do'
PRIMARY = HERE / 'output/sa_rooms_first_birth_household_aligned_v1'
DEFAULT_OUTPUT = ROOT / 'output/model/e5f_first_birth_measurement_review_20260905a'
STATA = Path('/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp')


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def generate(output: Path, mode: str) -> str:
    text = BUILDER.read_text()
    old = 'local outdir  "`outroot\'/sa_rooms_first_birth_household_aligned_v1"'
    assert text.count(old) == 1
    text = text.replace(old, f'local outdir  "{output}"')
    anchor = 'di as text "Running weighted Sun-Abraham event study:'
    assert text.count(anchor) == 1
    pre = r'''
* Supplemental input-support receipts; no estimation changes.
assert `input_obs' == 49872
assert `input_ids' == 4527
assert `treated_ids' == 2486
* These establish that the input-support table is the cohort-share sample.
assert !missing(ID, year, AGEREP, EDUYEAR, IW, rooms)
bysort ID: assert first_child_year == first_child_year[1]
egen byte audit_event_total = rowtotal(L*event F*event)
assert audit_event_total == (K != -2) if !never_treated
assert audit_event_total == 0 if never_treated
drop audit_event_total
preserve
    gen byte post1997 = year >= 1997
    gen long audit_one = 1
    collapse (sum) observations=audit_one weight_sum=IW, ///
        by(first_child_year K post1997 never_treated)
    export delimited using "`outdir'/input_cohort_event_support.csv", replace
restore
'''
    if mode == 'smoke':
        pre += '\ndi as result "FIRST_BIRTH_SUPPORT_SMOKE_PASS"\nlog close _all\nexit 0\n'
    text = text.replace(anchor, pre + '\n' + anchor)
    anchor2 = 'control_cohort(never_treated) covariates(i.AGEREP i.EDUYEAR)\n'
    assert text.count(anchor2) == 1
    capture_objects = r'''
* Capture original estimation objects before any preserve/clear/export.
gen byte audit_estimation_sample = e(sample)
bysort ID: egen byte audit_observed_kminus2 = ///
    max(audit_estimation_sample & K == -2)
foreach matrix_name in b_iw V_iw ff_w b_interact V_interact {
    matrix audit_`matrix_name' = e(`matrix_name')
}
'''
    text = text.replace(anchor2, anchor2 + capture_objects)
    extra = r'''
* Export captured objects after all primary estimator calculations.
foreach matrix_name in b_iw V_iw ff_w b_interact V_interact {
    matrix audit_matrix = audit_`matrix_name'
    preserve
        clear
        svmat2 audit_matrix, names(col) rnames(row_name)
        export delimited using "`outdir'/`matrix_name'.csv", replace
    restore
}
preserve
    keep if audit_estimation_sample
    gen byte post1997 = year >= 1997
    gen long audit_one = 1
    gen double audit_weighted_rooms = IW * rooms
    gen byte audit_nonpositive_rooms = rooms <= 0
    collapse (sum) observations=audit_one weight_sum=IW ///
        weighted_rooms_sum=audit_weighted_rooms ///
        nonpositive_room_rows=audit_nonpositive_rooms, ///
        by(first_child_year K post1997 never_treated audit_observed_kminus2)
    export delimited using "`outdir'/estimation_cohort_event_support.csv", replace
restore
di as result "FIRST_BIRTH_SUPPORT_OBJECTS_SAVED"
'''
    assert text.count('timer off 1') == 1
    return text.replace('timer off 1', extra + '\ntimer off 1')


def run(mode: str, outroot: Path, timeout: int) -> None:
    output = outroot / mode
    output.mkdir(parents=True, exist_ok=True)
    if (output / 'run_receipt.json').exists():
        raise RuntimeError('Refusing to overwrite an existing audit run')
    primary_hashes = {p.name: sha(p) for p in PRIMARY.iterdir() if p.is_file()}
    source_hash = sha(BUILDER)
    do_path = output / 'replay.do'
    do_path.write_text(generate(output, mode))
    receipt = {'mode': mode, 'status': 'running', 'started_at': time.time(),
               'source_builder': str(BUILDER), 'source_builder_sha256': source_hash,
               'generated_do_sha256': sha(do_path), 'timeout_seconds': timeout,
               'primary_output_hashes_before': primary_hashes,
               'regression_specification_changed': False,
               'microdata_exported': False}
    (output / 'run_receipt.json').write_text(json.dumps(receipt, indent=2))
    start = time.monotonic()
    try:
        p = subprocess.Popen([str(STATA), '-bq', 'do', str(do_path)], cwd=output)
        receipt['pid'] = p.pid
        (output / 'run_receipt.json').write_text(json.dumps(receipt, indent=2))
        while p.poll() is None:
            elapsed = time.monotonic() - start
            log = output / 'sa_rooms_first_birth_household_aligned_v1.log'
            heartbeat = {'status': 'running', 'elapsed_seconds': elapsed,
                         'log_bytes': log.stat().st_size if log.exists() else 0,
                         'updated_at': time.time()}
            (output / 'heartbeat.json').write_text(json.dumps(heartbeat, indent=2))
            if elapsed > timeout:
                p.terminate()
                p.wait(timeout=15)
                raise TimeoutError('Declared Stata audit budget exceeded')
            time.sleep(2)
        log_text = (output / 'sa_rooms_first_birth_household_aligned_v1.log').read_text()
        marker = 'FIRST_BIRTH_SUPPORT_SMOKE_PASS' if mode == 'smoke' else 'CORRECTED_FIRST_BIRTH_ROOMS_TARGET estimate='
        if p.returncode != 0 or marker not in log_text:
            raise RuntimeError('Stata did not complete the declared audit; inspect replay.log')
        if sha(BUILDER) != source_hash:
            raise RuntimeError('Primary builder changed during the audit')
        after = {p.name: sha(p) for p in PRIMARY.iterdir() if p.is_file()}
        if after != primary_hashes:
            raise RuntimeError('Primary target outputs changed during the audit')
        receipt['primary_outputs_unchanged'] = True
        receipt['sample_design_assertions'] = {
            'all_regression_inputs_nonmissing': True,
            'birth_cohort_constant_within_person': True,
            'event_indicators_partition_all_treated_rows_except_kminus2': True,
            'never_treated_event_indicators_zero': True,
        }
        if mode == 'full':
            old, new = rows(PRIMARY / 'target_receipt.csv')[0], rows(output / 'target_receipt.csv')[0]
            deltas = {}
            for k, value in old.items():
                if k == 'runtime_seconds':
                    continue
                try:
                    a, b = float(value), float(new[k])
                except ValueError:
                    if value != new[k]:
                        raise RuntimeError(f'Receipt text differs: {k}')
                    continue
                deltas[k] = b - a
                if not math.isclose(a, b, rel_tol=1e-11, abs_tol=2e-10):
                    raise RuntimeError(f'Exact replay differs on {k}: {a} versus {b}')
            receipt['numeric_receipt_deltas'] = deltas
            receipt['target_estimate'] = float(new['estimate'])
            receipt['target_standard_error'] = float(new['standard_error'])
        receipt['status'] = 'pass'
    except BaseException as error:
        receipt['status'] = 'failed'
        receipt['error'] = str(error)
        raise
    finally:
        receipt['elapsed_seconds'] = time.monotonic() - start
        receipt['source_builder_unchanged'] = sha(BUILDER) == source_hash
        receipt['primary_outputs_unchanged'] = (
            {p.name: sha(p) for p in PRIMARY.iterdir() if p.is_file()} == primary_hashes
        )
        (output / 'run_receipt.json').write_text(json.dumps(receipt, indent=2))
        (output / 'heartbeat.json').write_text(json.dumps(receipt, indent=2))
    print(json.dumps({k: receipt[k] for k in ('mode', 'status', 'elapsed_seconds')}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['smoke', 'full'])
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--timeout', type=int, default=900)
    args = parser.parse_args()
    run(args.mode, args.output.resolve(), args.timeout)
