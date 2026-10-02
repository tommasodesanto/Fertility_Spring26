"""Run one local search and its selected-point postcheck under one deadline."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def write(path: Path, payload: dict) -> None:
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
    temp.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--chain', type=int, choices=range(48, 58), required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--deadline-epoch', type=float, required=True)
    args = parser.parse_args()
    assert os.environ.get('ALLOW_LOCAL_CALIBRATION') == '1'
    out = args.out.resolve()
    assert out.is_relative_to(HERE.resolve()) and not out.exists()
    out.mkdir(parents=True)
    common = [sys.executable, str(HERE/'bootstrap.py'), '--chain', str(args.chain),
              '--deadline-epoch', str(args.deadline_epoch)]
    search = out/'search'
    postcheck = out/'postcheck'
    search_command = common + ['--out', str(search), '--fast-objective']
    search_exit = subprocess.run(search_command, check=False).returncode
    terminal = dict(chain=args.chain, deadline_epoch=args.deadline_epoch,
                    search_exit_code=search_exit, search_dir=str(search),
                    postcheck_dir=str(postcheck), no_automatic_retry=True)
    result = search/'completed.json'
    if search_exit != 0 or not result.is_file():
        terminal['status'] = 'search_failed'
    else:
        selected = json.loads(result.read_text()).get('selected')
        if selected is None:
            terminal['status'] = 'search_no_selected_candidate'
        elif time.time() >= args.deadline_epoch:
            terminal['status'] = 'selected_postcheck_deadline_elapsed'
        else:
            # Fresh process, same absolute four-hour cap, and exact selected point.
            check_command = common + ['--out', str(postcheck), '--verify-only', str(result)]
            check_exit = subprocess.run(check_command, check=False).returncode
            terminal['postcheck_exit_code'] = check_exit
            check_result = postcheck/'completed.json'
            terminal['status'] = (json.loads(check_result.read_text()).get('status')
                                  if check_exit == 0 and check_result.is_file()
                                  else 'selected_postcheck_failed')
    terminal['finished_epoch'] = time.time()
    write(out/'worker_terminal.json', terminal)
    return 0 if terminal['status'] == 'selected_numerically_verified' else 1


if __name__ == '__main__':
    raise SystemExit(main())
