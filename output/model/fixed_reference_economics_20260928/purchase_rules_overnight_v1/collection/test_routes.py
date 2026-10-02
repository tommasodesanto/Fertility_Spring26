"""Synthetic receipt routing test; no model or cluster access."""
from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path

os.environ['PURCHASE_PACKET_ROOT'] = str(Path(__file__).resolve().parent.parent)
import scan_remote as scan
from collect import select_arm


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + '\n')


class RestartRouteTests(unittest.TestCase):
    def test_local_worker_receipt_route_uses_local_schema(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            parent, restart = base / 'local10_v1/chain50', base / 'restart_v2/runs/chain50'
            old_parent = scan.PARENT_RESULTS
            scan.PARENT_RESULTS = base / 'local10_v1'
            try:
                write(scan.PARENT_RESULTS / 'pids.json', {'50': dict(started_epoch=1, deadline_epoch=14401)})
                write(parent / 'worker_terminal.json', dict(chain=50, deadline_epoch=14401,
                                                              status='search_no_selected_candidate'))
                write(parent / 'search/completed.json', dict(objective_calls=1, selected=None,
                    search_stop_reason='native_evaluation_budget_exhausted'))
                write(restart / 'search/completed.json', dict(objective_calls=2,
                    selected=dict(loss=3.)))
                write(restart / 'restart_contract.json', dict(chain=50, parent=str(parent),
                    parent_search_receipt_sha256=scan.sha(parent / 'search/completed.json'),
                    parent_worker_receipt_sha256=scan.sha(parent / 'worker_terminal.json'),
                    parent_postcheck_receipt_sha256=None, original_start_epoch=1,
                    original_deadline_epoch=14401, original_objective_calls=1,
                    remaining_objective_calls=249, no_clock_reset=True,
                    no_call_count_reset=True, no_native_cap_change=True,
                    reviewed_optimizer_function_sha256=
                        '9f72add8cced622aae2532ca231231a27386d22f601cdcf3f2e6f2fa3ac96c1b',
                    optimizer_source_sha256='a'*64,
                    original_target_contract_sha256=scan.canonical(scan.CONTRACT)))
                write(restart / 'restart_summary.json', dict(chain=50, winner='restart',
                    status='restart_search_finished',
                    cumulative_objective_calls=3, original_deadline_epoch=14401,
                    parent_selected_loss=None, restart_selected_loss=3.,
                    selected_requires_fresh_postcheck=True))
                self.assertEqual(scan.local_restart_provenance(restart, 50)['restart_state'], 'postcheck_running')
                write(restart / 'postcheck/completed.json', dict(status='selected_numerically_verified'))
                summary = json.loads((restart / 'restart_summary.json').read_text())
                summary['status'] = 'restart_selected_numerically_verified'
                summary['postcheck_exit_code'] = 0
                write(restart / 'restart_summary.json', summary)
                self.assertEqual(scan.local_restart_provenance(restart, 50)['restart_state'], 'verified')
                summary['status'] = 'restart_selected_postcheck_failed'
                write(restart / 'restart_summary.json', summary)
                self.assertEqual(scan.local_restart_provenance(restart, 50)['restart_state'],
                                 'restart_selected_postcheck_failed')
            finally:
                scan.PARENT_RESULTS = old_parent

    def test_only_fresh_postchecked_improvement_can_win(self):
        rows = [dict(chain=2, arm='hard', status='postchecked', loss=4., source_run='original'),
                dict(chain=2, arm='hard', status='restart_pending', loss=1., source_run='restart'),
                dict(chain=3, arm='quarter', status='postchecked', loss=5., source_run='original')]
        self.assertEqual(select_arm(rows, 'hard')['source_run'], 'original')
        rows[1].update(status='postchecked', loss=3.)
        self.assertEqual(select_arm(rows, 'hard')['source_run'], 'restart')
        self.assertEqual(select_arm(rows, 'quarter')['chain'], 3)

    def test_restart_parent_hashes_budget_and_strict_improvement(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            parent, restart = base / 'parent/chain_2', base / 'restart/chain_2'
            old_parent = scan.PARENT_RESULTS
            scan.PARENT_RESULTS = base / 'parent'
            try:
                write(parent / 'launcher_start.json', dict(chain=2, start_epoch=1,
                                                          deadline_epoch=14401))
                write(parent / 'search/search_completed.json', dict(objective_calls=10,
                    selected=dict(loss=4.), search_stop_reason='native_evaluation_budget_exhausted'))
                write(parent / 'postcheck/completed.json', dict(status='selected_numerically_verified'))
                write(restart / 'search/search_completed.json', dict(objective_calls=5,
                    selected=dict(loss=3.)))
                write(restart / 'restart_contract.json', dict(chain=2, parent=str(parent),
                    parent_search_receipt_sha256=scan.sha(parent / 'search/search_completed.json'),
                    parent_launcher_receipt_sha256=scan.sha(parent / 'launcher_start.json'),
                    parent_postcheck_receipt_sha256=scan.sha(parent / 'postcheck/completed.json'),
                    original_start_epoch=1, original_deadline_epoch=14401,
                    original_objective_calls=10, remaining_objective_calls=240,
                    no_clock_reset=True, no_call_count_reset=True,
                    optimizer_source_sha256='a'*64,
                    original_target_contract_sha256=scan.canonical(scan.CONTRACT)))
                write(restart / 'restart_summary.json', dict(chain=2, winner='restart',
                    cumulative_objective_calls=15, original_deadline_epoch=14401,
                    parent_selected_loss=4., restart_selected_loss=3.,
                    selected_requires_fresh_postcheck=True))
                provenance = scan.restart_provenance(restart, 2)
                self.assertEqual(provenance['parent_remote_root'], str(parent))
                self.assertEqual(provenance['restart_winner'], 'restart')
                contract = json.loads((restart / 'restart_contract.json').read_text())
                contract['parent'] = '/work/parent'
                write(restart / 'restart_contract.json', contract)
                self.assertEqual(scan.restart_provenance(restart, 2)['restart_winner'], 'restart')
                summary = json.loads((restart / 'restart_summary.json').read_text())
                summary['cumulative_objective_calls'] = 251
                write(restart / 'restart_summary.json', summary)
                with self.assertRaisesRegex(ValueError, 'budget drift'):
                    scan.restart_provenance(restart, 2)
            finally:
                scan.PARENT_RESULTS = old_parent


if __name__ == '__main__':
    unittest.main()
