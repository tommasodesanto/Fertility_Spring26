"""Zero-model-solve check of local source pins, receipts, and censored loop."""
import json
import sys
import tempfile
import time
from pathlib import Path

from controller import LOCAL, original_contract, install_censored_optimizer


def main():
    checked = []
    for chain in (48, 50, 52):
        parent, launch, terminal, search, selected, target = original_contract(chain)
        assert parent.is_dir() and terminal['deadline_epoch'] == launch['deadline_epoch']
        checked.append(dict(chain=chain, calls=search['objective_calls'],
                            valid_selected=selected is not None))
    sys.argv = ['run_local_psi.py', '--chain', '48']
    sys.path.insert(0, str(LOCAL))
    import run_local_psi as original
    source_hash = install_censored_optimizer(original, 8)
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)/'toy'
        calls = [0]
        def evaluate(label, point, deadline):
            calls[0] += 1
            if calls[0] == 2:
                return dict(status='budget_exhausted', reason='uncomputed_bounded_budget',
                            lifecycle_solves=31)
            return dict(status='passed', residual=[float(calls[0])]+[0.]*9,
                        lifecycle_solves=1, report=str(out/label))
        result = original.optimize(out, {'x': .5}, {'x': (0., 1.)}, ('x',),
                                   evaluate, time.time()+3600, toy=True)
        cases = json.loads((out/'cases.json').read_text())
        assert len(cases) >= 3 and cases[0]['status'] == 'passed'
        assert cases[1]['optimizer_only_censor'] and cases[1]['loss'] is None
        assert cases[1]['objective'] is None and cases[1]['computed_valid_loss'] is False
        assert cases[2]['status'] == 'passed' and result['selected']['status'] == 'passed'
        assert result['objective_calls'] <= 8
        late = original.optimize(Path(tmp)/'late', {'x': .5}, {'x': (0., 1.)},
                                 ('x',), evaluate, time.time()+1, toy=True)
        assert late['objective_calls'] == 0
        assert late['search_stop_reason'] == 'four_hour_actual_start_final_reserve'
    print(json.dumps(dict(status='passed_zero_model_solves', checked=checked,
                          optimizer_source_sha256=source_hash,
                          censored_case_did_not_select=True)))


if __name__ == '__main__':
    main()
