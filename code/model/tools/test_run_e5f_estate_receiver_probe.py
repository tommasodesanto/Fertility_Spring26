"""Torch-only smoke for the exact three-case driver loop; intentionally not run locally."""
from __future__ import annotations
import importlib.util
import tempfile
import time
import pytest
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).with_name("run_e5f_estate_receiver_probe.py")
spec = importlib.util.spec_from_file_location("estate_probe", HERE)
probe = importlib.util.module_from_spec(spec); spec.loader.exec_module(probe)

class FakeAdapter:
    CASES = probe.CASES
    def configure(self, P, case, transfer=0.): return SimpleNamespace(estate_receiver="none", case=case, transfer=0.)
    def solve_case(self, **kw):
        kw["native_solver"]()
        P = kw["parameters"]
        P.transfer = 2. if P.case == "net_valuation_transfer" else 0.
        return object(), P, [1.], {"passed":True}, {"native_solves":1}
    def estate_accounts(self, sol, P, prices):
        paid = P.transfer*.5
        return {"generated_net_period":1.,"generated_gross_period":1.,"recipient_mass":.5,"paid_period":paid,"residual":paid-1.,"transfer":P.transfer,"recipient_ages":[46,50,54,58,62]}

def test_exact_three_case_loop_with_fake_native_solver():
    with tempfile.TemporaryDirectory() as tmp:
        seen=[]
        def report(path, case, *_): seen.append(case); return {"loss":float(len(seen))}
        records=probe.run_cases(adapter=FakeAdapter(),model=object(),native_solver=lambda:None,parameters=object(),b_grid=[],initial_prices=[1.],output=Path(tmp),deadline_epoch=time.time()+10,report_case=report)
        assert [r["case"] for r in records] == list(probe.CASES) == seen
        assert records[-1]["estate_accounts"]["transfer"] == 2.
        assert records[0]["estate_accounts"]["residual"] == -1.
        assert (Path(tmp)/"latest_completed.json").is_file() and (Path(tmp)/"best_so_far.json").is_file()

def test_unfunded_receiver_fails_but_control_outflow_is_allowed(tmp_path):
    class BrokenReceiver(FakeAdapter):
        def estate_accounts(self, sol, P, prices):
            accounts = super().estate_accounts(sol,P,prices)
            accounts["residual"] = -.1
            return accounts
    with pytest.raises(RuntimeError,match="paid/generated"):
        probe.run_cases(adapter=BrokenReceiver(),model=object(),native_solver=lambda:None,
            parameters=object(),b_grid=[],initial_prices=[1.],output=tmp_path,
            deadline_epoch=time.time()+10,report_case=lambda *a:{"loss":1.})
    assert (tmp_path/"control").exists()
    assert (tmp_path/"net_valuation").exists()

def test_case_stops_at_shared_native_solve_cap(tmp_path, monkeypatch):
    monkeypatch.setattr(probe,"MAX_SOLVES",2)
    with pytest.raises(TimeoutError,match="shared native"):
        probe.run_cases(adapter=FakeAdapter(),model=object(),native_solver=lambda:None,
            parameters=object(),b_grid=[],initial_prices=[1.],output=tmp_path,
            deadline_epoch=time.time()+10,report_case=lambda *a:{"loss":1.})
