import json, os, signal, sys, tempfile, time, unittest
from pathlib import Path
import importlib.util

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("mini", HERE / "run_local_mini_calibration.py")
mini = importlib.util.module_from_spec(spec); spec.loader.exec_module(mini)

def ref():
    return {"score": {"loss": 1.0, "normalization": {"target": 2.1, "completed_fertility": 2.1},
        "parameters": [{"parameter": p, "estimate": 1.0, "lower": 0.0, "upper": 2.0} for p in mini.PARAMETERS],
        "target_fit": [{"label": "x", "target": 1.0, "model": 1.0}]}}

class MiniTests(unittest.TestCase):
  def test_proposals_are_bounded_and_unique(self):
    r = ref(); bounds = mini.restrictions({"parameters": [{"parameter": p, "lower": 0, "upper": 2} for p in mini.PARAMETERS]}, r)
    original = {"initial_psi": .1, "parameters": {p: (0.98 if p == "beta_annual" else 1.0) for p in mini.PARAMETERS}}
    got = mini.proposals(original["parameters"], bounds)
    self.assertEqual(len(got), 12)
    self.assertEqual(len({tuple(x["parameters"][p] for p in mini.PARAMETERS) for x in got}), 12)
    self.assertTrue(all(bounds[p][0] <= x["parameters"][p] <= bounds[p][1] for x in got for p in mini.PARAMETERS))

  def test_score_comparison_ignores_paths_and_timing(self):
    a = ref(); b = ref(); b["score"]["path"] = "/tmp/a"; b["score"]["runtime_seconds"] = 900
    b["score"]["loss"] += 1e-9
    self.assertTrue(mini.compare_score(a, b))
    b["score"]["loss"] = 1.1
    self.assertFalse(mini.compare_score(a, b))

  def test_timeout_kills_process_group(self):
    with tempfile.TemporaryDirectory() as d:
        out = mini.run_process([sys.executable, "-c", "import time; time.sleep(30)"], Path(d), timeout=.05)
        self.assertEqual(out["status"], "timeout")

  def test_timeout_kills_detached_descendant(self):
    with tempfile.TemporaryDirectory() as d:
      pidfile=Path(d)/"child.pid"
      script="import subprocess,sys,time; from pathlib import Path; p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)'],start_new_session=True); Path(sys.argv[1]).write_text(str(p.pid)); time.sleep(30)"
      result=mini.run_process([sys.executable,"-c",script,str(pidfile)],Path(d)/"case",timeout=.5)
      self.assertEqual(result["status"],"timeout")
      pid=int(pidfile.read_text()); time.sleep(.05)
      import subprocess
      state=subprocess.run(["ps","-o","stat=","-p",str(pid)],capture_output=True,text=True).stdout.strip()
      self.assertTrue(not state or state.startswith("Z"),state)

  def test_atomic_write_and_fingerprint(self):
    with tempfile.TemporaryDirectory() as d:
      p = Path(d) / "x.json"; mini.atomic(p, {"x": 1}); self.assertEqual(json.loads(p.read_text()), {"x": 1})
      self.assertEqual(len(mini.sha(p)), 64)

if __name__ == "__main__":
    raise SystemExit(unittest.main(verbosity=1))
