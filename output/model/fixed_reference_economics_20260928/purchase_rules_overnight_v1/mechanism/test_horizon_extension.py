"""Zero-solve checks for the isolated 48/64-date numerical extension."""
from __future__ import annotations

import difflib
import re
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEPLOY = HERE.parent / "mechanism_deployment"


class ExtensionTests(unittest.TestCase):
    def test_driver_changes_only_horizon_routing(self):
        base = (HERE / "run_case.py").read_text().splitlines()
        extended = (HERE / "run_case_extended.py").read_text().splitlines()
        changes = [line for line in difflib.unified_diff(base, extended, n=0)
                   if line.startswith(("+", "-")) and not line.startswith(("+++", "---"))]
        self.assertEqual(changes, [
            '-"""One bounded, source-authenticated dated financing experiment.',
            '-',
            '-Run separately for each purchase rule, policy duration, and 12/16-date horizon.',
            '+"""One bounded, source-authenticated longer-horizon financing diagnostic.',
            '+',
            '+Run separately for each purchase rule, policy duration, and 48/64-date horizon.',
            '-    ap.add_argument("--horizon", type=int, choices=(1,12,16), required=True)',
            '+    ap.add_argument("--horizon", type=int, choices=(1,12,16,48,64), required=True)',
            '-            or (not args.smoke_one_date and args.horizon in (12, 16)),',
            '-            "One-date smoke must be a one-period control; production horizons are 12 or 16")',
            '+            or (not args.smoke_one_date and args.horizon in (12, 16, 48, 64)),',
            '+            "One-date smoke must be a one-period control; diagnostic horizons are 12, 16, 48, or 64")',
            '-                if args.horizon == 16: J = extend_measured_jacobian(derivative_receipt, 16)',
            '+                if args.horizon > 12: J = extend_measured_jacobian(derivative_receipt, args.horizon)',
        ])
        self.assertTrue(any('write(folder / "terminal_checks.json", terminal_check)' in line
                            for line in extended))

    def test_exact_loop_and_budget_contract(self):
        launcher = (DEPLOY / "launch_torch_extended.sh").read_text()
        submit = (DEPLOY / "submit_torch_extended.sh").read_text()
        def vector(name):
            return re.search(rf"^{name}=\(([^)]*)\)$", launcher, re.M).group(1).split()
        self.assertEqual(vector("arms"), ["hard"]*6 + ["quarter"]*6)
        self.assertEqual(vector("kinds"), ["control","temporary","permanent"]*4)
        self.assertEqual(vector("horizons"), ["48"]*3+["64"]*3+["48"]*3+["64"]*3)
        for item in ("cap=1024", "deadline_epoch=$((start_epoch+14400))",
                     "deadline_epoch=1790949600", "--mem=32G", "--cpus-per-task=1",
                     "run_case_extended.py"):
            self.assertIn(item, launcher)
        self.assertIn('--dependency="afterok:$smoke"', submit)
        self.assertIn('--array=0-11%12', submit)


if __name__ == "__main__":
    unittest.main()
