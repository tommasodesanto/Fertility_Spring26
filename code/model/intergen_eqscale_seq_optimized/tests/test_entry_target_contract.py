"""Regression tests for the empirical outside-entry target contract."""

from pathlib import Path
from types import SimpleNamespace
import sys
import unittest


TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from entry_target_contract import entry_target_from_outside_origin_share


class EntryTargetContractTests(unittest.TestCase):
    def test_outside_origin_share_maps_to_candidate_specific_entry_probability(self) -> None:
        solution = SimpleNamespace(
            entry_rate=0.06173345618074889,
            entrants_mature_total=0.05292942076730084,
        )
        parameters = SimpleNamespace(local_birth_entry_weight=1.0)

        contract = entry_target_from_outside_origin_share(
            solution,
            parameters,
            outside_origin_share=0.169,
        )

        self.assertAlmostEqual(contract["target_qbar"], 0.9692247023019598)
        self.assertAlmostEqual(contract["outside_origin_entrant_share_identity_check"], 0.169)

    def test_infeasible_empirical_share_is_rejected(self) -> None:
        solution = SimpleNamespace(entry_rate=1.0, entrants_mature_total=0.2)
        parameters = SimpleNamespace(local_birth_entry_weight=1.0)

        with self.assertRaisesRegex(RuntimeError, "infeasible model entry probability"):
            entry_target_from_outside_origin_share(
                solution,
                parameters,
                outside_origin_share=0.169,
            )


if __name__ == "__main__":
    unittest.main()
