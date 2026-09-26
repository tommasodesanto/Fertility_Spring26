import importlib.util
import math
from pathlib import Path
import unittest


MODULE_PATH = Path(__file__).with_name("build_inheritance_receipt_profile.py")
SPEC = importlib.util.spec_from_file_location("inheritance_profile", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
PROFILE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PROFILE)


def figure3_html(rows: str) -> str:
    return f"""
    <h5 id="fig3">Figure 3</h5><table><tr><th>Age</th></tr>{rows}</table>
    <h5 id="fig4">Figure 4</h5><table><tr><td>25</td><td>.9</td><td>.9</td><td>.9</td><td>9</td><td>9</td><td>9</td></tr></table>
    """


class Figure3ParserTests(unittest.TestCase):
    def test_parses_only_figure3_table_without_full_support_requirement(self) -> None:
        records = PROFILE.parse_figure3_html(figure3_html(
            "<tr><td>25</td><td>.1</td><td>.2</td><td>.3</td><td></td><td>1</td><td>2</td><td>3</td></tr>"
        ))
        self.assertEqual(len(records), 3)
        self.assertEqual(records[0], {"age": 25, "income_group": "bottom50", "probability_3y": .1, "conditional_amount_3y": 1.0})
        self.assertEqual(records[-1]["conditional_amount_3y"], 3.0)

    def test_duplicate_rows_rejected(self) -> None:
        row = "<tr><td>25</td><td>.1</td><td>.2</td><td>.3</td><td></td><td>1</td><td>2</td><td>3</td></tr>"
        with self.assertRaisesRegex(ValueError, "duplicate"):
            PROFILE.parse_figure3_html(figure3_html(row + row))

    def test_invalid_probability_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "outside"):
            PROFILE.parse_figure3_html(figure3_html(
                "<tr><td>25</td><td>1.1</td><td>.2</td><td>.3</td><td></td><td>1</td><td>2</td><td>3</td></tr>"
            ))

    def test_complete_support_is_separate_validation(self) -> None:
        records = PROFILE.parse_figure3_html(figure3_html(
            "<tr><td>25</td><td>.1</td><td>.2</td><td>.3</td><td></td><td>1</td><td>2</td><td>3</td></tr>"
        ))
        with self.assertRaisesRegex(ValueError, "incomplete"):
            PROFILE.validate_complete_support(records)

    def test_figure4_cannot_substitute_for_figure3(self) -> None:
        figure4_only = "<h5 id=\"fig4\">Figure 4</h5><table><tr><td>25</td><td>.9</td><td>.9</td><td>.9</td><td>9</td><td>9</td><td>9</td></tr></table>"
        with self.assertRaisesRegex(ValueError, "fig3"):
            PROFILE.parse_figure3_html(figure4_only)


class MappingTests(unittest.TestCase):
    def test_expected_amount_conservation_uses_independent_annual_comparison(self) -> None:
        record = {"age": 25, "income_group": "bottom50", "probability_3y": .3, "conditional_amount_3y": 12.0}
        mapped = PROFILE.map_three_to_four_years([record])[0]
        annual_published_expected_amount = (.3 * 12.0) / 3.0
        annual_mapped_expected_amount = (mapped["probability_4y"] * mapped["conditional_amount_4y"]) / 4.0
        self.assertTrue(math.isclose(annual_mapped_expected_amount, annual_published_expected_amount, rel_tol=1e-14))
        self.assertTrue(math.isclose(mapped["mean_amount_4y"], mapped["probability_4y"] * mapped["conditional_amount_4y"], rel_tol=1e-14))

    def test_probability_limits(self) -> None:
        zero = PROFILE.map_three_to_four_years([{"age": 25, "income_group": "bottom50", "probability_3y": 0.0, "conditional_amount_3y": 17.0}])[0]
        one = PROFILE.map_three_to_four_years([{"age": 25, "income_group": "bottom50", "probability_3y": 1.0, "conditional_amount_3y": 17.0}])[0]
        self.assertEqual((zero["probability_4y"], zero["mean_amount_4y"], zero["conditional_amount_4y"]), (0.0, 0.0, 0.0))
        self.assertEqual(one["probability_4y"], 1.0)
        self.assertTrue(math.isclose(one["conditional_amount_4y"], 68.0 / 3.0))

    def test_mapping_rejects_probability_outside_range(self) -> None:
        with self.assertRaisesRegex(ValueError, "outside"):
            PROFILE.map_three_to_four_years([{"age": 25, "income_group": "bottom50", "probability_3y": -0.01, "conditional_amount_3y": 1.0}])


if __name__ == "__main__":
    unittest.main()
