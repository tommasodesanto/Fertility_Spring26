"""Small accounting tests for the selectable transition entry clocks."""

import unittest

import run_e5f_open_population_transition as transition
from intergen_eqscale_seq_optimized.adult_entry import SplitBirthEntryQueue


class AdultEntryCallerTests(unittest.TestCase):
    def test_split_cohort_injects_half_at_dates_four_and_five(self):
        queue = SplitBirthEntryQueue.constant_prehistory(0.0)
        entry_next = []
        for t in range(7):
            due, queue = transition.advance_adult_entry_clock(
                queue, 2.1 if t == 0 else 0.0, 1 / 2.1, "split-16-20"
            )
            entry_next.append(due)
        self.assertEqual(entry_next, [0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0])

    def test_legacy_single_queue_remains_twenty_years(self):
        queue = [0.0] * 4
        entry_next = []
        for t in range(7):
            due, queue = transition.advance_adult_entry_clock(
                queue, 2.1 if t == 0 else 0.0, 1 / 2.1, "legacy-20"
            )
            entry_next.append(due)
        self.assertEqual(entry_next, [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0])

    def test_constant_prehistory_no_transitional_entry_jump(self):
        queue = SplitBirthEntryQueue.constant_prehistory(21.0)
        for _ in range(8):
            due, queue = transition.advance_adult_entry_clock(
                queue, 21.0, 1 / 2.1, "split-16-20"
            )
            self.assertAlmostEqual(due, 10.0)


if __name__ == "__main__":
    unittest.main()
