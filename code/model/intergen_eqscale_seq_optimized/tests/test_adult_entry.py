"""Accounting tests for aggregate entry, independent of the household DP."""

import unittest

from intergen_eqscale_seq_optimized.adult_entry import (
    SplitBirthEntryQueue,
    adjusted_births,
    potential_entry_households,
    require_closed_stationary_renewal,
)


class AdultEntryTests(unittest.TestCase):
    def test_top_bin_and_single_conversion(self):
        children = adjusted_births(100.0, 10.0, 3.6)
        self.assertAlmostEqual(children, 106.0)
        self.assertAlmostEqual(potential_entry_households(children), 106.0 / 2.1)

    def test_impulse_arrives_at_16_and_20_years(self):
        queue = SplitBirthEntryQueue.constant_prehistory(0.0)
        due = []
        for t in range(7):
            flow, queue = queue.step(2.1 if t == 0 else 0.0)
            due.append(flow)
        self.assertEqual(due, [0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0])
        # Each returned flow is injected at t+1, hence dates 4 and 5.

    def test_constant_births_and_prehistory(self):
        queue = SplitBirthEntryQueue.constant_prehistory(21.0)
        self.assertAlmostEqual(queue.stock, 35.0)  # 3*5 + 4*5 pending households
        for _ in range(9):
            due, queue = queue.step(21.0)
            self.assertAlmostEqual(due, 10.0)
            self.assertAlmostEqual(queue.stock, 35.0)

    def test_queue_mass_and_parental_death_independence(self):
        queue = SplitBirthEntryQueue.constant_prehistory(0.0)
        parent_alive = True
        for t, births in enumerate((2.1, 0.0, 0.0, 0.0, 0.0)):
            old_stock = queue.stock
            if t == 1:
                parent_alive = False
            due, queue = queue.step(births)
            self.assertAlmostEqual(queue.stock, old_stock + births / 2.1 - due)
            if t == 3:
                self.assertFalse(parent_alive)
                self.assertAlmostEqual(due, 0.5)
            if t == 4:
                self.assertAlmostEqual(due, 0.5)

    def test_closed_stationary_basis_and_fertility_gate(self):
        E = 0.06173345618094337
        births = 0.12964029045028488
        receipt = require_closed_stationary_renewal(E, births, 1e-6)
        self.assertAlmostEqual(receipt["potential_B"], births / 2.1)
        self.assertLess(abs(receipt["entry_residual"]), 2e-8)
        with self.assertRaises(ValueError):
            require_closed_stationary_renewal(E, births * 0.9, 1e-6)


if __name__ == "__main__":
    unittest.main()
