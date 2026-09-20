"""Tiny deterministic contract checks for the ACS second-birth proxy.

These tests use hand-built rows only.  They intentionally do not load ACS data,
run R, or implement the estimator.
"""

from dataclasses import dataclass
from typing import Optional, Sequence
import unittest


@dataclass(frozen=True)
class MotherRow:
    year: int
    age: int
    nchild: int
    linked_ages: Sequence[int]
    fertyr: Optional[int] = None


def proxy_clock(row: MotherRow, *, strict: bool = True):
    """Return the clock; age-band/FERTYR eligibility is separate below."""
    ages = sorted(row.linked_ages, reverse=True)
    if any(age < 0 for age in ages):
        return None, "negative linked child age"
    if row.nchild < 2 or len(ages) < 2 or ages[0] < 1:
        return None, "no older child plus two linked children"
    if strict and (row.nchild != 2 or len(ages) != 2):
        return None, "strict sample excludes third/mismatched roster"
    if ages[0] <= ages[1] or (len(ages) >= 3 and ages[1] <= ages[2]):
        return None, "ambiguous age tie"
    t = ages[1]
    return {
        "event_time": t,
        "event_year": row.year - t,
        "age_at_event": row.age - t,
        "birth_gap": ages[0] - ages[1],
    }, "ok"


def strict_eligibility(row: MotherRow):
    """Apply the declared strict sample rules to the clock helper output."""
    clock, reason = proxy_clock(row, strict=True)
    if clock is None:
        return None, reason
    if row.nchild != 2 or len(row.linked_ages) != 2:
        return None, "strict sample excludes third/mismatched roster"
    if not 25 <= clock["age_at_event"] <= 45:
        return None, "age at proxy birth outside 25-45"
    if clock["event_time"] == 0 and row.fertyr is not None and row.fertyr != 2:
        return None, "observed FERTYR is not yes"
    if clock["birth_gap"] < 1:
        return None, "nonpositive birth gap"
    return clock, "ok"


def anchor_from_event0(row: MotherRow):
    """Create a donor anchor only from a strict event-time-zero row."""
    clock, reason = strict_eligibility(row)
    if clock is None:
        return None, reason
    if clock["event_time"] != 0:
        return None, "donor anchor requires event time zero"
    return dict(clock, is_event0_anchor=True), "ok"


def donor_target(anchor, k: int):
    """Target one-child donor cell for an anchor and relative time k."""
    if not anchor.get("is_event0_anchor", False):
        raise ValueError("donor targets require an event-0 anchor")
    return {
        "year": anchor["event_year"] + k,
        "mother_age": anchor["age_at_event"] + k,
        "sole_child_age": anchor["birth_gap"] + k,
    }


def linked_children(rows, mother):
    """MOMLOC linkage scoped to the complete household key."""
    key = (mother["year"], mother["sample"], mother["serial"])
    return [r for r in rows if
            (r["year"], r["sample"], r["serial"]) == key and
            r["momloc"] == mother["pernum"]]


class SecondBirthProxyContractTests(unittest.TestCase):
    def test_post_clock_and_full_pre_support_for_gap_five(self):
        anchor_row = MotherRow(2022, 34, 2, [5, 0], fertyr=2)
        clock, reason = anchor_from_event0(anchor_row)
        self.assertEqual(reason, "ok")
        self.assertEqual(clock["event_time"], 0)
        self.assertEqual(clock["birth_gap"], 5)
        # g=5 is the boundary for every negative target -5,...,-1.
        self.assertEqual([donor_target(clock, k)["sole_child_age"]
                          for k in range(-5, 0)], [0, 1, 2, 3, 4])
        # A later cross-section maps independently and cannot create an anchor.
        post, post_reason = proxy_clock(MotherRow(2024, 36, 2, [7, 2]))
        self.assertEqual(post_reason, "ok")
        self.assertEqual(post["event_time"], 2)
        self.assertEqual(post["event_year"], 2022)

    def test_gap_one_cannot_supply_reference_minus_two(self):
        row = MotherRow(2022, 30, 2, [1, 0], fertyr=2)
        clock, reason = anchor_from_event0(row)
        self.assertEqual(reason, "ok")
        self.assertEqual(clock["birth_gap"], 1)
        self.assertEqual(donor_target(clock, -2)["sole_child_age"], -1)
        self.assertLess(donor_target(clock, -2)["sole_child_age"], 0)

    def test_negative_child_age_is_impossible_not_fabricated(self):
        row = MotherRow(2024, 34, 2, [5, 0], fertyr=2)  # gap five, not negative
        clock, _ = anchor_from_event0(row)
        self.assertEqual(donor_target(clock, -5)["sole_child_age"], 0)
        # A separate gap-two anchor cannot be fabricated into t=-3 support.
        gap_two, _ = anchor_from_event0(MotherRow(2024, 34, 2, [2, 0], fertyr=2))
        self.assertEqual(donor_target(gap_two, -3)["sole_child_age"], -1)

    def test_third_child_is_separate_sensitivity(self):
        row = MotherRow(2024, 35, 3, [8, 5, 1])
        strict, reason = strict_eligibility(row)
        wide, wide_reason = proxy_clock(row, strict=False)
        self.assertIsNone(strict)
        self.assertEqual(reason, "strict sample excludes third/mismatched roster")
        self.assertEqual(wide_reason, "ok")
        self.assertEqual(wide["event_time"], 5)  # second-oldest, not youngest

    def test_twin_or_close_birth_tie_is_ambiguous(self):
        row = MotherRow(2020, 31, 3, [4, 0, 0], fertyr=2)
        clock, reason = proxy_clock(row, strict=False)
        self.assertIsNone(clock)
        self.assertEqual(reason, "ambiguous age tie")

    def test_observed_fertyr_must_be_yes_at_event_zero(self):
        row = MotherRow(2021, 30, 2, [4, 0], fertyr=1)
        clock, reason = strict_eligibility(row)
        self.assertIsNone(clock)
        self.assertEqual(reason, "observed FERTYR is not yes")

    def test_first_birth_twin_has_no_older_child(self):
        row = MotherRow(2020, 27, 2, [0, 0], fertyr=2)
        clock, reason = proxy_clock(row)
        self.assertIsNone(clock)
        self.assertEqual(reason, "no older child plus two linked children")

    def test_age_band_and_nonnegative_age_are_eligibility_gates(self):
        too_young, reason = strict_eligibility(MotherRow(2020, 24, 2, [4, 0], fertyr=2))
        self.assertIsNone(too_young)
        self.assertEqual(reason, "age at proxy birth outside 25-45")
        negative, reason = strict_eligibility(MotherRow(2020, 30, 2, [4, -1], fertyr=2))
        self.assertIsNone(negative)
        self.assertEqual(reason, "negative linked child age")

    def test_nonzero_post_row_cannot_become_donor_anchor(self):
        post = MotherRow(2024, 36, 2, [7, 2])
        anchor, reason = anchor_from_event0(post)
        self.assertIsNone(anchor)
        self.assertEqual(reason, "donor anchor requires event time zero")
        post_clock, _ = proxy_clock(post)
        with self.assertRaises(ValueError):
            donor_target(post_clock, -1)

    def test_momloc_is_scoped_by_household_not_pernum_alone(self):
        mother = {"year": 2020, "sample": 1, "serial": 10, "pernum": 1}
        rows = [
            {"year": 2020, "sample": 1, "serial": 10, "pernum": 2,
             "momloc": 1, "age": 0},
            # Same PERNUM and MOMLOC in a different household: must not link.
            {"year": 2020, "sample": 1, "serial": 11, "pernum": 2,
             "momloc": 1, "age": 4},
        ]
        self.assertEqual([r["age"] for r in linked_children(rows, mother)], [0])


if __name__ == "__main__":
    unittest.main()
