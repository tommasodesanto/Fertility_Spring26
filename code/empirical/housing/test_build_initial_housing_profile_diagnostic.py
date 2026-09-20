from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_initial_housing_profile_diagnostic as d


def _fixture():
    return pd.DataFrame([
        dict(year=2005, sample=200501, met2013=1, gq=1, pernum=1, relate=1, hhwt=2., age=30, ownershp=1, rooms=10., unitsstr=3, nchild=3, yngch=5, eldch=2),
        dict(year=2005, sample=200501, met2013=2, gq=1, pernum=1, relate=1, hhwt=1., age=19, ownershp=2, rooms=4., unitsstr=2, nchild=1, yngch=99, eldch=99),
        dict(year=2005, sample=200501, met2013=1, gq=1, pernum=1, relate=1, hhwt=3., age=30, ownershp=1, rooms=8., unitsstr=10, nchild=0, yngch=99, eldch=99),
        dict(year=2005, sample=200501, met2013=1, gq=1, pernum=1, relate=1, hhwt=4., age=31, ownershp=2, rooms=6., unitsstr=2, nchild=3, yngch=5, eldch=10),
        dict(year=2005, sample=200501, met2013=1, gq=1, pernum=1, relate=1, hhwt=5., age=32, ownershp=2, rooms=5., unitsstr=2, nchild=2, yngch=5, eldch=10),
    ])


def test_active_national_split_and_due_filter():
    acc = {}
    d._aggregate_chunk(_fixture(), {1}, acc)
    active = sum(v["hhwt"] for k, v in acc.items() if k[0] == "active42" and k[1] == "all_structures" and k[2] == "annual")
    national = sum(v["hhwt"] for k, v in acc.items() if k[0] == "national" and k[1] == "all_structures" and k[2] == "annual")
    assert active == 14.0
    assert national == 15.0
    due = sum(v["hhwt"] for k, v in acc.items() if k[0] == "active42" and k[1] == "DUE" and k[2] == "annual")
    assert due == 5.0


def test_cap_before_aggregate_and_current_child_bins():
    acc = {}
    d._aggregate_chunk(_fixture(), {1}, acc)
    capped = sum(v["rooms_capped9_sum"] for k, v in acc.items() if k[0] == "active42" and k[1] == "all_structures" and k[2] == "annual")
    uncapped = sum(v["rooms_sum"] for k, v in acc.items() if k[0] == "active42" and k[1] == "all_structures" and k[2] == "annual")
    assert capped == 91.0
    assert uncapped == 93.0
    assert any(k[5] == "3+" and k[6] == "with_young_child" for k in acc)
    assert any(k[5] == "0" and k[6] == "without_young_child" for k in acc)


def test_four_year_age_bins_start_at_18():
    acc = {}
    d._aggregate_chunk(_fixture(), {1}, acc)
    assert any(k[2] == "four_year" and k[3] == 18 for k in acc)
    assert any(k[2] == "four_year" and k[3] == 30 for k in acc)


def test_exact_scalar_groups_do_not_use_proxy_definitions():
    acc = {}
    stats = {scope: {name: {"hhwt": 0.0, "rooms_capped9_sum": 0.0, "owner_hhwt": 0.0}
                     for name in ("mean_rooms", "ownership", "newparent", "nochild", "family_large", "family_small")}
                     for scope in ("active42", "national")}
    d._aggregate_chunk(_fixture(), {1}, acc, stats)
    # The non-DUE 3-child row enters the family-room contrast but not ownership.
    assert stats["active42"]["family_large"]["hhwt"] == 6.0
    assert stats["active42"]["family_small"]["hhwt"] == 5.0
    # ELDCH<4 identifies the recent parent even when YNGCH=99; the adult-oldest
    # row with a young child is excluded from recent-parent status.
    assert stats["active42"]["newparent"]["hhwt"] == 2.0
    assert stats["active42"]["ownership"]["hhwt"] == 5.0
