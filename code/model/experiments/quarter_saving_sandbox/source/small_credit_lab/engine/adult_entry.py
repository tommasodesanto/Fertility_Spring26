"""Extracted from code/model/intergen_eqscale_seq_optimized/adult_entry.py (sha256 1e16e29cec889a78e651572a586cec3c4384cc72a1e945de4b68e0cd8f8b5a1e).

Mechanical copy by refactor_lab/materialize.py: only reachable top-level
definitions, bodies byte-identical; import edits listed in the receipt.
"""
from __future__ import annotations
from dataclasses import dataclass
import math


REPLACEMENT_FERTILITY = 2.1


def adjusted_births(raw_births: float, top_bin_entries: float,
                    top_bin_weight: float, top_state: int = 3) -> float:
    """Return child units after adding the excess represented by the top bin."""
    raw = float(raw_births)
    top = float(top_bin_entries)
    weight = float(top_bin_weight)
    if not all(map(math.isfinite, (raw, top, weight))):
        raise ValueError("Birth accounting inputs must be finite")
    if raw < 0 or top < 0 or top > raw + 1e-12 or weight < top_state:
        raise ValueError("Invalid birth or top-bin accounting inputs")
    return raw + (weight - top_state) * top


def potential_entry_households(adjusted_birth_children: float) -> float:
    """Apply the retained birth-to-household conversion exactly once."""
    births = float(adjusted_birth_children)
    if not math.isfinite(births) or births < 0:
        raise ValueError("Adjusted births must be finite and nonnegative")
    return births / REPLACEMENT_FERTILITY
