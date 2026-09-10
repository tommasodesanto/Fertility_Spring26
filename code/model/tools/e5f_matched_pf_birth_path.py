"""Calendar-aligned aggregate birth-flow diagnostics; never an empirical TFR fit.

The four-year decision at t supplies births in t+1,...,t+4. Historical
comparisons therefore use decisions 2007,2011,2015,2019, not the 2023 choice.
No parameters, observation weights or production moments are changed here.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

HISTORICAL_DECISIONS = (2007, 2011, 2015, 2019)


def number(row, name):
    value = float(row[name])
    if not math.isfinite(value):
        raise ValueError(f'Nonfinite {name}')
    return value


def read_rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def compare_birth_path(model_rows, empirical_rows, *, anchor_first_block_births=None):
    """Return four block comparisons plus transparent diagnostic scores.

    Each candidate's own first block normalizes the shape index. The optional
    common anchor is the unchanged candidate's first-block birth flow, allowing
    its announcement/level movement to remain visible across alternatives.
    Scores are unweighted mean squared index gaps, not SMM losses. Missing,
    duplicate, wrongly dated or incomplete sets of empirical blocks are rejected.
    The source builder must separately verify all four annual observations.
    """
    model, empirical = {}, {}
    for row in model_rows:
        year_value = number(row, 'calendar_year')
        year = int(year_value)
        if year_value != year or year in model:
            raise ValueError('Noninteger or duplicate model date')
        model[year] = row
    for row in empirical_rows:
        year_value = number(row, 'decision_year')
        year = int(year_value)
        if year_value != year or year in empirical:
            raise ValueError('Noninteger or duplicate empirical decision date')
        if (number(row, 'birth_year_start') != year + 1
                or number(row, 'birth_year_end') != year + 4):
            raise ValueError('Empirical birth years must match t+1,...,t+4')
        empirical[year] = row
    if set(empirical) != set(HISTORICAL_DECISIONS):
        raise ValueError('Exactly the four observed historical blocks are required')
    if not set(HISTORICAL_DECISIONS).issubset(model):
        raise ValueError('Missing historical model decision')
    first = number(model[2007], 'birth_children_topcode_adjusted')
    first_data = number(empirical[2007], 'live_births_total')
    anchor = first if anchor_first_block_births is None else float(anchor_first_block_births)
    if not all(math.isfinite(v) and v > 0 for v in (first, first_data, anchor)):
        raise ValueError('Initial birth counts and common anchor must be positive')
    rows = []
    for year in HISTORICAL_DECISIONS:
        m, e = model[year], empirical[year]
        raw = number(m, 'birth_children')
        adjusted = number(m, 'birth_children_topcode_adjusted')
        heads = number(m, 'adult_population')
        data = number(e, 'live_births_total')
        if min(raw, adjusted, heads, data) <= 0 or adjusted + 1e-12 < raw:
            raise ValueError('Nonpositive birth/head count or invalid top-bin addition')
        if m.get('household_heads') not in ('', None):
            other = number(m, 'household_heads')
            if not math.isclose(other, heads, rel_tol=0, abs_tol=2e-9):
                raise ValueError('Model household mass columns disagree')
        empirical_index = data / first_data
        own_index, common_index = adjusted / first, adjusted / anchor
        rows.append(dict(decision_year=year, birth_year_start=year+1, birth_year_end=year+4,
            observed_live_births=data, observed_birth_index=empirical_index,
            model_raw_births=raw, model_adjusted_births=adjusted,
            model_top_bin_addition=adjusted-raw, model_start_households=heads,
            model_shape_index=own_index, shape_gap=own_index-empirical_index,
            model_common_anchor_index=common_index,
            common_anchor_gap=common_index-empirical_index,
            model_annualized_births_per_start_household=adjusted/(4*heads)))
    score = sum(row['shape_gap']**2 for row in rows[1:]) / 3
    common_score = sum(row['common_anchor_gap']**2 for row in rows) / 4
    return dict(status='calendar_aligned_aggregate_birth_shape_diagnostic', rows=rows,
        shape_mean_squared_gap=score, informative_shape_blocks=3,
        common_anchor_mean_squared_gap=common_score,
        common_anchor_supplied=anchor_first_block_births is not None,
        first_block_births=first, first_block_change_from_anchor=first/anchor-1,
        female_tfr_comparable=False, estimated_preference_path=False,
        production_promoted=False,
        scope='Four-year aggregate birth shape; inherited parameters and explicit source/price provenance required',
        caveats=[
            'Own-block normalization discards the initial birth level; display the first-block change separately.',
            'The common anchor is a diagnostic scale from the unchanged candidate, not an estimated biological conversion.',
            'Top-bin additional children are imputed at entry into parity 3+, not separately timed births.',
            'Births per starting household are not female TFR or annual household-exposure rates.',
            'Historical household totals and age margins are externally conditioned; their fit is not an independent success.',
            'The 2023 decision supplies 2024–2027 births and is excluded from observed 2008–2023 comparisons.',
            'National empirical births and the model housing calibration geography remain a maintained approximation.',
        ])
