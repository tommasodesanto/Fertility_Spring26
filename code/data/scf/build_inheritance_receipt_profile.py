#!/usr/bin/env python3
"""Build an auditable SCF Figure 3 inheritance-receipt profile.

This module intentionally contains no model mapping, interpolation, or data
download.  It extracts the published Figure 3 table and applies an explicitly
diagnostic three-to-four-year timing conversion.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import socket
from html.parser import HTMLParser
from pathlib import Path


SOURCE_URL = (
    "https://www.federalreserve.gov/econres/notes/feds-notes/"
    "how-does-intergenerational-wealth-transmission-affect-wealth-"
    "concentration-accessible-20180601.htm"
)
INCOME_GROUPS = ("bottom50", "middle40", "top10")
AGES = tuple(range(25, 81))
REFERENCE_CELLS = {
    (25, "bottom50"): (0.0263, 11.1566),
    (25, "middle40"): (0.0410, 26.4286),
    (25, "top10"): (0.0540, 283.5551),
    (60, "bottom50"): (0.0420, 113.6140),
    (60, "middle40"): (0.0725, 160.8385),
    (60, "top10"): (0.0919, 420.6898),
    (80, "bottom50"): (0.0186, 71.5786),
    (80, "middle40"): (0.0315, 106.3383),
    (80, "top10"): (0.0483, 389.7972),
}


class Figure3TableParser(HTMLParser):
    """Extract data rows from the first table following ``h5#fig3``."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self._saw_fig3 = False
        self._in_table = False
        self._table_depth = 0
        self._in_cell = False
        self._row: list[str] | None = None
        self._cell_parts: list[str] = []
        self.rows: list[list[str]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = dict(attrs)
        if tag == "h5" and attributes.get("id") == "fig3":
            self._saw_fig3 = True
            return
        if tag == "table":
            if self._in_table:
                self._table_depth += 1
            elif self._saw_fig3:
                self._in_table = True
                self._table_depth = 1
                self._saw_fig3 = False
            return
        if self._in_table and tag == "tr":
            self._row = []
        elif self._in_table and tag in ("td", "th") and self._row is not None:
            self._in_cell = True
            self._cell_parts = []

    def handle_endtag(self, tag: str) -> None:
        if tag == "table" and self._in_table:
            self._table_depth -= 1
            if self._table_depth == 0:
                self._in_table = False
            return
        if not self._in_table:
            return
        if tag in ("td", "th") and self._in_cell and self._row is not None:
            self._row.append("".join(self._cell_parts).strip())
            self._in_cell = False
        elif tag == "tr" and self._row is not None:
            self.rows.append(self._row)
            self._row = None

    def handle_data(self, data: str) -> None:
        if self._in_cell:
            self._cell_parts.append(data)


def _number(value: str, label: str) -> float:
    try:
        number = float(value.strip().replace(",", ""))
    except ValueError as exc:
        raise ValueError(f"invalid {label}: {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"non-finite {label}: {value!r}")
    return number


def _validate_record(record: dict[str, object]) -> None:
    probability = float(record["probability_3y"])
    amount = float(record["conditional_amount_3y"])
    if not math.isfinite(probability) or not math.isfinite(amount):
        raise ValueError("probability_3y and conditional_amount_3y must be finite")
    if not 0.0 <= probability <= 1.0:
        raise ValueError(f"probability_3y outside [0, 1]: {probability}")
    if amount < 0.0:
        raise ValueError(f"conditional_amount_3y is negative: {amount}")


def parse_figure3_html(html: str) -> list[dict[str, object]]:
    """Parse Figure 3 rows without imposing the expected age/group support."""
    parser = Figure3TableParser()
    parser.feed(html)
    parser.close()
    if not parser.rows:
        raise ValueError("could not find the table following h5#fig3")

    records: list[dict[str, object]] = []
    seen: set[tuple[int, str]] = set()
    for row in parser.rows:
        # The published data rows have age, three probabilities, a blank cell,
        # and three conditional amounts. Header rows are ignored.
        if len(row) == 7:
            probability_cells, amount_cells = row[1:4], row[4:7]
        elif len(row) == 8 and not row[4]:
            probability_cells, amount_cells = row[1:4], row[5:8]
        else:
            continue
        try:
            age = int(row[0].strip())
        except ValueError:
            continue
        probabilities = [_number(value, "probability") for value in probability_cells]
        amounts = [_number(value, "conditional amount") for value in amount_cells]
        for group, probability, amount in zip(INCOME_GROUPS, probabilities, amounts):
            record: dict[str, object] = {
                "age": age,
                "income_group": group,
                "probability_3y": probability,
                "conditional_amount_3y": amount,
            }
            _validate_record(record)
            key = (age, group)
            if key in seen:
                raise ValueError(f"duplicate age/income-group cell: {key}")
            seen.add(key)
            records.append(record)
    if not records:
        raise ValueError("Figure 3 table contains no parseable data rows")
    return records


def validate_complete_support(records: list[dict[str, object]]) -> None:
    """Require the complete published 25--80 by income-group table."""
    expected = {(age, group) for age in AGES for group in INCOME_GROUPS}
    observed: set[tuple[int, str]] = set()
    for record in records:
        _validate_record(record)
        key = (int(record["age"]), str(record["income_group"]))
        if key in observed:
            raise ValueError(f"duplicate age/income-group cell: {key}")
        observed.add(key)
    if observed != expected:
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        raise ValueError(f"incomplete Figure 3 support; missing={missing}, extra={extra}")
    lookup = {(int(r["age"]), str(r["income_group"])): r for r in records}
    for key, (expected_p, expected_amount) in REFERENCE_CELLS.items():
        record = lookup[key]
        if record["probability_3y"] != expected_p or record["conditional_amount_3y"] != expected_amount:
            raise ValueError(f"reference cell differs from published value: {key}")


def map_three_to_four_years(records: list[dict[str, object]]) -> list[dict[str, object]]:
    """Map prior-three-year receipt profiles to a four-year Poisson diagnostic."""
    mapped: list[dict[str, object]] = []
    for record in records:
        _validate_record(record)
        probability = float(record["probability_3y"])
        amount = float(record["conditional_amount_3y"])
        mean_four_year = (4.0 / 3.0) * probability * amount
        if probability == 0.0:
            probability_four_year = 0.0
            conditional_amount_four_year = 0.0
        elif probability == 1.0:
            probability_four_year = 1.0
            conditional_amount_four_year = (4.0 / 3.0) * amount
        else:
            probability_four_year = -math.expm1((4.0 / 3.0) * math.log1p(-probability))
            conditional_amount_four_year = mean_four_year / probability_four_year
        mapped.append(
            {
                **record,
                "probability_4y": probability_four_year,
                "mean_amount_4y": mean_four_year,
                "conditional_amount_4y": conditional_amount_four_year,
            }
        )
    return mapped


def _write_csv(path: Path, records: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def build_profile(source_html: Path, output_dir: Path, overwrite: bool = False) -> dict[str, Path]:
    """Build the two CSVs and their JSON source receipt from local HTML."""
    html_bytes = source_html.read_bytes()
    records = parse_figure3_html(html_bytes.decode("utf-8"))
    validate_complete_support(records)
    mapped = map_three_to_four_years(records)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "published_csv": output_dir / "inheritance_receipt_profile_published_3y.csv",
        "mapped_csv": output_dir / "inheritance_receipt_profile_mapped_4y.csv",
        "receipt_json": output_dir / "inheritance_receipt_profile_receipt.json",
    }
    existing = [path for path in outputs.values() if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(f"refusing to overwrite existing output(s): {existing}")
    _write_csv(outputs["published_csv"], records, ["age", "income_group", "probability_3y", "conditional_amount_3y"])
    _write_csv(outputs["mapped_csv"], mapped, ["age", "income_group", "probability_3y", "conditional_amount_3y", "probability_4y", "mean_amount_4y", "conditional_amount_4y"])
    receipt = {
        "source_url": SOURCE_URL,
        "source_sha256": hashlib.sha256(html_bytes).hexdigest(),
        "source_method": "Parsed the first HTML table after h5#fig3; Figure 4 living gifts are excluded.",
        "published_measure": "probability of any inheritance during the prior three years and total receipt conditional on receipt, pooled SCF 1995--2016, by usual-income percentile within age.",
        "amount_unit_caveat": "The accessible HTML amount header says 'Fraction'. Numeric amounts are retained as published but treated only as relative amount weights; they are not labeled dollars.",
        "mapping_assumption": "Diagnostic constant-Poisson arrival intensity within each age/income cell: p4=-expm1((4/3)*log1p(-p3)); mean4=(4/3)*p3*amount3; conditional_amount4=mean4/p4, with explicit p3=0 and p3=1 limits.",
        "undefined_tails": "No interpolation or extrapolation outside ages 25--80; no age-18 or age-82 assumption.",
        "uncertainty": "Unavailable in the published HTML table.",
        "income_concept_caveat": "Usual-income percentile within age does not directly match the model's labor-earnings income states.",
        "record_count": len(records),
    }
    with outputs["receipt_json"].open("w", encoding="utf-8") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return outputs


def _torch_or_slurm() -> bool:
    return bool(os.environ.get("SLURM_JOB_ID")) or "torch" in socket.gethostname().lower()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-html", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if not _torch_or_slurm():
        parser.error("CLI is restricted to Torch or a SLURM allocation")
    build_profile(args.source_html, args.output_dir, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
