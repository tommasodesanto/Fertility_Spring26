#!/usr/bin/env python3
"""Bounded, read-only check of whether shelf FID and AGEREP are contemporaneous
with their row-labelled survey wave, or shifted to the next wave the way
ACTUALROOMS_ was already found to be shifted.

Reader class copied/adapted from FixedWidthStataReader in
code/data/psid_followup_mar2026/review_first_birth_event_study_corrections.py
(the class used by that script's --check-raw-timing mode), extended with a
labels() method so raw variable identities are confirmed from the file's own
variable-label dictionary rather than trusted from any hard-coded name. Never
loads either .dta fully: only the small header/dictionary blocks, plus one
full row per sampled person (raw file) or a short forward scan of a person's
row block (shelf file, sorted by ID/year).

No economic interpretation. Aggregate counts only; two small person examples
with masked IDs.
"""
from __future__ import annotations

import json
import re
import struct
import time
from collections import defaultdict
from pathlib import Path

SHELF_PATH = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta")
# The path named in the task text and hard-coded in review_first_birth_event_study_corrections.py
# (Fertility_Spring26/../PSID/Construction_Files/Data/Users/DD/Dropbox (University of
# Michigan)/Data/PSID/PSID_CMS/PSID_COMPLETE_MAIN_STUDY_1968_2019.dta) no longer resolves --
# only an unopened 9.58GB zip sits at the Construction_Files/Data level now. Located the live
# file via `mdfind` under a *different* project directory with the identical relative subpath
# and a matching size (9,579,079,729 bytes ~= the task's "~9.6 GB"). Treated as the same file,
# just relocated; flagged as an anomaly in the report.
RAW_PATH = Path(
    "/Users/tommasodesanto/Desktop/Projects/Datasets/PSID-SHELF/Construction_Files/Data/Users/DD/"
    "Dropbox (University of Michigan)/Data/PSID/PSID_CMS/PSID_COMPLETE_MAIN_STUDY_1968_2019.dta"
)
OUT_DIR = Path(__file__).resolve().parent


class FixedWidthStataReader:
    """Read selected scalar cells from local release-118/119 Stata files."""

    def __init__(self, path):
        self.path = Path(path)
        self.file = self.path.open("rb")
        header = self.file.read(1000)
        release = int(re.search(b"<release>(.*?)</release>", header).group(1))
        assert release in (118, 119)
        self.release = release
        self.endian = "<" if b"<byteorder>LSF" in header else ">"
        k_format = "I" if release == 119 else "H"
        self.k = struct.unpack(self.endian + k_format, re.search(b"<K>(.*?)</K>", header, re.S).group(1))[0]
        self.n = struct.unpack(self.endian + "Q", re.search(b"<N>(.*?)</N>", header, re.S).group(1))[0]
        map_start = header.index(b"<map>") + len(b"<map>")
        self.positions = struct.unpack(self.endian + "14Q", header[map_start:map_start + 112])
        positions = self.positions
        self.data_start = positions[9] + len(b"<data>")
        self.file.seek(positions[2] + len(b"<variable_types>"))
        types = struct.unpack(self.endian + str(self.k) + "H", self.file.read(self.k * 2))
        self.file.seek(positions[3] + len(b"<varnames>"))
        names = self.file.read(self.k * 129)
        self.varnames = [names[i * 129:(i + 1) * 129].split(b"\0")[0].decode() for i in range(self.k)]
        numeric = {65530: (1, "b"), 65529: (2, "h"), 65528: (4, "i"),
                   65527: (4, "f"), 65526: (8, "d"), 32768: (8, None)}
        self.fields = {}
        offset = 0
        for i, kind in enumerate(types):
            name = self.varnames[i]
            if kind in numeric:
                width, fmt = numeric[kind]
            else:
                assert 1 <= kind <= 2045
                width, fmt = kind, None
            self.fields[name] = (offset, width, fmt)
            offset += width
        self.row_width = offset
        assert positions[10] == self.data_start + self.n * offset + len(b"</data>")
        self._labels = None

    def value(self, row, variable):
        assert 0 <= row < self.n
        offset, width, fmt = self.fields[variable]
        assert fmt is not None
        self.file.seek(self.data_start + row * self.row_width + offset)
        raw = struct.unpack(self.endian + fmt, self.file.read(width))[0]
        return self._nan(raw, fmt)

    @staticmethod
    def _nan(value, fmt):
        thresholds = {"b": 101, "h": 32741, "i": 2147483621, "f": 8e36, "d": 8e307}
        return None if value >= thresholds[fmt] else value

    def read_row(self, row, variables):
        """Read a whole row once and decode several fields from the buffer
        (avoids one seek+read per field; used for the shelf forward scan)."""
        assert 0 <= row < self.n
        self.file.seek(self.data_start + row * self.row_width)
        buf = self.file.read(self.row_width)
        out = {}
        for v in variables:
            offset, width, fmt = self.fields[v]
            raw = struct.unpack(self.endian + fmt, buf[offset:offset + width])[0]
            out[v] = self._nan(raw, fmt)
        return out

    def first_id_row(self, person_id):
        low, high = 0, self.n
        while low < high:
            middle = (low + high) // 2
            value = self.value(middle, "ID")
            if value is not None and value < person_id:
                low = middle + 1
            else:
                high = middle
        return low

    def labels(self):
        if self._labels is None:
            self.file.seek(self.positions[7] + len(b"<variable_labels>"))
            raw = self.file.read(self.k * 321)
            self._labels = {
                name: raw[i * 321:(i + 1) * 321].split(b"\0")[0].decode("latin-1")
                for i, name in enumerate(self.varnames)
            }
        return self._labels

    def close(self):
        self.file.close()


# --- Raw variable identities, confirmed against the file's own variable_labels() dict below ---
RAW_ROOMS = {1983: "V8969", 1984: "V10432", 1985: "V11614", 1986: "V13019",
             2013: "ER53028", 2015: "ER60029", 2017: "ER66029", 2019: "ER72029"}
RAW_FID = {1984: "V10002", 1985: "V11102", 1986: "V12502", 2017: "ER66002", 2019: "ER72002"}
RAW_AGE = {1984: "ER30432", 1985: "ER30466", 1986: "ER30501", 2017: "ER34504", 2019: "ER34704"}
NEXT_WAVE = {1984: 1985, 1985: 1986, 2017: 2019, 2019: None}
WAVES = (1984, 1985, 2017, 2019)


def room_cap(year):
    return 8 if year <= 1984 else 20


def verify_raw_labels(raw):
    """Cross-check every raw variable name used below against the file's own
    variable-label dictionary; fail loudly rather than silently trusting a name."""
    labels = raw.labels()
    checks = {}
    for year, var in RAW_FID.items():
        label = labels[var]
        assert str(year) in label, (year, var, label)
        checks[var] = label
    for year, var in RAW_AGE.items():
        label = labels[var]
        short = f"{year % 100:02d}"
        assert re.search(r"AGE OF INDIVIDUAL", label, re.I) and label.strip().endswith(short), (year, var, label)
        checks[var] = label
    for year, var in RAW_ROOMS.items():
        checks[var] = labels[var]  # rooms label doesn't carry the year; trusted from the
        # existing, already-established audit (review_first_birth_event_study_corrections.py
        # raw_variables dict / check_all_waves crosswalk), reported here for transparency only.
    return checks


def verify_shelf_labels(shelf):
    labels = shelf.labels()
    used = ["FID", "HHID", "AGEREP", "ACTUALROOMS_", "CURRENT", "ID", "year"]
    return {v: labels.get(v, "") for v in used}


def sampled_persons(raw):
    """Same construction as check_raw_timing(): 512 deterministic, evenly spaced
    raw row indices -> person IDs."""
    return [(j, round(j * (raw.n - 1) / 511)) for j in range(512)]


def mask(person_id, seen={}):
    """Stable but non-reversible per-run label for a person ID (report only; never the ID)."""
    if person_id not in seen:
        seen[person_id] = f"person_{chr(65 + len(seen) % 26)}{len(seen) // 26 or ''}"
    return seen[person_id]


def main():
    started = time.monotonic()
    raw = FixedWidthStataReader(RAW_PATH)
    shelf = FixedWidthStataReader(SHELF_PATH)
    assert raw.n == 82573
    raw_label_check = verify_raw_labels(raw)
    shelf_label_check = verify_shelf_labels(shelf)

    shelf_row_vars = ["ID", "year", "CURRENT", "FID", "AGEREP", "ACTUALROOMS_"]

    # counts[wave][variable] -> Counter-like dict
    counts = {w: {v: defaultdict(int) for v in ("FID", "AGEREP", "rooms")} for w in WAVES}
    matched_ids = 0
    annual_candidates = []    # (decisiveness, j, person, records, raw_r)
    biennial_candidates = []
    anomalies = []

    for j, index in sampled_persons(raw):
        person = raw.value(index, "ID")
        assert person is not None
        position = shelf.first_id_row(person)
        if not (position < shelf.n and shelf.read_row(position, ["ID"])["ID"] == person):
            anomalies.append({"person_index": j, "issue": "no shelf row found for sampled raw person"})
            continue
        matched_ids += 1

        records, years_seen = {}, []
        row = position
        while row < shelf.n:
            data = shelf.read_row(row, shelf_row_vars)
            if data["ID"] != person:
                break
            year = int(data["year"])
            years_seen.append(year)
            records[year] = data
            row += 1
            if row - position > 120:
                anomalies.append({"person_index": j, "issue": "person block exceeded 120 rows; scan stopped"})
                break
        assert years_seen == sorted(set(years_seen)), (person, years_seen)

        raw_full = raw.read_row(index, ["ID"] + list(RAW_ROOMS.values()) + list(RAW_FID.values()) + list(RAW_AGE.values()))
        assert raw_full["ID"] == person

        for wave in WAVES:
            current = records.get(wave)
            if current is None or current["CURRENT"] != 1:
                continue
            nxt_wave = NEXT_WAVE[wave]

            # --- rooms (anchors / reproduces the earlier finding) ---
            shelf_rooms = current["ACTUALROOMS_"]
            raw_rooms_same = raw_full.get(RAW_ROOMS.get(wave))
            raw_rooms_next = raw_full.get(RAW_ROOMS.get(nxt_wave)) if nxt_wave else None
            _tally(counts[wave]["rooms"], shelf_rooms, raw_rooms_same, raw_rooms_next,
                   nxt_wave, valid=lambda v, y=wave: v is not None and 1 <= v <= room_cap(y),
                   valid_next=(lambda v: v is not None and 1 <= v <= room_cap(nxt_wave)) if nxt_wave else None)

            # --- FID (family interview number) ---
            shelf_fid = current["FID"]
            raw_fid_same = raw_full.get(RAW_FID.get(wave))
            raw_fid_next = raw_full.get(RAW_FID.get(nxt_wave)) if nxt_wave else None
            _tally(counts[wave]["FID"], shelf_fid, raw_fid_same, raw_fid_next, nxt_wave,
                   valid=lambda v: v is not None and v > 0,
                   valid_next=(lambda v: v is not None and v > 0) if nxt_wave else None)

            # --- AGEREP ---
            shelf_age = current["AGEREP"]
            raw_age_same = raw_full.get(RAW_AGE.get(wave))
            raw_age_next = raw_full.get(RAW_AGE.get(nxt_wave)) if nxt_wave else None
            _tally(counts[wave]["AGEREP"], shelf_age, raw_age_same, raw_age_next, nxt_wave,
                   valid=lambda v: v is not None and 0 <= v <= 120,
                   valid_next=(lambda v: v is not None and 0 <= v <= 120) if nxt_wave else None)

        # --- Task 2 candidate scan (reuse the same 512-person sample; collect every
        # qualifying candidate and pick the clearest afterwards, rather than the first) ---
        for years, bucket in (((1983, 1984, 1985, 1986), annual_candidates),
                               ((2013, 2015, 2017, 2019), biennial_candidates)):
            block = {y: records[y]["ACTUALROOMS_"] for y in years if y in records
                     and records[y]["CURRENT"] == 1 and records[y]["ACTUALROOMS_"] is not None
                     and 1 <= records[y]["ACTUALROOMS_"] <= 20}
            if len(block) != 4:
                continue
            raw_r = {y: raw_full.get(RAW_ROOMS[y]) for y in years}
            if not all(v is not None and 1 <= v <= room_cap(y) for y, v in raw_r.items()):
                continue
            transitions = sum(raw_r[years[i]] != raw_r[years[i + 1]] for i in range(len(years) - 1))
            if transitions == 0:
                continue
            bucket.append((transitions, j, person, years, raw_r, records))

        assert time.monotonic() - started < 180, "Bounded read-only probe exceeded 180 seconds"

    assert matched_ids >= 500  # a handful of edge misses is fine; a mass failure is not

    # Pick the clearest annual candidate (most visible raw-rooms transitions, ties broken by
    # earliest sample index for reproducibility), then the clearest biennial candidate that is
    # a different masked person where possible, so the two example tables are not redundant.
    annual_candidates.sort(key=lambda c: (-c[0], c[1]))
    biennial_candidates.sort(key=lambda c: (-c[0], c[1]))
    task2_annual = task2_biennial = None
    if annual_candidates:
        transitions, j, person, years, raw_r, records = annual_candidates[0]
        task2_annual = _example_table(person, records, years, raw_r)
    other_biennial = [c for c in biennial_candidates if c[2] != (annual_candidates[0][2] if annual_candidates else None)]
    pick_from = other_biennial or biennial_candidates
    if pick_from:
        transitions, j, person, years, raw_r, records = pick_from[0]
        task2_biennial = _example_table(person, records, years, raw_r)

    task1 = {}
    for wave in WAVES:
        task1[str(wave)] = {"next_wave": NEXT_WAVE[wave],
                             "variables": {v: dict(counts[wave][v]) for v in ("FID", "AGEREP", "rooms")}}

    receipt = {
        "scope": "512 deterministic evenly spaced raw person indices (same construction as "
                 "review_first_birth_event_study_corrections.py --check-raw-timing); aggregate "
                 "counts and two masked example tables only",
        "matched_person_ids": matched_ids,
        "waves": WAVES,
        "next_wave_map": NEXT_WAVE,
        "usability_gate": "CURRENT==1 on the shelf row labelled year=wave; shelf value observed; "
                           "raw same-wave value observed (and raw next-wave value observed when a "
                           "next wave exists); rooms additionally restricted to 1..8 (year<=1984) or "
                           "1..20 (year>1984) on the raw side, mirroring the established rooms audit "
                           "and the project's 0/99-are-not-room-counts convention; FID>0; 0<=AGEREP<=120",
        "raw_variables_used": {"family_interview_number": RAW_FID, "age": RAW_AGE, "rooms": RAW_ROOMS},
        "raw_label_check": raw_label_check,
        "shelf_label_check": shelf_label_check,
        "task1_counts": task1,
        "task2_examples": {"annual_1983_1986": task2_annual, "biennial_2013_2019": task2_biennial},
        "anomalies": anomalies,
        "elapsed_seconds": time.monotonic() - started,
        "files": {
            str(RAW_PATH): {"bytes": RAW_PATH.stat().st_size, "mtime_ns": RAW_PATH.stat().st_mtime_ns,
                             "observations": raw.n, "note": "located via mdfind; path differs from "
                             "the task text / review_first_birth_event_study_corrections.py hard-code"},
            str(SHELF_PATH): {"bytes": SHELF_PATH.stat().st_size, "mtime_ns": SHELF_PATH.stat().st_mtime_ns,
                               "observations": shelf.n},
        },
    }
    raw.close()
    shelf.close()
    (OUT_DIR / "row_year_check.json").write_text(json.dumps(receipt, indent=2, default=str) + "\n")
    print(json.dumps({k: v for k, v in receipt.items() if k != "files"}, indent=2, default=str))


def _tally(counter, shelf_val, raw_same, raw_next, nxt_wave, valid, valid_next):
    if shelf_val is None or not valid(raw_same):
        return
    if nxt_wave is not None:
        if not valid_next(raw_next):
            return
        same = shelf_val == raw_same
        nxt = shelf_val == raw_next
        counter["usable_persons"] += 1
        counter["same_wave_matches"] += int(same)
        counter["next_wave_matches"] += int(nxt)
        counter["both"] += int(same and nxt)
        counter["neither"] += int((not same) and (not nxt))
        if raw_same != raw_next:
            counter["changers_n"] += 1
            counter["changers_same_wave_matches"] += int(same)
            counter["changers_next_wave_matches"] += int(nxt)
    else:
        same = shelf_val == raw_same
        counter["usable_persons"] += 1
        counter["same_wave_matches"] += int(same)
        counter["next_wave_matches"] = None
        counter["both"] = None
        counter["neither"] += int(not same)


def _example_table(person, records, years, raw_rooms):
    label = mask(person)
    rows = []
    sorted_years = sorted(records)
    for y in years:
        idx_in_block = sorted_years.index(y)
        preceding_year = sorted_years[idx_in_block - 1] if idx_in_block > 0 else None
        preceding_val = records[preceding_year]["ACTUALROOMS_"] if preceding_year is not None else None
        rows.append({
            "year_label": y,
            "raw_rooms_official_variable": raw_rooms[y],
            "shelf_ACTUALROOMS__on_row_labelled_this_year": records[y]["ACTUALROOMS_"],
            "preceding_row_year_label": preceding_year,
            "shelf_ACTUALROOMS__on_preceding_row": preceding_val,
        })
    return {"person": label, "rows": rows}


if __name__ == "__main__":
    main()
