#!/usr/bin/env python3
"""Audit the saved rooms curve's reference support; no regression or target edit.

Reconstruct every primary point from the September 5 cohort coefficients and
the estimator's input-sample cohort shares. Then demonstrate a null
reparameterization: for recent cohorts without K=-2, set their observed K=-1
coefficient to zero and offset their person fixed effects. This is deliberately
NOT a replacement estimator, alternative causal effect, or standard-error claim.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import struct
import time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PRIMARY = ROOT / "code/data/psid_followup_mar2026/output/sa_rooms_first_birth_household_aligned_v1"
REVIEW = ROOT / "output/model/e5f_first_birth_measurement_review_20260905a"
PAIRED = REVIEW / "reference_cluster/full"
OUTPUT = ROOT / "code/data/psid_followup_mar2026/output/first_birth_correction_review"


class FixedWidthStataReader:
    """Read selected scalar cells from local release-118/119 Stata files.

    Avoids a 6-GB shelf scan and a 9-GB raw-file load for a small matched audit.
    The reader was checked against pandas on 144 raw and 240 shelf cells.
    """

    def __init__(self, path):
        self.path = path
        self.file = path.open("rb")
        header = self.file.read(1000)
        release = int(re.search(b"<release>(.*?)</release>", header).group(1))
        assert release in (118, 119)
        self.endian = "<" if b"<byteorder>LSF" in header else ">"
        k_format = "I" if release == 119 else "H"
        k = struct.unpack(self.endian + k_format, re.search(b"<K>(.*?)</K>", header, re.S).group(1))[0]
        self.n = struct.unpack(self.endian + "Q", re.search(b"<N>(.*?)</N>", header, re.S).group(1))[0]
        map_start = header.index(b"<map>") + len(b"<map>")
        positions = struct.unpack(self.endian + "14Q", header[map_start:map_start + 112])
        self.data_start = positions[9] + len(b"<data>")
        self.file.seek(positions[2] + len(b"<variable_types>"))
        types = struct.unpack(self.endian + str(k) + "H", self.file.read(k * 2))
        self.file.seek(positions[3] + len(b"<varnames>"))
        names = self.file.read(k * 129)
        numeric = {65530: (1, "b"), 65529: (2, "h"), 65528: (4, "i"),
                   65527: (4, "f"), 65526: (8, "d"), 32768: (8, None)}
        self.fields = {}
        offset = 0
        for i, kind in enumerate(types):
            name = names[i * 129:(i + 1) * 129].split(b"\0")[0].decode()
            if kind in numeric:
                width, fmt = numeric[kind]
            else:
                assert 1 <= kind <= 2045
                width, fmt = kind, None
            self.fields[name] = (offset, width, fmt)
            offset += width
        self.row_width = offset
        assert positions[10] == self.data_start + self.n * offset + len(b"</data>")

    def value(self, row, variable):
        assert 0 <= row < self.n
        offset, width, fmt = self.fields[variable]
        assert fmt is not None
        self.file.seek(self.data_start + row * self.row_width + offset)
        value = struct.unpack(self.endian + fmt, self.file.read(width))[0]
        thresholds = {"b": 101, "h": 32741, "i": 2147483621, "f": 8e36, "d": 8e307}
        return None if value >= thresholds[fmt] else value

    def first_id_row(self, person_id):
        # The preserved shelf is ordered by ID/year. Each successful result is
        # checked against its preceding ID and its entire local person block.
        low, high = 0, self.n
        while low < high:
            middle = (low + high) // 2
            value = self.value(middle, "ID")
            if value is not None and value < person_id:
                low = middle + 1
            else:
                high = middle
        return low


def check_raw_timing():
    started = time.monotonic()
    data = ROOT.parent / "PSID"
    raw_path = data / "Construction_Files/Data/Users/DD/Dropbox (University of Michigan)/Data/PSID/PSID_CMS/PSID_COMPLETE_MAIN_STUDY_1968_2019.dta"
    shelf_path = data / "PSIDSHELF_MOBILITY.dta"
    raw, shelf = FixedWidthStataReader(raw_path), FixedWidthStataReader(shelf_path)
    raw_variables = {1983: "V8969", 1984: "V10432", 1985: "V11614",
                     2015: "ER60029", 2017: "ER66029", 2019: "ER72029"}
    statistics = defaultdict(Counter)
    evidence_hash = hashlib.sha256()
    matched_ids = 0
    for j in range(512):
        index = round(j * (raw.n - 1) / 511)
        person = raw.value(index, "ID")
        position = shelf.first_id_row(person)
        assert position < shelf.n and shelf.value(position, "ID") == person
        assert position == 0 or shelf.value(position - 1, "ID") < person
        matched_ids += 1
        records = {}
        years = []
        for row in range(position, min(position + 80, shelf.n)):
            if shelf.value(row, "ID") != person:
                break
            year = shelf.value(row, "year")
            years.append(year)
            if year in (1982, 1983, 1984, 1985, 2013, 2015, 2017, 2019):
                records[year] = (shelf.value(row, "ACTUALROOMS_"), shelf.value(row, "CURRENT"))
        if len(years) >= 80:
            raise AssertionError("Person block unexpectedly exceeds 80 rows")
        assert years == sorted(set(years))
        raw_rooms = {year: raw.value(index, variable) for year, variable in raw_variables.items()}
        for year in (1984, 1985, 2017, 2019):
            gap = 1 if year < 1997 else 2
            before, current = records.get(year - gap), records.get(year)
            if not before or not current or current[1] != 1:
                continue
            value, previous_value = raw_rooms[year], raw_rooms[year - gap]
            if value is None or previous_value is None or before[0] is None or current[0] is None:
                continue
            if not (1 <= value <= (8 if year <= 1984 else 20)
                    and 1 <= previous_value <= (8 if year - gap <= 1984 else 20)):
                continue
            counts = statistics[year]
            counts["usable_person_waves"] += 1
            counts["same_label_year_equal_to_raw"] += current[0] == value
            counts["preceding_label_year_equal_to_raw"] += before[0] == value
            if value != previous_value:
                counts["raw_rooms_changed_since_previous_wave"] += 1
                counts["changers_same_label_equal"] += current[0] == value
                counts["changers_preceding_label_equal"] += before[0] == value
            evidence_hash.update(repr((person, year, value, previous_value, before, current)).encode())
        assert time.monotonic() - started < 100, "Bounded timing probe exceeded 100 seconds"
    assert matched_ids == 512
    assert set(statistics) == {1984, 1985, 2017, 2019}
    for counts in statistics.values():
        assert counts["usable_person_waves"] == counts["preceding_label_year_equal_to_raw"]
        assert counts["raw_rooms_changed_since_previous_wave"] > 0
    OUTPUT.mkdir(parents=True, exist_ok=True)
    table = [{"raw_interview_year": year, **dict(statistics[year])} for year in sorted(statistics)]
    write_csv(OUTPUT / "raw_timing_match.csv", table)
    receipt = {
        "scope": "512 deterministic evenly spaced raw person indices; aggregate exports only",
        "sample": "current in destination year; raw rooms 1..8 through 1984 and 1..20 later in both adjacent waves; both shelf room fields observed",
        "matched_person_ids": matched_ids, "comparisons": table,
        "interpretation": "In the checked person-waves, the custom shelf series is labelled one interview early. Moving it forward one observed interview recovers contemporaneous raw rooms.",
        "limitations": "This is a deterministic matched sample in four waves, not a full-panel/all-vintage validation. It does not identify the exact executable that produced the May graph.",
        "selected_evidence_sha256": evidence_hash.hexdigest(),
        "elapsed_seconds": time.monotonic() - started,
        "files": {str(reader.path): {"bytes": reader.path.stat().st_size,
                   "mtime_ns": reader.path.stat().st_mtime_ns, "observations": reader.n,
                   "data_start": reader.data_start, "record_width": reader.row_width,
                   "fields": {v: reader.fields[v] for v in variables}}
                  for reader, variables in ((raw, ["ID", *raw_variables.values()]),
                                            (shelf, ["ID", "year", "CURRENT", "ACTUALROOMS_"]))},
    }
    raw.file.close()
    shelf.file.close()
    (OUTPUT / "raw_timing_verification.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({k: v for k, v in receipt.items() if k != "files"}, indent=2))


def validate_reader():
    """Cross-check the selected-cell reader against the bundled pandas reader."""
    import pandas as pd
    receipt = json.loads((OUTPUT / "raw_timing_verification.json").read_text())
    checked = {}
    for filename, metadata in receipt["files"].items():
        path = Path(filename)
        is_raw = "COMPLETE_MAIN_STUDY" in path.name
        columns = ["ID", "V8969", "V10432", "V11614", "ER66029", "ER72029"] if is_raw else ["ID", "year", "CURRENT", "ACTUALROOMS_"]
        count = 24 if is_raw else 60
        custom = FixedWidthStataReader(path)
        with pd.read_stata(path, columns=columns, convert_categoricals=False, chunksize=count) as reader:
            frame = next(reader)
        for index, row in frame.iterrows():
            for variable in columns:
                actual = custom.value(index, variable)
                expected = row[variable]
                assert (actual is None and pd.isna(expected)) or actual == expected
        custom.file.close()
        checked["raw" if is_raw else "shelf"] = len(frame) * len(columns)
    result = {"status": "pass", "cells_matching_pandas": checked, "builder_sha256": sha(Path(__file__))}
    (OUTPUT / "reader_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


def rows(path):
    with path.open(newline="") as file:
        return list(csv.DictReader(file))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, records):
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(records[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def event(k):
    if k <= -7:
        return "F7event"
    if k >= 11:
        return "L11event"
    return f"F{-k}event" if k < 0 else f"L{k}event"


def support(records):
    weights = defaultdict(float)
    counts = defaultdict(int)
    for row in records:
        if int(row["never_treated"]):
            continue
        key = int(float(row["first_child_year"])), event(int(float(row["K"])))
        weights[key] += float(row["weight_sum"])
        counts[key] += int(row["observations"])
    return weights, counts


def main():
    metadata = json.loads((PRIMARY / "metadata.json").read_text())
    verification = json.loads((REVIEW / "reference_followup_verification.json").read_text())
    expected = {PRIMARY / name: digest for name, digest in metadata["outputs"].items()}
    expected[Path(metadata["source"]["do_file"])] = metadata["source"]["do_file_sha256"]
    for relative, digest in verification["aggregate_evidence"].items():
        if Path(relative).name in {
            "coefficients_original.csv", "input_support.csv",
            "estimation_support_original.csv", "estimation_support_reference.csv",
            "fit_receipt_original.csv", "fit_receipt_reference.csv",
        }:
            expected[ROOT / relative] = digest
    for path, digest in expected.items():
        assert sha(path) == digest, f"Historical input changed: {path}"

    inputs = rows(PAIRED / "input_support.csv")
    fitted = rows(PAIRED / "estimation_support_original.csv")
    assert sum(int(r["observations"]) for r in inputs) == 49872
    assert sum(int(r["observations"]) for r in fitted) == 49457
    assert fitted == rows(PAIRED / "estimation_support_reference.csv")
    weights, _ = support(inputs)
    fitted_weights, fitted_counts = support(fitted)
    coefficients = rows(PAIRED / "coefficients_original.csv")
    b = {(int(r["cohort"]), r["event"]): float(r["coefficient"]) for r in coefficients}
    variance = {(int(r["cohort"]), r["event"]): float(r["variance"]) for r in coefficients}
    cohorts = sorted({g for g, _ in fitted_weights})
    reference = {g for g, e in fitted_weights if e == "F2event"}
    assert reference == {g for g, e in weights if e == "F2event"}
    missing_reference = set(cohorts) - reference
    # These cohorts have only odd calendar-year interviews around childbirth.
    shifted = sorted(g for g in missing_reference if g >= 2000 and g % 2 == 0)
    assert shifted == list(range(2000, 2020, 2))
    assert all(fitted_weights[g, "F1event"] > 0 for g in shifted)
    shifts = {g: -b[g, "F1event"] for g in shifted}

    # Check the exact null identity on every observed final-sample cohort/event
    # cell. Each individual belongs to a fixed first-birth cohort (also asserted
    # by the historical paired-regression builder). Every nonreference row has
    # exactly one included event indicator, including the pooled tails.
    max_prediction_shift = 0.0
    for row in fitted:
        if int(row["never_treated"]):
            continue
        g, k = int(float(row["first_child_year"])), int(float(row["K"]))
        delta = shifts.get(g, 0.0)
        if delta:
            assert k != -2
            assert (g, event(k)) in b
        prediction_shift = delta * (k != -2) - delta
        max_prediction_shift = max(max_prediction_shift, abs(prediction_shift))
    assert max_prediction_shift == 0.0

    cohort_rows = []
    for g in cohorts:
        omitted = [e for gg, e in fitted_weights if gg == g and e != "F2event"
                   and b[g, e] == 0.0 and variance[g, e] == 0.0]
        cohort_rows.append({
            "first_birth_year": g, "observed_minus2": g in reference,
            "observations_minus2": fitted_counts.get((g, "F2event"), 0),
            "observations_minus1": fitted_counts.get((g, "F1event"), 0),
            "supported_zero_and_zero_variance_columns": ";".join(sorted(omitted)),
            "diagnostic_coefficient_shift": shifts.get(g, 0.0),
            "opposite_person_fixed_effect_shift": -shifts.get(g, 0.0),
        })
    assert all("F7event" in r["supported_zero_and_zero_variance_columns"]
               for r in cohort_rows if r["first_birth_year"] in shifted)

    published = {int(float(r["relative_time"])): float(r["b"])
                 for r in rows(PRIMARY / "event_study_estimates.csv")}
    curve = []
    for k in range(-7, 12):
        e = event(k)
        shares = {g: weight for (g, ev), weight in weights.items() if ev == e}
        total = sum(shares.values())
        shares = {g: weight / total for g, weight in shares.items()}
        original = sum(w * b.get((g, e), 0.0) for g, w in shares.items())
        change = sum(w * shifts.get(g, 0.0) for g, w in shares.items())
        curve.append({
            "event_time": k, "published_point": published[k],
            "reconstructed_point": original,
            "same_fit_diagnostic_point": original + change,
            "normalization_only_change": change,
            "share_without_minus2": sum(w for g, w in shares.items() if g in missing_reference),
            "share_recent_even_without_minus2": sum(w for g, w in shares.items() if g in shifts),
            "cohorts_in_input_share": len(shares),
        })
    max_reconstruction_gap = max(abs(r["published_point"] - r["reconstructed_point"]) for r in curve)
    assert max_reconstruction_gap < 1e-6  # published CSV used Stata float storage
    by_time = {r["event_time"]: r for r in curve}
    def teeth(column):
        return sum(by_time[k][column] - (by_time[k-1][column] + by_time[k+1][column])/2
                   for k in (1, 3, 5, 7, 9)) / 5
    def contrast(column):
        return by_time[3][column] - by_time[-1][column]

    OUTPUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUTPUT / "cohort_reference_audit.csv", cohort_rows)
    write_csv(OUTPUT / "event_curve_audit.csv", curve)
    receipt = {
        "scope": "saved aggregate evidence only; no new regression, target, or standard errors",
        "input_observations": 49872, "fitted_observations": 49457,
        "cohorts_without_reference": sorted(missing_reference),
        "diagnostic_shifted_cohorts": shifted,
        "diagnostic_shift_rule": "set each recent even-birth cohort's K=-1 coefficient to zero; offset its person fixed effects",
        "max_primary_curve_reconstruction_gap": max_reconstruction_gap,
        "max_algebraic_prediction_shift_on_final_sample_support": max_prediction_shift,
        "reported_four_year_contrast": contrast("reconstructed_point"),
        "same_fit_diagnostic_four_year_contrast": contrast("same_fit_diagnostic_point"),
        "mean_odd_peak_above_adjacent_even_points_original": teeth("reconstructed_point"),
        "mean_odd_peak_above_adjacent_even_points_diagnostic": teeth("same_fit_diagnostic_point"),
        "limitation": "The alternative curve still mixes K=-2 and K=-1 references. Its purpose is a non-identification witness, not a corrected curve or a bound on bias. Remaining fluctuations are not diagnosed here.",
        "input_hashes": {str(p.relative_to(ROOT)): digest for p, digest in expected.items()},
        "builder_sha256": sha(Path(__file__)),
    }
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.6), gridspec_kw={"width_ratios": [1.3, 1]})
    x = [r["event_time"] for r in curve]
    axes[0].plot(x, [r["reconstructed_point"] for r in curve], "o-", color="#214d70", ms=4, label="Current reported curve")
    axes[0].plot(x, [r["same_fit_diagnostic_point"] for r in curve], "s--", color="#bb5934", ms=3.5, label="Same fit, different arbitrary zeros")
    axes[0].set(title="The annual curve depends on arbitrary reference choices", ylabel="Rooms coefficient", xlabel="Years relative to first birth (tails pooled)", xticks=range(-7,12,2))
    axes[0].axhline(0, color=".7", lw=.7)
    axes[0].axvline(0, color=".7", lw=.7)
    axes[0].legend(frameon=False, fontsize=9, loc="lower right")
    axes[1].bar(x, [100*r["share_without_minus2"] for r in curve], color="#a7b7c5", label="All cohorts missing year -2")
    axes[1].bar(x, [100*r["share_recent_even_without_minus2"] for r in curve], color="#bb5934", label="First births in 2000, 2002, …, 2018")
    axes[1].set(title="Missing reference observations vary by event year", ylabel="Share of coefficient aggregation weight (%)", xlabel="Years relative to first birth (tails pooled)", xticks=range(-7,12,2), ylim=(0,57))
    axes[1].legend(frameon=False, fontsize=8.5, loc="upper right")
    fig.suptitle("Diagnostic only — changing the zeros leaves every fitted observation unchanged", fontsize=13)
    fig.text(.5, .02, "The orange line is not a replacement estimate. Alternative uncertainty has not been computed.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0,.055,1,.93))
    fig.savefig(OUTPUT / "normalization_diagnostic.png", dpi=180)
    plt.close(fig)
    for path, digest in expected.items():
        assert sha(path) == digest, f"Input mutated: {path}"
    (OUTPUT / "verification.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({k: v for k, v in receipt.items() if k != "input_hashes"}, indent=2))


def check_all_waves():
    """Validate every prepared row and every survey wave, with aggregate exports."""
    import mmap
    import numpy as np
    import pandas as pd
    started = time.monotonic()
    data = ROOT.parent/'PSID'
    raw_path = data/'Construction_Files/Data/Users/DD/Dropbox (University of Michigan)/Data/PSID/PSID_CMS/PSID_COMPLETE_MAIN_STUDY_1968_2019.dta'
    shelf_path = data/'PSIDSHELF_MOBILITY.dta'
    prepared_path = Path('/tmp/psid_original_timing_20260912b/analysis_sample.dta')
    waves = list(range(1968,1998))+list(range(1999,2020,2))
    variables = ('V102 V592 V1263 V1966 V2565 V3107 V3521 V3937 V4448 V5362 V5862 V6477 '
                 'V7080 V7671 V8360 V8969 V10432 V11614 V13019 V14122 V15138 V16639 '
                 'V18070 V19370 V20670 V22425 ER2029 ER5028 ER7028 ER10032 ER13037 '
                 'ER17040 ER21039 ER25027 ER36027 ER42028 ER47328 ER53028 ER60029 ER66029 ER72029').split()
    assert len(waves)==len(variables)==41

    def heartbeat(phase):
        elapsed=time.monotonic()-started
        (OUTPUT/'all_wave_heartbeat.json').write_text(json.dumps({'phase':phase,'elapsed_seconds':elapsed})+'\n')
        assert elapsed<600, 'Ten-minute read-only audit budget exhausted'

    def columns(path,names):
        cache_dir=Path('/tmp/psid_all_wave_cache')
        cache_dir.mkdir(mode=0o700,exist_ok=True)
        cached=cache_dir/(path.name+'.npz');meta=cache_dir/(path.name+'.json')
        signature={'size':path.stat().st_size,'mtime_ns':path.stat().st_mtime_ns,'columns':names}
        if cached.exists() and meta.exists() and json.loads(meta.read_text())==signature:
            with np.load(cached) as z:return {n:z[n] for n in names}
        r=FixedWidthStataReader(path)
        mm=mmap.mmap(r.file.fileno(),0,access=mmap.ACCESS_READ)
        types={'b':'i1','h':'i2','i':'i4','f':'f4','d':'f8'}
        thresholds={'b':101,'h':32741,'i':2147483621,'f':8e36,'d':8e307}
        out={n:np.empty(r.n,dtype='f8') for n in names}
        for start in range(0,r.n,1024):
            stop=min(start+1024,r.n)
            for n in names:
                offset,width,fmt=r.fields[n]
                a=np.ndarray((r.n,),dtype=r.endian+types[fmt],buffer=mm,
                             offset=r.data_start+offset,strides=(r.row_width,))
                out[n][start:stop]=a[start:stop]
                del a
            if start%16384==0:heartbeat(path.name)
        for n in names:
            out[n][out[n]>=thresholds[r.fields[n][2]]]=np.nan
        mm.close();r.file.close()
        # Independent pandas parsing of all requested columns on the first four rows.
        with pd.io.stata.StataReader(path,columns=names,convert_categoricals=False) as reader:
            first=reader.read(nrows=4)
        for n in names:
            assert np.array_equal(out[n][:4],first[n].to_numpy(dtype=float),equal_nan=True)
        np.savez(cached,**out);cached.chmod(0o600)
        meta.write_text(json.dumps(signature)+'\n')
        return out

    heartbeat('source crosswalk')
    with pd.io.stata.StataReader(raw_path,convert_categoricals=False) as reader:
        labels=reader.variable_labels()
    ordered=list(labels);crosswalk=[]
    for year,var in zip(waves,variables):
        i=ordered.index(var)
        release=next(j for j in range(i,-1,-1) if labels[ordered[j]].strip()=='RELEASE NUMBER')
        identifier=ordered[release+1];label=labels[identifier]
        assert str(year) in label or re.search(r'(?<!\d)'+str(year)[2:]+r'(?!\d)',label), (year,var,identifier,label)
        crosswalk.append({'year':year,'rooms_variable':var,'rooms_label':labels[var],
                          'family_year_identifier':identifier,'family_year_label':label})
    write_csv(OUTPUT/'all_wave_variable_crosswalk.csv',crosswalk)
    raw=columns(raw_path,['ID']+variables)
    shelf=columns(shelf_path,['ID','year','CURRENT','ACTUALROOMS_'])
    prepared=columns(prepared_path,['ID','year','rooms','rooms_aligned','AGEREP','EDUYEAR','f_c_y','lastcohort'])
    shelf_order=np.lexsort((shelf['year'],shelf['ID']))
    shelf={n:a[shelf_order] for n,a in shelf.items()}
    order=np.argsort(raw['ID']);raw_ids=raw['ID'][order]
    assert np.all(np.diff(raw_ids)>0)
    keys=shelf['ID'].astype('i8')*10000+shelf['year'].astype('i8')
    assert np.all(np.diff(keys)>0)
    donor=np.full(len(keys),np.nan)
    valid=(shelf['ID'][1:]==shelf['ID'][:-1]) & np.isin(shelf['year'][1:]-shelf['year'][:-1],[1,2])
    donor[1:]=np.where(valid,shelf['ACTUALROOMS_'][:-1],np.nan)
    pkeys=prepared['ID'].astype('i8')*10000+prepared['year'].astype('i8')
    spos=np.searchsorted(keys,pkeys)
    assert np.all(spos<len(keys)) and np.array_equal(keys[spos],pkeys)
    equal=lambda a,b:(a==b)|(np.isnan(a)&np.isnan(b))
    original_bad=~equal(prepared['rooms'],shelf['ACTUALROOMS_'][spos])
    adjusted_bad=~equal(prepared['rooms_aligned'],donor[spos])
    rows=[];pairs=[];uncovered=[]
    common=np.isfinite(prepared['rooms']) & np.isfinite(prepared['rooms_aligned']) & np.isfinite(prepared['AGEREP']) & np.isfinite(prepared['EDUYEAR'])
    report_population=[('current_shelf',shelf['ID'],shelf['year'],shelf['ACTUALROOMS_'],donor,shelf['CURRENT']==1),
                       ('prepared_common',prepared['ID'],prepared['year'],prepared['rooms'],prepared['rooms_aligned'],common)]
    for population,ids,years,original,adjusted,eligible in report_population:
        for year in np.unique(years[eligible & ~np.isin(years,waves)]):
            mask=eligible & (years==year)
            uncovered.append({'population':population,'year':int(year),'rows':int(mask.sum()),
                              'adjusted_observed':int(np.isfinite(adjusted[mask]).sum())})
        idx=np.searchsorted(raw_ids,ids)
        idx=np.minimum(idx,len(raw_ids)-1)
        matched=raw_ids[idx]==ids
        idx=order[idx]
        for year,var in zip(waves,variables):
            requested=eligible & (years==year)
            mask=requested & matched
            observed=raw[var][idx[mask]];old=original[mask];new=adjusted[mask]
            oldeq=equal(observed,old);neweq=equal(observed,new)
            both=np.isfinite(observed)&np.isfinite(new)
            substantive=both & (observed>=1) & (observed<=(8 if year<=1984 else 20))
            rows.append({'population':population,'year':year,'rooms_variable':var,'rows':int(mask.sum()),
                         'eligible_rows':int(requested.sum()),'unmatched_source_person':int((requested&~matched).sum()),
                         'original_exact':int(oldeq.sum()),'adjusted_exact':int(neweq.sum()),
                         'raw_missing':int(np.isnan(observed).sum()),'adjusted_missing':int(np.isnan(new).sum()),
                         'both_observed':int(both.sum()),'both_observed_mismatch':int((both&~neweq).sum()),
                         'raw_observed_adjusted_missing':int((np.isfinite(observed)&np.isnan(new)).sum()),
                         'raw_missing_adjusted_observed':int((np.isnan(observed)&np.isfinite(new)).sum()),
                         'substantive_observed':int(substantive.sum()),'substantive_mismatch':int((substantive&~neweq).sum())})
            if np.any(both&~neweq):
                values,counts=np.unique(np.column_stack([observed[both&~neweq],new[both&~neweq]]),axis=0,return_counts=True)
                pairs.extend({'population':population,'year':year,'raw_value':float(v[0]),'adjusted_value':float(v[1]),'count':int(c)} for v,c in zip(values,counts))
        heartbeat(population)
    write_csv(OUTPUT/'all_wave_timing_validation.csv',rows)
    if pairs:write_csv(OUTPUT/'all_wave_mismatch_values.csv',pairs)
    receipt={'status':'completed audit; inspect mismatch counts before claiming validation',
             'waves':waves,'source_rows':len(raw_ids),'shelf_rows':len(keys),'prepared_rows':len(pkeys),
             'prepared_original_construction_mismatches':int(original_bad.sum()),
             'prepared_adjusted_construction_mismatches':int(adjusted_bad.sum()),
             'independent_pandas_cells':4*(42+4+8),'elapsed_seconds':time.monotonic()-started,
             'years_without_source_wave':uncovered,
             'scope':'Every source person, every survey wave, all current shelf rows and all prepared common-sample rows; no value recoding',
             'source_files':{str(p):{'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in [raw_path,shelf_path,prepared_path]},
             'comparison':rows}
    (OUTPUT/'all_wave_validation_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['comparison','source_files']},indent=2))
    print(json.dumps([r for r in rows if r['both_observed_mismatch'] or r['raw_observed_adjusted_missing'] or r['raw_missing_adjusted_observed']],indent=2))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check-raw-timing", action="store_true", help="Run only the bounded direct raw-to-shelf timing comparison")
    mode.add_argument("--validate-reader", action="store_true", help="Cross-check selected-cell reads against pandas; requires pandas")
    mode.add_argument("--check-all-waves", action="store_true", help="Validate every wave and prepared row against source survey columns")
    args = parser.parse_args()
    if args.check_all_waves:
        check_all_waves()
    elif args.check_raw_timing:
        check_raw_timing()
    elif args.validate_reader:
        validate_reader()
    else:
        main()
