"""Focused pure-fixture checks. Execute on Torch; no native solve or launch.

The fixtures use small arrays and SimpleNamespace packets, not native model
imports. Production export/rendering is separately verified by the lead.
"""
from __future__ import annotations

import copy
import csv
import gzip
import json
import os
from pathlib import Path
import pickle
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import collect_e5f_utility_comparison as collector


TARGETS = (
    "initial_normalization", "cps_childlessness", "cps_exactly_one", "nchs_mean_age",
    "nchs_share30", "wealth_earnings", "bequest_wealth", "old_dispersion", "mean_rooms",
    "ownership_30_55", "first_birth_rooms", "family_rooms", "recent_parent_ownership",
)
COMMON = ("H0", "beta_annual", "chi", "first_birth_fixed_cost", "kappa_fert",
          "kappa_fert_continuation", "theta0")


def tables(housing="floor"):
    free = COMMON + (("h_P",) if housing == "floor" else ("delta_alpha_jump", "delta_alpha"))
    restrictions = [dict(parameter=name, lower=0., upper=10.) for name in free]
    target_rows = [dict(restriction_id=name, target=2.1 if i == 0 else 1.,
                        actual_weight=None if i == 0 else 1.) for i, name in enumerate(TARGETS)]
    objective = dict(target_rows=target_rows, parameter_restrictions=restrictions)
    fits = [dict(moment=row["restriction_id"], target=str(row["target"]), model=str(row["target"]),
                 gap="0.0", weight="" if i == 0 else "1.0", loss_contribution="" if i == 0 else "0.0")
            for i, row in enumerate(target_rows)]
    params = [dict(parameter=name, estimate="5.0", lower="0.0", upper="10.0",
                   near_bound="False", status="experimental free coordinate") for name in free]
    params += [dict(parameter=name, estimate="1.0", lower="", upper="", near_bound="", status="fixed")
               for name in collector.FIXED_NAMES]
    receipt = dict(point={name: 5. for name in free}, loss=0.)
    return objective, fits, params, receipt


def write_table(path, rows, fields):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def case_fixture(root, housing="floor"):
    arm = housing + "_linear"
    objective, fits, params, receipt = tables(housing)
    objective_path = root / "objective.json"
    objective_path.write_text(json.dumps(objective))
    contract = dict(arms={arm: dict(objective=dict(path=str(objective_path),
        sha256=collector.sha256(objective_path)), free_count=len(receipt["point"]), benefit_exponent=1.)},
        parent_source_inventory=dict(sha256="s" * 64), reference_rent=1., pension_ratio=1.)
    P = SimpleNamespace(adult_entry_clock="split_birth_vintage", utility_comparison_arm=arm,
        utility_child_benefit_exponent=1., utility_reference_rent=1., theta1=1., psi_child=1.,
        tau_pay=1., pension=1., xi_supply=[1.], tenure_choice_kappa=1., alpha_cons=1., sigma=1.,
        psi=1., phi=[1.], delta=1., tau_H=1., entrant_conversion_factor=1.)
    for row in params:
        if row["parameter"] == "income_process": row["estimate"] = "15.0"
        if row["parameter"] == "adult_entry_birth_to_household_conversion": row["estimate"] = str(1 / 2.1)
    packet = dict(parameters=P,
        solution=SimpleNamespace(p_eq=np.array([1.]), V=np.arange(3.), g=np.array([.2, .3, .5])),
        evaluation=SimpleNamespace(policy=SimpleNamespace(V=np.arange(3.), price=np.array([1.])), g_current=np.array([.2, .3, .5])),
        stationary_g_pre=np.array([.1, .2, .7]))
    case = root / "case"
    case.mkdir()
    with gzip.open(case / "initial_state.pkl.gz", "wb") as stream:
        pickle.dump(packet, stream)
    receipt.update(status="verified_experimental_point", utility_comparison_arm=arm,
        comparison_contract_sha256="c" * 64, target_system_sha256=collector.sha256(objective_path),
        source_manifest_sha256="s" * 64, free_count=len(receipt["point"]), weighted_count=12,
        selected_checkpoint_sha256="b" * 64,
        display_count=13, benefit_exponent=1., utility_reference_rent=1.,
        case_checkpoint_sha256=collector.sha256(case / "initial_state.pkl.gz"),
        annual_depreciation=1., annual_property_tax=1., price=1.,
        normalization=dict(psi_child=1., target=2.1, stationary_solve_seconds=999., status="derived_intercept"))
    (case / "receipt.json").write_text(json.dumps(receipt))
    write_table(case / "target_fit.csv", fits, collector.TARGET_FIELDS)
    write_table(case / "parameters.csv", params, collector.PARAMETER_FIELDS)
    runtime = dict(chain=SimpleNamespace(extract_moments=lambda sol, parameters: dict(tfr=2.1, all_other_moment=3.)))
    tax = SimpleNamespace(actual_parameters=lambda parameters: dict(receipt["point"]), CHECKPOINT_SHA="b" * 64)
    return contract, arm, case, runtime, tax


class NumericComparisonTests(unittest.TestCase):
    def test_frozen_array_semantics_and_no_relative_tolerance(self):
        self.assertTrue(collector.same_numerics(dict(a=np.array([np.nan, np.inf, 1.])),
                                                dict(a=np.array([np.nan, np.inf, 1.]))))
        self.assertFalse(collector.same_numerics(np.array([1e9]), np.array([1e9 + 1e-6])))
        self.assertFalse(collector.same_numerics(np.ones((1, 2)), np.ones((2, 1))))
        self.assertFalse(collector.same_numerics(dict(a=1.), dict(a=1., b=0.)))

    def test_duration_is_excluded_but_normalization_is_not(self):
        receipt = dict(normalization=dict(psi_child=.2, stationary_solves=4, stationary_solve_seconds=10.),
                       chosen_solve_seconds=2., case_checkpoint_sha256="first", loss=1.)
        other = copy.deepcopy(receipt)
        other.update(chosen_solve_seconds=20., case_checkpoint_sha256="second")
        other["normalization"]["stationary_solve_seconds"] = 100.
        self.assertTrue(collector.same_numerics(collector._receipt_science(receipt), collector._receipt_science(other)))
        other["normalization"]["stationary_solves"] += 1
        self.assertFalse(collector.same_numerics(collector._receipt_science(receipt), collector._receipt_science(other)))

    def test_original_is_compared_to_both_repeats(self):
        root = Path("/dummy")
        def record(path, *args):
            value = 1. if "original" in str(path) else 2.
            return dict(signature=dict(moment=value), pins={}, case=path)
        with patch.object(collector, "scientific_checkpoint", side_effect=record):
            with self.assertRaisesRegex(RuntimeError, "original vs repeat 0"):
                collector._verify_loaded({}, "arm", root / "original", [root / "r1", root / "r2"], 2, {}, None)

    def test_second_repeat_cannot_be_skipped(self):
        def record(path, *args):
            return dict(signature=dict(moment=2. if "r2" in str(path) else 1.), pins={}, case=path)
        with patch.object(collector, "scientific_checkpoint", side_effect=record):
            with self.assertRaisesRegex(RuntimeError, "original vs repeat 1"):
                collector._verify_loaded({}, "arm", "/dummy/original", ["/dummy/r1", "/dummy/r2"], 2, {}, None)

    def test_count_and_distinct_paths_fail_closed(self):
        for others, count in (([], 2), (["/dummy/r1"], 2), (["/dummy/r1", "/dummy/r2"], 1)):
            with self.assertRaisesRegex(RuntimeError, "require exactly"):
                collector._verify_loaded({}, "arm", "/dummy/original", others, count, {}, None)
        with self.assertRaisesRegex(RuntimeError, "distinct"):
            collector._verify_loaded({}, "arm", "/dummy/original", ["/dummy/original"], 1, {}, None)

    def test_smoke_is_not_final_repeat_claim(self):
        record = dict(signature=dict(moment=1.), pins={})
        with patch.dict(os.environ, {collector.ENV_PIN: "c" * 64}), \
                patch.object(collector, "scientific_checkpoint", return_value=record):
            smoke, _ = collector._verify_loaded({}, "arm", "/dummy/original", ["/dummy/r1"], 1, {}, None)
            final, _ = collector._verify_loaded({}, "arm", "/dummy/original", ["/dummy/r1", "/dummy/r2"], 2, {}, None)
        self.assertFalse(smoke["exact_repeat_claim"])
        self.assertTrue(final["exact_repeat_claim"])
        self.assertEqual(final["repetitions_verified"], 2)


class CompleteTableTests(unittest.TestCase):
    def test_all_floor_and_share_rows_and_pages_are_retained(self):
        for housing, count in (("floor", 8), ("shares", 9)):
            objective, fits, params, receipt = tables(housing)
            collector.validate_tables(fits, params, objective, receipt, count)
            pages = collector.parameter_pages(params)
            self.assertEqual([row for page in pages for row in page], params)
            self.assertEqual(len(params), count + 20)
            self.assertEqual(len(pages), 3)

    def test_omitted_rows_and_zero_weight_are_rejected(self):
        objective, fits, params, receipt = tables()
        with self.assertRaises(RuntimeError):
            collector.validate_tables(fits[:-1], params, objective, receipt, 8)
        with self.assertRaises(RuntimeError):
            collector.validate_tables(fits, params[:-1], objective, receipt, 8)
        fits[1]["weight"] = "0.0"
        with self.assertRaises(RuntimeError):
            collector.validate_tables(fits, params, objective, receipt, 8)

    def test_changed_weight_bound_flag_or_payload_is_rejected(self):
        objective, fits, params, receipt = tables()
        for field, value in (("lower", "1.0"), ("near_bound", "True"), ("estimate", "4.0")):
            changed = copy.deepcopy(params)
            changed[0][field] = value
            with self.assertRaises(RuntimeError):
                collector.validate_tables(fits, changed, objective, receipt, 8)
        changed = copy.deepcopy(fits)
        changed[1]["target"] = "2.0"
        with self.assertRaises(RuntimeError):
            collector.validate_tables(changed, params, objective, receipt, 8)


class SavedPacketTests(unittest.TestCase):
    def test_native_arrays_moments_and_all_tables_are_in_signature(self):
        with tempfile.TemporaryDirectory() as temp, patch.dict(os.environ, {collector.ENV_PIN: "c" * 64}):
            contract, arm, case, runtime, tax = case_fixture(Path(temp))
            item = collector.scientific_checkpoint(case, contract, arm, runtime, tax)
            self.assertEqual(item["signature"]["moments"]["all_other_moment"], 3.)
            self.assertTrue({"native_V", "native_g", "V", "g", "stationary_g_pre", "native_price",
                             "target_table", "parameter_table"}.issubset(item["signature"]))
            self.assertEqual(len(item["parameters"]), 28)

    def test_checkpoint_hash_is_verified_before_unpickle(self):
        with tempfile.TemporaryDirectory() as temp, patch.dict(os.environ, {collector.ENV_PIN: "c" * 64}):
            contract, arm, case, runtime, tax = case_fixture(Path(temp))
            (case / "initial_state.pkl.gz").write_bytes(b"corrupt")
            with patch.object(collector.pickle, "load") as load:
                with self.assertRaisesRegex(RuntimeError, "checkpoint hash"):
                    collector.scientific_checkpoint(case, contract, arm, runtime, tax)
                load.assert_not_called()

    def test_wrong_contract_or_parameter_table_fails(self):
        with tempfile.TemporaryDirectory() as temp, patch.dict(os.environ, {collector.ENV_PIN: "c" * 64}):
            contract, arm, case, runtime, tax = case_fixture(Path(temp))
            receipt = collector.read_json(case / "receipt.json")
            receipt["comparison_contract_sha256"] = "changed"
            (case / "receipt.json").write_text(json.dumps(receipt))
            with self.assertRaisesRegex(RuntimeError, "comparison_contract_sha256"):
                collector.scientific_checkpoint(case, contract, arm, runtime, tax)
            receipt["comparison_contract_sha256"] = "c" * 64
            (case / "receipt.json").write_text(json.dumps(receipt))
            params = collector.read_table(case / "parameters.csv", collector.PARAMETER_FIELDS, "parameter")
            next(row for row in params if row["parameter"] == "psi_child")["estimate"] = "2.0"
            write_table(case / "parameters.csv", params, collector.PARAMETER_FIELDS)
            with self.assertRaisesRegex(RuntimeError, "native parameter object"):
                collector.scientific_checkpoint(case, contract, arm, runtime, tax)

    def test_atomic_file_never_overwrites(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "receipt.json"
            collector.write_new_json(path, dict(value=1))
            with self.assertRaises(FileExistsError):
                collector.write_new_json(path, dict(value=2))
            self.assertEqual(collector.read_json(path), dict(value=1))

    def test_collection_failure_is_published_without_certification(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            arm_root = root / "run/floor_linear"
            arm_root.mkdir(parents=True)
            (arm_root / "controller_summary.json").write_text('{"failed":1}')
            with patch.object(collector, "_load_contract", side_effect=RuntimeError("missing second repeat")):
                result = collector.collect_arm(root / "contract.json", "floor_linear", root / "run", root / "export")
            self.assertEqual(result["status"], "comparison_incomplete")
            self.assertFalse(result["exact_repeat_claim"])
            self.assertTrue((root / "export/failure.txt").is_file())
            self.assertFalse(list((root / "export").glob("*.pdf")))
            with self.assertRaises(FileExistsError):
                collector.collect_arm(root / "contract.json", "floor_linear", root / "run", root / "export")

    def test_layout_fixture_is_separate_and_retains_all_29_rows(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            _, fits, params, _ = tables()
            # The historical report has 25 rows, before the three new utility /
            # pension descriptors. No native imports or PDF rendering in test.
            historical = [row for row in params if row["parameter"] not in
                          {"child_benefit_exponent", "utility_reference_rent", "pension_to_gross_worker_earnings"}]
            write_table(root / "targets.csv", fits, collector.TARGET_FIELDS)
            write_table(root / "parameters.csv", historical, collector.PARAMETER_FIELDS)
            graphs = root / "graphs"
            graphs.mkdir()
            for name in collector.STANDARD_NAMES:
                (graphs / (name + ".png")).write_bytes(b"placeholder fixture image")
            with patch.object(collector, "render_report", return_value=dict(pdf_name="layout_fixture.pdf")) as render:
                result = collector.render_layout_fixture(root / "targets.csv", root / "parameters.csv", graphs, root / "layout")
            self.assertTrue(render.call_args.kwargs["layout_fixture"])
            selected = render.call_args.args[2]
            self.assertEqual(len(selected["fits"]), 13)
            self.assertEqual(len(selected["parameters"]), 29)
            self.assertFalse(result["scientific_repeat_claim"])
            self.assertTrue(result["layout_fixture"])


if __name__ == "__main__":
    unittest.main()
