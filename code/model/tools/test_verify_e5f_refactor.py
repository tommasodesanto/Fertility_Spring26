from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from verify_e5f_refactor import ComparisonError, compare, main


SOURCE_IDS = ("source-reference", "source-candidate")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_bundle(root: Path, source: str = SOURCE_IDS[0], input_fp: str = "input-a") -> Path:
    root.mkdir()
    np.save(root / "values.npy", np.array([1.0, np.nan, np.inf], dtype=np.float64))
    np.savez(
        root / "decisions.npz",
        branch=np.array([0, 1], dtype=np.int8),
        gate=np.array([True, False], dtype=bool),
    )
    (root / "payload.json").write_text(
        json.dumps(
            {
                "gate": True,
                "timestamp": "2026-09-16T12:00:00Z",
                "duration_seconds": 1.25,
                "run_metadata": {"timestamp": "scientific-value", "duration_seconds": 2.5},
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    (root / "table.csv").write_text("state,value\n0,1.0\n1,2.0\n", encoding="utf-8")
    manifest = {
        "schema": "e5f_refactor_evidence_v1",
        "provenance": {"source_identity": source, "input_fingerprint": input_fp},
        "branch": "branch-a",
        "numerical_gate": "converged",
        "run_metadata": {
            "timestamp": "run-specific",
            "elapsed_seconds": 1.0,
            "worker": "fixed-worker",
        },
        "artifacts": {},
    }
    for name, units in (
        ("values.npy", "index"),
        ("decisions.npz", "boolean"),
        ("payload.json", "metadata"),
        ("table.csv", "table"),
    ):
        manifest["artifacts"][name] = {"path": name, "sha256": _sha256(root / name), "units": units}
    (root / "manifest.json").write_text(json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8")
    return root


def read_manifest(root: Path) -> dict:
    return json.loads((root / "manifest.json").read_text(encoding="utf-8"))


def write_manifest(root: Path, manifest: dict) -> None:
    (root / "manifest.json").write_text(json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8")


class VerifyTests(unittest.TestCase):
    def assert_bundle_failure(self, mutate, message: str | None = None) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = make_bundle(root / "reference", SOURCE_IDS[0])
            candidate = make_bundle(root / "candidate", SOURCE_IDS[1])
            mutate(reference, candidate)
            with self.subTest(message=message):
                result = compare(reference, candidate, expected_source_ids=SOURCE_IDS)
                self.assertEqual(result["status"], "fail")
                if message is not None:
                    self.assertIn(message, "\n".join(result["failures"]))

    def assert_bundle_error(self, mutate, message: str | None = None) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = make_bundle(root / "reference", SOURCE_IDS[0])
            candidate = make_bundle(root / "candidate", SOURCE_IDS[1])
            mutate(reference, candidate)
            with self.subTest(message=message):
                with self.assertRaises(ComparisonError) as caught:
                    compare(reference, candidate, expected_source_ids=SOURCE_IDS)
                if message is not None:
                    self.assertIn(message, str(caught.exception))

    def test_positive_bundle_allows_distinct_sources_and_manifest_timings(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = make_bundle(root / "reference", SOURCE_IDS[0])
            candidate = make_bundle(root / "candidate", SOURCE_IDS[1])
            manifest = read_manifest(candidate)
            manifest["run_metadata"]["timestamp"] = "different-run"
            manifest["run_metadata"]["elapsed_seconds"] = 999.0
            write_manifest(candidate, manifest)
            result = compare(reference, candidate, expected_source_ids=SOURCE_IDS)
            self.assertEqual(result["status"], "pass")
            self.assertEqual(result["scope"], "artifact_bundle")

    def test_wrong_source_and_input_fingerprint_are_rejected(self) -> None:
        self.assert_bundle_error(
            lambda _reference, candidate: write_manifest(
                candidate,
                {**read_manifest(candidate), "provenance": {"source_identity": "wrong", "input_fingerprint": "input-a"}},
            ),
            "source identity",
        )
        self.assert_bundle_error(
            lambda _reference, candidate: write_manifest(
                candidate,
                {**read_manifest(candidate), "provenance": {"source_identity": SOURCE_IDS[1], "input_fingerprint": "other"}},
            ),
            "input fingerprints differ",
        )

    def test_array_values_shape_dtype_and_tiny_numeric_changes_are_exact(self) -> None:
        cases = (
            (lambda _r, c: np.save(c / "values.npy", np.array([1.00000000000001, np.nan, np.inf])), "array values differ"),
            (lambda _r, c: np.save(c / "values.npy", np.array([[1.0, np.nan, np.inf]])), "shape/dtype differ"),
            (lambda _r, c: np.save(c / "values.npy", np.array([1, 0, 2], dtype=np.int64)), "shape/dtype differ"),
        )
        def with_fresh_hash(mutate):
            def wrapped(_reference, candidate):
                mutate(_reference, candidate)
                manifest = read_manifest(candidate)
                manifest["artifacts"]["values.npy"]["sha256"] = _sha256(candidate / "values.npy")
                write_manifest(candidate, manifest)
            return wrapped
        for mutate, message in cases:
            self.assert_bundle_failure(with_fresh_hash(mutate), message)

    def test_discrete_branch_bool_integer_and_numerical_gate_changes(self) -> None:
        def mutate_branch(_reference, candidate):
            np.savez(candidate / "decisions.npz", branch=np.array([1, 1], dtype=np.int8), gate=np.array([True, False]))
            manifest = read_manifest(candidate)
            manifest["artifacts"]["decisions.npz"]["sha256"] = _sha256(candidate / "decisions.npz")
            write_manifest(candidate, manifest)

        def mutate_bool(_reference, candidate):
            (candidate / "payload.json").write_text(
                json.dumps({"gate": 1, "timestamp": "2026-09-16T12:00:00Z", "duration_seconds": 1.25, "run_metadata": {"timestamp": "scientific-value", "duration_seconds": 2.5}}),
                encoding="utf-8",
            )
            manifest = read_manifest(candidate)
            manifest["artifacts"]["payload.json"]["sha256"] = _sha256(candidate / "payload.json")
            write_manifest(candidate, manifest)

        def mutate_gate(_reference, candidate):
            manifest = read_manifest(candidate)
            manifest["numerical_gate"] = "failed"
            write_manifest(candidate, manifest)

        for mutate, message in ((mutate_branch, "array values differ"), (mutate_bool, "JSON values differ"), (mutate_gate, "manifest scientific metadata differ")):
            self.assert_bundle_failure(mutate, message)

    def test_json_timestamp_duration_and_run_metadata_are_exact(self) -> None:
        for field, value in (("timestamp", "different"), ("duration_seconds", 1.2500000000001), ("run_metadata", {"timestamp": "changed", "duration_seconds": 2.5})):
            def mutate(_reference, candidate, field=field, value=value):
                payload = json.loads((candidate / "payload.json").read_text(encoding="utf-8"))
                payload[field] = value
                (candidate / "payload.json").write_text(json.dumps(payload), encoding="utf-8")
                manifest = read_manifest(candidate)
                manifest["artifacts"]["payload.json"]["sha256"] = _sha256(candidate / "payload.json")
                write_manifest(candidate, manifest)

            self.assert_bundle_failure(mutate, "JSON values differ")

    def test_manifest_unknown_metadata_and_units_are_exact(self) -> None:
        self.assert_bundle_failure(
            lambda _r, c: write_manifest(c, {**read_manifest(c), "new_scientific_key": 1}),
            "manifest scientific metadata differ",
        )
        def units(_r, c):
            manifest = read_manifest(c)
            manifest["artifacts"]["values.npy"]["units"] = "wrong-units"
            write_manifest(c, manifest)
        self.assert_bundle_failure(units, "manifest scientific metadata differ")

        def unknown_run_metadata(_r, c):
            manifest = read_manifest(c)
            manifest["run_metadata"]["worker"] = "different-worker"
            write_manifest(c, manifest)
        self.assert_bundle_failure(unknown_run_metadata, "manifest scientific metadata differ")

        def reserved_name(_r, c):
            manifest = read_manifest(c)
            manifest["artifacts"]["manifest.json"] = {
                "path": "payload.json",
                "sha256": _sha256(c / "payload.json"),
                "units": "reserved",
            }
            write_manifest(c, manifest)
        self.assert_bundle_error(reserved_name, "logical name manifest.json is reserved")

    def test_missing_artifact_key_file_stale_hash_and_path_escape_fail_closed(self) -> None:
        def missing_key(_r, c):
            manifest = read_manifest(c)
            del manifest["artifacts"]["table.csv"]
            write_manifest(c, manifest)
        self.assert_bundle_error(missing_key, "artifact keys differ")

        def missing_file(_r, c):
            (c / "table.csv").unlink()
        self.assert_bundle_error(missing_file, "artifact missing")

        def stale_hash(_r, c):
            (c / "table.csv").write_text("state,value\n0,99.0\n1,2.0\n", encoding="utf-8")
        self.assert_bundle_error(stale_hash, "sha256 claim")

        def escape(_r, c):
            manifest = read_manifest(c)
            manifest["artifacts"]["table.csv"]["path"] = "../table.csv"
            write_manifest(c, manifest)
        self.assert_bundle_error(escape, "escapes bundle")

    def test_self_check_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            bundle_root = make_bundle(Path(directory) / "bundle")
            with self.assertRaisesRegex(ComparisonError, "self-check"):
                compare(bundle_root, bundle_root, expected_source_ids=("bundle", "bundle"))

    def test_object_npy_npz_and_corrupt_npz_are_rejected(self) -> None:
        def object_npy(_r, c):
            np.save(c / "values.npy", np.array([{"x": 1}], dtype=object), allow_pickle=True)
            manifest = read_manifest(c)
            manifest["artifacts"]["values.npy"]["sha256"] = _sha256(c / "values.npy")
            write_manifest(c, manifest)
        self.assert_bundle_error(object_npy, "pickle/object")

        def object_npz(_r, c):
            np.savez(c / "decisions.npz", branch=np.array([{"x": 1}], dtype=object), gate=np.array([True, False]))
            manifest = read_manifest(c)
            manifest["artifacts"]["decisions.npz"]["sha256"] = _sha256(c / "decisions.npz")
            write_manifest(c, manifest)
        self.assert_bundle_error(object_npz, "object/pickle")

        def corrupt_npz(_r, c):
            (c / "decisions.npz").write_bytes(b"not an npz archive")
            manifest = read_manifest(c)
            manifest["artifacts"]["decisions.npz"]["sha256"] = _sha256(c / "decisions.npz")
            write_manifest(c, manifest)
        self.assert_bundle_error(corrupt_npz, "NPY/NPZ load failed")

    def test_valid_csv_and_isolated_changed_cell(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.csv"
            candidate = root / "candidate.csv"
            reference.write_text("a,b\n1,2\n", encoding="utf-8")
            candidate.write_text("a,b\n1,2\n", encoding="utf-8")
            self.assertEqual(compare(reference, candidate)["status"], "pass")
            candidate.write_text("a,b\n1,3\n", encoding="utf-8")
            self.assertEqual(compare(reference, candidate)["failures"], ["reference.csv: CSV cells differ"])

    def test_artifact_only_manifest_json_is_compared_exactly(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference" / "manifest.json"
            candidate = root / "candidate" / "manifest.json"
            reference.parent.mkdir()
            candidate.parent.mkdir()
            reference.write_text('{"timestamp":"a","value":1}\n', encoding="utf-8")
            candidate.write_text('{"timestamp":"b","value":1}\n', encoding="utf-8")
            result = compare(reference, candidate)
            self.assertEqual(result["status"], "fail")
            self.assertEqual(result["failures"], ["manifest.json: JSON values differ"])

    def test_main_writes_structured_errors_for_malformed_json_manifest_and_paths(self) -> None:
        cases = ("malformed", "missing_path", "wrongtype_path", "invalid_utf8", "invalid_artifact_utf8", "corrupt_npz")
        for case in cases:
            with self.subTest(case=case), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                reference = make_bundle(root / "reference")
                candidate = make_bundle(root / "candidate", SOURCE_IDS[1])
                manifest = read_manifest(candidate)
                if case == "malformed":
                    (candidate / "manifest.json").write_text("{bad", encoding="utf-8")
                elif case == "missing_path":
                    del manifest["artifacts"]["table.csv"]["path"]
                    write_manifest(candidate, manifest)
                elif case == "wrongtype_path":
                    manifest["artifacts"]["table.csv"]["path"] = 17
                    write_manifest(candidate, manifest)
                elif case == "invalid_artifact_utf8":
                    (candidate / "payload.json").write_bytes(b"{\xff")
                    manifest["artifacts"]["payload.json"]["sha256"] = _sha256(candidate / "payload.json")
                    write_manifest(candidate, manifest)
                elif case == "corrupt_npz":
                    (candidate / "decisions.npz").write_bytes(b"not an npz archive")
                    manifest["artifacts"]["decisions.npz"]["sha256"] = _sha256(candidate / "decisions.npz")
                    write_manifest(candidate, manifest)
                else:
                    (candidate / "manifest.json").write_bytes(b"{\xff")
                output = root / "receipt"
                self.assertEqual(
                    main([
                        str(reference),
                        str(candidate),
                        "--output-dir",
                        str(output),
                        "--reference-source-id",
                        SOURCE_IDS[0],
                        "--candidate-source-id",
                        SOURCE_IDS[1],
                    ]),
                    1,
                )
                receipt = json.loads((output / "verification.json").read_text(encoding="utf-8"))
                self.assertEqual(receipt["status"], "error")
                self.assertTrue(receipt["failures"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
