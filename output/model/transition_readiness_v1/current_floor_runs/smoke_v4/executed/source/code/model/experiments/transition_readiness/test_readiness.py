"""Zero-solve fail-closed and deployed-source identity checks."""
import json
from pathlib import Path
import unittest

import legacy_changed_psi as diagnostic
import selected_adapter as adapter


class ReadinessTests(unittest.TestCase):
    def test_exact_three_mapping_loop_cap(self):
        for count in range(3):
            diagnostic.reserve_mapping(count)
        with self.assertRaisesRegex(adapter.ReadinessBlocked, 'Three-map'):
            diagnostic.reserve_mapping(3)

    def test_actual_frozen_reference_and_deployed_source_pins(self):
        proof = diagnostic.preflight()
        self.assertEqual(proof['model_calls'], 0)
        prior = json.loads((diagnostic.ROOT / 'output/model/fixed_reference_transition_20260928/four_shock_v1/budget_diagnostic_v2/configs/smoke.json').read_text())
        ours = json.loads((diagnostic.HERE / 'legacy_source_pins.json').read_text())
        self.assertEqual(ours, prior['source_pins'])

    def test_missing_selected_artifact_never_calls_native_setup(self):
        called = []
        with self.assertRaises(adapter.ReadinessBlocked):
            adapter.setup_authenticated({'path': '/nonexistent-selected-calibration.json', 'sha256': '0'*64},
                engine_contract={}, setup=lambda state: called.append(state))
        self.assertEqual(called, [])

    def test_existing_source_with_wrong_hash_rejected(self):
        with self.assertRaisesRegex(adapter.ReadinessBlocked, 'hash mismatch'):
            adapter.pinned({'path': str(diagnostic.HERE / 'legacy_source_pins.json'), 'sha256': '0'*64})

    def test_relative_artifacts_rejected(self):
        with self.assertRaises(adapter.ReadinessBlocked):
            adapter.pinned({'path': 'legacy_source_pins.json', 'sha256': '0'*64})

    def test_legacy_manifest_cannot_be_selected_adapter_manifest(self):
        manifest = diagnostic.ROOT / 'output/model/fertility_identification_20260928/fixed_reference_manifest.json'
        with self.assertRaisesRegex(adapter.ReadinessBlocked, 'schema'):
            adapter.authenticate({'path': str(manifest), 'sha256': diagnostic.REFERENCE_SHA})


if __name__ == '__main__':
    unittest.main()
