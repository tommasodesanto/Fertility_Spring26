"""Isolated contract test for fixed fertility-scale probes."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path


MODEL_ROOT = Path(__file__).resolve().parents[2]


def run_contract(environment: dict[str, str]) -> dict[str, object]:
    environment["PYTHONPATH"] = os.pathsep.join(
        filter(
            None,
            (str(MODEL_ROOT), environment.get("PYTHONPATH", "")),
        )
    )
    command = (
        "import json; "
        "import intergen_eqscale_seq_optimized.run_e1_chain as c; "
        "print(json.dumps({'domain': [r[0] for r in c.DOMAIN], "
        "'fixed': c.FIXED, "
        "'seed': c.apply_probe_seed_overrides({'psi_child': 0.0})}))"
    )
    completed = subprocess.run(
        [sys.executable, "-c", command],
        check=True,
        capture_output=True,
        text=True,
        env=environment,
    )
    return json.loads(completed.stdout)


def test_e5f_probe_can_fix_both_fertility_scales() -> None:
    environment = dict(os.environ)
    environment.update(
        E5="1",
        E5_MATURATION_REPAIR="1",
        E5F="1",
        E5_PROBE_FIX_KE="0.15",
        E5_PROBE_FIX_KC="0.08",
        E5_PROBE_SEED_PSI="0.10",
    )
    payload = run_contract(environment)
    assert "kappa_fert" not in payload["domain"]
    assert "kappa_fert_continuation" not in payload["domain"]
    assert payload["fixed"]["kappa_fert"] == 0.15
    assert payload["fixed"]["kappa_fert_continuation"] == 0.08
    assert payload["seed"]["psi_child"] == 0.10


def test_e5f_probe_is_default_off() -> None:
    environment = dict(os.environ)
    environment.update(E5="1", E5_MATURATION_REPAIR="1", E5F="1")
    for name in ("E5_PROBE_FIX_KE", "E5_PROBE_FIX_KC", "E5_PROBE_SEED_PSI"):
        environment.pop(name, None)
    payload = run_contract(environment)
    assert "kappa_fert" in payload["domain"]
    assert "kappa_fert_continuation" in payload["domain"]
    assert "kappa_fert" not in payload["fixed"]
    assert "kappa_fert_continuation" not in payload["fixed"]
    assert payload["seed"]["psi_child"] == 0.0
