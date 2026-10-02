"""Read first-birth flow from a saved native stationary one-step packet; no solve.

Matches run_e5f_transition_calibration.first_birth_accounting_by_age, including
fecundity, the one-child first-birth choice, and the all-age childless risk set.
Run separately for each isolated purchase-rule engine in a fresh interpreter.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pickle
import sys
from pathlib import Path

import numpy as np


PARAMETERS_SHA256 = "f8871ea7f026b0900f2cd89e2ecbe7f24fbc8b08e83ec46227533a11dea2e163"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine-root", type=Path, required=True)
    parser.add_argument("--frozen-model-root", type=Path, required=True)
    parser.add_argument("--packet", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    parameters_source = args.engine_root / "small_credit_lab/engine/parameters.py"
    if sha(parameters_source) != PARAMETERS_SHA256:
        raise ValueError("Isolated fertility parameter source changed")
    sys.path.insert(0, str(args.engine_root))
    from small_credit_lab.engine.parameters import get_fecundity_by_age, readiness_settled_state

    sys.path.extend([str(args.frozen_model_root), str(args.frozen_model_root / "tools")])
    with gzip.open(args.packet, "rb") as stream:
        packet = pickle.load(stream)
    evaluation, P = packet["evaluation"], packet["parameters"]
    fecundity = get_fecundity_by_age(P)
    settled = readiness_settled_state(P)
    flows, risk = [], []
    for age in range(int(P.J)):
        childless = evaluation.g_pre[:, :, :, age, :, 0, settled]
        risk.append(float(childless.sum()))
        flow = 0.0
        if int(P.A_f_start) <= age + 1 <= int(P.A_f_end):
            for income in range(evaluation.g_pre.shape[4]):
                attempt = evaluation.policy.fert_probs[:, :, :, age, income, 1]
                flow += float(np.sum(float(fecundity[age]) * childless[:, :, :, income] * attempt))
        flows.append(flow)
    saved = json.loads(args.mapping.read_text())["fertility"][0]
    saved_flow = np.fromstring(saved["birth_flow_first"].strip("[]"), sep=" ")
    saved_risk = np.fromstring(saved["childless_at_risk_mass"].strip("[]"), sep=" ")
    if saved_flow.shape != (int(P.J),) or saved_risk.shape != (int(P.J),):
        raise ValueError("Saved observer age vectors have the wrong shape")
    max_printed_difference = max(float(np.max(np.abs(np.asarray(flows) - saved_flow))),
                                 float(np.max(np.abs(np.asarray(risk) - saved_risk))))
    if max_printed_difference > 1e-8:
        raise ValueError("Reconstructed first-birth accounting disagrees with saved observer")
    total_flow, total_risk = float(sum(flows)), float(sum(risk))
    result = dict(
        first_birth_flow=total_flow, childless_at_risk_mass_all_ages=total_risk,
        first_birth_hazard_all_childless=total_flow / total_risk,
        maximum_absolute_difference_from_saved_printed_age_vectors=max_printed_difference,
        packet_sha256=sha(args.packet), mapping_sha256=sha(args.mapping),
        isolated_parameters_sha256=PARAMETERS_SHA256,
        source_observer="run_e5f_transition_calibration.py:first_birth_accounting_by_age",
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x") as stream:
        stream.write(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
