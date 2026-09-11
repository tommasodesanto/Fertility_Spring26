#!/usr/bin/env python3
"""Read two pinned checkpoints; export existing graph arrays, never solve a model.

Run in the c6dd3508 snapshot with one CPU/16 GB and an external five-minute
wall cap. Output must be a new directory. Inputs and original graphs are read only.
"""
from __future__ import annotations

import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_name] = "1"
import argparse
import csv
import gc
import gzip
import hashlib
import json
from pathlib import Path
import pickle
import sys
import time

SNAPSHOT = Path("/scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508")
RELATIVE = Path("output/utility_fiscal_decomposition_20260911a")
CONTRACT_SHA = "65c7b48a7d26a8c0ca67632acf499aeb53705d2b7436469833e61b6df8ece6a9"
CHECKPOINTS = {
    "old_old": "60ccfff3d2a16a007f093c997c4a9de0f1ef81b4ae16d286ec6003edb5539cf1",
    "new_balanced": "abe1d627e3df9b5af422514b36b9d4c4dbb932bb834f31d68c7d201c90050152",
}


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verify(path, expected):
    actual = digest(path)
    if actual != expected:
        raise RuntimeError(f"Hash mismatch: {path}: {actual}")


def write_csv(path, rows):
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    contract_path = SNAPSHOT / RELATIVE / "contract.json"
    verify(contract_path, CONTRACT_SHA)
    contract = json.loads(contract_path.read_text())
    # Validate named immutable files, without a directory scan or model call.
    for name, expected in contract["source_sha256"].items():
        source = (SNAPSHOT / name).resolve()
        if not source.is_relative_to(SNAPSHOT):
            raise RuntimeError("Source path escapes pinned snapshot")
        verify(source, expected)
    for case, expected in CHECKPOINTS.items():
        verify(SNAPSHOT / RELATIVE / case / "repetition_01/initial_state.pkl.gz", expected)
    args.output.mkdir(parents=True, exist_ok=False)
    sys.path[:0] = [str(SNAPSHOT / "code/model"), str(SNAPSHOT / "code/model/tools")]
    import numpy as np
    import run_e5f_open_population_transition as transition
    # This only loads the sequential runtime needed by checkpoint classes.
    # No parameter initialization, Bellman, equilibrium, or population operator.
    transition.configure_sequential_model()
    from intergen_eqscale_seq_optimized.parameters import readiness_settled_state
    receipts = []
    for case, expected in CHECKPOINTS.items():
        path = SNAPSHOT / RELATIVE / case / "repetition_01/initial_state.pkl.gz"
        with gzip.open(path, "rb") as stream:
            packet = pickle.load(stream)
        P, e = packet["parameters"], packet["evaluation"]
        p, bg = e.policy, np.asarray(packet["b_grid"]).reshape(-1)
        if bool(getattr(P, "joint_nested_choice", False)) or int(P.I) != 1:
            raise RuntimeError("This reduction is pinned to sequential one-market graphs")
        z = np.asarray(P.z_grid).reshape(-1)
        groups = np.asarray(P.permanent_income_group_index, dtype=int)
        base_states = np.asarray(P.permanent_income_base_state_index, dtype=int)
        levels = np.asarray(P.permanent_income_level_values)
        if len(z) != 15 or len(bg) != p.V.shape[0]:
            raise RuntimeError("Unexpected income/wealth grid")
        cs = readiness_settled_state(P)
        pools = {"pre_fertility": e.g_pre, "post_fertility_pre_tenure": e.g_post_fertility,
                 "current_post_tenure": e.g_current}
        rows, screens, age_rows = [], [], []
        for age in (30, 42):
            j = int(round((age - float(P.age_start)) / float(P.da)))
            if float(P.age_start + j * P.da) != age:
                raise RuntimeError("Requested policy age is not an exact grid age")
            for zz in range(len(z)):
                idx = (slice(None), 0, 0, j, zz, 0, cs)
                values = np.asarray(p.V[idx])
                valid = values > -1e9  # exact plotting mask, not an occupancy test
                choice = np.asarray(p.tenure_choice[idx], dtype=int)
                tp = np.asarray(p.tenure_probs[idx])
                owner = np.sum(tp[:, 1:], axis=1)
                hr = np.asarray(p.hR_pol[idx])
                h_graph = np.where(choice <= 0, hr,
                    np.asarray(P.H_own)[np.maximum(choice - 1, 0)])
                fertility = np.asarray(p.fert_probs[:, 0, 0, j, zz, :]) @ np.arange(P.n_parity)
                mass = {name: np.asarray(g[idx]) for name, g in pools.items()}
                adjacent = valid[:-1] & valid[1:]
                drops = adjacent & (np.diff(owner) < -1e-7)
                for ib, wealth in enumerate(bg):
                    row = dict(case=case, age=age, z_index=zz, z=float(z[zz]),
                        permanent_index=int(groups[zz]), permanent_level=float(levels[groups[zz]]),
                        changing_income_state=int(base_states[zz]), wealth_index=ib,
                        wealth=float(wealth), renter=0, location=0, parity=0, child_state=int(cs),
                        graph_value_mask=bool(valid[ib]), value=float(values[ib]),
                        consumption_conditional_renter=float(p.c_pol[idx][ib]),
                        housing_conditional_renter=float(hr[ib]),
                        selected_tenure=int(choice[ib]), housing_graph_selected_tenure=float(h_graph[ib]),
                        owner_entry_probability=float(owner[ib]),
                        graph_fertility_scalar=float(fertility[ib]),
                        tenure_probability_sum=float(tp[ib].sum()))
                    for name, gm in mass.items():
                        row[name + "_mass"] = float(gm[ib])
                    for destination in range(tp.shape[1]):
                        row[f"tenure_probability_{destination}"] = float(tp[ib, destination])
                    rows.append(row)
                for timing, gm in mass.items():
                    occupied_lower = gm[:-1] > 1e-12
                    screens.append(dict(case=case, age=age, z_index=zz, z=float(z[zz]),
                        timing=timing, state_mass=float(gm.sum()), valid_nodes=int(valid.sum()),
                        negative_owner_steps_all_valid=int(drops.sum()),
                        negative_owner_steps_occupied_lower=int((drops & occupied_lower).sum()),
                        lower_node_mass_at_negative_owner_steps=float(gm[:-1][drops].sum()),
                        largest_owner_drop_all_valid=float(max(0., -np.min(np.diff(owner)[adjacent])))
                            if np.any(adjacent) else 0.,
                        largest_owner_drop_occupied_lower=float(max(0., -np.min(np.diff(owner)[adjacent & occupied_lower])))
                            if np.any(adjacent & occupied_lower) else 0.,
                        masked_node_mass=float(gm[~valid].sum()),
                        occupied_threshold=1e-12, probability_step_threshold=1e-7))
        g = np.asarray(e.g_current)
        for j in range(int(P.J)):
            for scope, ids in [("combined", [zz]) for zz in range(len(z))] + [
                    ("permanent", np.flatnonzero(groups == k).tolist()) for k in range(len(levels))]:
                denominator = sum(float(g[:, :, :, j, zz, :, :].sum()) for zz in ids)
                numerator = sum(float(g[:, 1:, :, j, zz, :, :].sum()) for zz in ids)
                age_rows.append(dict(case=case, age=float(P.age_start+j*P.da), scope=scope,
                    group_index=ids[0] if scope == "combined" else int(groups[ids[0]]),
                    permanent_level=float(levels[groups[ids[0]]]),
                    z=float(z[ids[0]]) if scope == "combined" else "",
                    mass=denominator, owner_mass=numerator,
                    owner_rate=numerator/denominator if denominator > 1e-14 else ""))
        write_csv(args.output / f"{case}_conditional_policy.csv", rows)
        write_csv(args.output / f"{case}_owner_probability_steps.csv", screens)
        write_csv(args.output / f"{case}_ownership_age_type.csv", age_rows)
        receipts.append(dict(case=case, input=str(path), checkpoint_sha256=expected,
            wealth_points=len(bg), income_states=len(z), readiness_state=int(cs),
            price=np.asarray(p.price).tolist(), parameter_phi=np.asarray(P.phi).tolist(),
            owner_products=np.asarray(P.H_own).tolist(), rows=len(rows),
            distribution_timing="g_pre: before birth; g_post_fertility: before tenure; g_current: after tenure/location",
            policy_condition="renter, location 0, parity 0, readiness-settled; tenure choice is conditional on the post-fertility state",
            housing_graph_definition="selected tenure product or conditional rental housing; not probability-weighted housing",
            elapsed_seconds=time.monotonic()-started))
        del packet, e, p, g, pools, rows, screens, age_rows
        gc.collect()
    receipt = dict(status="completed_read_only_checkpoint_reduction", model_solves=0,
        contract_sha256=CONTRACT_SHA, source_files_verified=len(contract["source_sha256"]),
        cases=receipts, elapsed_seconds=time.monotonic()-started,
        output_sha256={q.name:digest(q) for q in args.output.glob("*.csv")})
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status":receipt["status"], "seconds":receipt["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
