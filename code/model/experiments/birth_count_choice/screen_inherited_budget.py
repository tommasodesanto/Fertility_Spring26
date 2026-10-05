#!/usr/bin/env python3
"""Inert CURRENT-BUDGET stay-or-exit witness screen (not a feasibility classifier).

This module imports no model code.  It reads only a saved packet and evaluates
the one-period, same-tenure resource inequality documented in the frozen
``household.py`` helpers ``native_due_owner_floor`` and
``native_due_death_floor`` and ``kernels.py::full_owner_block_kernel``.
Passing means only that a current-budget *stay witness* exists at that cell.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pickle
from pathlib import Path
from typing import Any

import numpy as np

PACKET_SHA256 = "a2c6b2b266ef524fa8d020da1623b5d38e55caa4acfd364a07cf2260762244fe"
CLOSURE_SHA256 = "61861173b391aefafa36bae97040951f7f100cd99e4ed989f5af14785d12d6ad"
EXPECTED_SCALE = 1.0009489339264241
TOL = 1e-10


class _InertRecord:
    """Pickle placeholder accepting arbitrary constructor/state conventions."""
    def __new__(cls, *args: Any, **kwargs: Any) -> "_InertRecord":
        return object.__new__(cls)

    def __setstate__(self, state: Any) -> None:
        if isinstance(state, dict):
            self.__dict__.update(state)
        elif isinstance(state, tuple) and len(state) == 2 and isinstance(state[1], dict):
            self.__dict__.update(state[1])
        else:
            self.__dict__["_pickle_state"] = state


class _InertUnpickler(pickle.Unpickler):
    def find_class(self, module: str, name: str) -> Any:
        if module == "builtins" or module == "types" or module.startswith("numpy"):
            return super().find_class(module, name)
        return _InertRecord


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def array_sha256(a: Any) -> str:
    x = np.asarray(a)
    return hashlib.sha256(x.tobytes(order="C")).hexdigest()


def load_packet(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rb") as f:
        packet = _InertUnpickler(f).load()
    if not isinstance(packet, dict):
        raise ValueError("packet must be a dictionary")
    need = {"parameters", "b_grid", "shared", "stationary_g_pre"}
    if need - packet.keys():
        raise ValueError(f"packet missing required keys: {sorted(need - packet.keys())}")
    return packet


def _bool(P: Any, name: str) -> bool:
    return bool(getattr(P, name, False))


def _attr(P: Any, name: str) -> Any:
    if not hasattr(P, name):
        raise NotImplementedError(f"packet lacks required current-budget field {name}")
    return getattr(P, name)


def _matrix(shared: Any, base: str, nchild: int, nstate: int) -> np.ndarray:
    direct = getattr(shared, base, None)
    if direct is not None:
        x = np.asarray(direct, dtype=float)
        if x.shape == (nchild, nstate):
            return x
    flat = getattr(shared, base + "_flat", None)
    if flat is None:
        raise NotImplementedError(f"shared lacks {base} matrix/F-order flat array")
    return np.asarray(flat, dtype=float).reshape((nchild, nstate), order="F")


def _death_possible(P: Any, j: int) -> bool:
    if j == int(_attr(P, "J")) - 1:
        return True
    if _bool(P, "use_age_survival"):
        s = np.asarray(_attr(P, "survival_probs"), dtype=float)
        return j >= s.size or float(s[j]) < 1.0
    return False


def income(P: Any, i: int, j: int, z: float) -> float:
    """Frozen shared.py income logic; P.income is already after tax."""
    y = float(np.asarray(_attr(P, "income"), dtype=float)[i, j])
    if j < int(_attr(P, "J_R")):
        return y * z
    return y * (1.0 + float(_attr(P, "retirement_income_z_scale")) * (z - 1.0))


def check_supported(P: Any, shared: Any) -> dict[str, Any]:
    """Require exactly the reviewed native accounting branch; never approximate."""
    if not _bool(P, "native_purchase_income") or not _bool(P, "native_due_stayer_credit"):
        raise NotImplementedError("screen requires native_purchase_income and native_due_stayer_credit")
    forbidden = ("native_solvency_credit", "natural_credit", "experimental_natural_solvency", "use_pti_constraint", "estate_receiver")
    active = [x for x in forbidden if _bool(P, x)]
    if active:
        raise NotImplementedError("unsupported active budget flag(s): " + ", ".join(active))
    active = [x for x in ("child_earnings_penalty", "rental_wedge") if _bool(P, x)]
    if active:
        raise NotImplementedError("unsupported active budget flag(s): " + ", ".join(active))
    if float(getattr(P, "unsecured_credit_limit", 0.0) or 0.0) != 0.0:
        raise NotImplementedError("screen requires zero unsecured credit")
    if np.any(np.asarray(getattr(shared, "gb_flat", 0.0), dtype=float) != 0.0):
        raise NotImplementedError("screen requires zero transfers")
    if float(getattr(P, "property_tax_lump_sum_transfer", 0.0) or 0.0) != 0.0:
        raise NotImplementedError("screen does not support a fiscal transfer")
    ltv = np.asarray(_attr(P, "owner_ltv_multipliers"), dtype=float)
    if not np.all(np.isfinite(ltv)) or not np.all(ltv == 1.0):
        raise NotImplementedError("screen requires all owner_ltv_multipliers to equal one")
    if float(_attr(P, "owner_size_cost")) != 0.0:
        raise NotImplementedError("screen requires owner_size_cost=0")
    return {x: _bool(P, x) for x in ("native_purchase_income", "native_due_stayer_credit", *forbidden,
                                      "child_earnings_penalty", "rental_wedge")}


def screen(packet: dict[str, Any], *, q: float, r: float, pension: float, scale: float) -> dict[str, Any]:
    """Evaluate same-tenure, no-transaction current-budget witnesses only."""
    P, S = packet["parameters"], packet["shared"]
    flags = check_supported(P, S)
    if not all(np.isfinite(x) and x > 0.0 for x in (q, r, pension, scale)):
        raise ValueError("q, r, pension, and scale must be finite and positive")
    serialized_pension = float(_attr(P, "pension"))
    if pension != serialized_pension:
        raise ValueError(f"requested pension {pension!r} differs from serialized P.pension {serialized_pension!r}")
    bg, g = np.asarray(packet["b_grid"], float), np.asarray(packet["stationary_g_pre"], float)
    if bg.ndim != 1 or bg.size == 0 or not np.all(np.isfinite(bg)) or not np.all(np.diff(bg) >= 0) or g.ndim != 7:
        raise ValueError("invalid b_grid or stationary_g_pre shape")
    if not np.all(np.isfinite(g)) or np.any(g < 0): raise ValueError("population must be finite and nonnegative")
    nb, nt, I, J, nz, nc, ns = g.shape
    if (nb, nt, I, J, nz) != (bg.size, int(_attr(P, "n_house")) + 1, int(_attr(P, "I")), int(_attr(P, "J")), int(_attr(P, "Nz"))):
        raise ValueError("population shape disagrees with serialized parameters")
    cb, hb = _matrix(S, "cb", nc, ns), _matrix(S, "hb", nc, ns)
    phi = np.asarray(_attr(S, "phi_choice"), float)
    if phi.shape != (I, nt, nc, ns): raise ValueError("unexpected phi_choice shape")
    H, z = np.asarray(_attr(P, "H_own"), float), np.asarray(_attr(P, "z_grid"), float)
    if H.size != nt - 1 or z.size != nz: raise ValueError("unexpected H_own/z_grid shape")
    hRmax = float(_attr(P, "hR_max"))
    renter_limit_active = bool(np.isfinite(hRmax))
    strict_room = bool(_attr(P, "child_room_floor")) and (float(_attr(P, "hbar_child_rooms")) > 0.0 or float(_attr(P, "hbar_first_child_jump")) > 0.0)
    room_scale = float(getattr(P, "owner_h_bar_scale", 1.0))
    R, delta, tauH, psi = map(float, (_attr(P, "R_gross"), _attr(P, "delta"), _attr(P, "tau_H"), _attr(P, "psi")))
    g_hash_before = array_sha256(g)
    occupied = np.argwhere(g > 0.0)
    stay_mass = exit_mass = unresolved_mass = 0.0
    minmargin = np.inf
    unresolved: list[dict[str, Any]] = []
    for bidx, old, i, j, iz, child, state in occupied:
        mass = float(g[bidx, old, i, j, iz, child, state]) * scale
        b, cbar, hbar = float(bg[bidx]), float(cb[child, state]), float(hb[child, state])
        inc = income(P, int(i), int(j), float(z[iz]))
        if old == 0:
            floor, flow, adequate = max(float(bg[0]), 0.0), cbar + r * hbar, (not renter_limit_active or hbar < hRmax)
            x = b
        else:
            house = float(H[old - 1]); collateral = -float(phi[i, old, child, state]) * q * house
            death = -(1.0 - psi) * q * house if _death_possible(P, int(j)) else -np.inf
            floor = max(float(bg[0]), min(b, collateral), death)
            flow, adequate = cbar + (delta + tauH) * q * house, (not strict_room or house - room_scale * hbar > 0.0)
            x = b + (1.0 - psi) * q * house / R
        margin = R * b + inc - floor - flow
        minmargin = min(minmargin, margin)
        stay_ok = bool(np.isfinite(margin) and margin > TOL and adequate)
        exit_ok = False
        exit_margin = None
        if old > 0 and bg[0] <= x <= bg[-1]:
            next_floor = max(float(bg[0]), 0.0)
            exit_margin = R * x + inc - next_floor - cbar - r * hbar
            exit_ok = bool(exit_margin > TOL and hbar < hRmax)
        if stay_ok:
            stay_mass += mass
        elif exit_ok:
            exit_mass += mass
        else:
            unresolved_mass += mass
            unresolved.append({"mass": mass, "margin": margin, "exit_margin": exit_margin, "wealth": b, "origin_tenure": int(old), "location": int(i), "age_index": int(j), "income_index": int(iz), "children_ever_born_index": int(child), "child_state_index": int(state), "housing_adequate": bool(adequate)})
    if g_hash_before != array_sha256(g): raise RuntimeError("packet population array mutated")
    unresolved.sort(key=lambda x: (-x["mass"], x["margin"]))
    return {"label": "ANALYTICAL CURRENT-BUDGET witnesses only; no interpolation, Bellman, native equilibrium, choices, or population changes. Unresolved is not infeasible; downsizing/other choices are omitted.", "threshold": TOL, "occupied_count": int(len(occupied)), "population_mass": float(g.sum() * scale), "witnessed_stay_mass": stay_mass, "witnessed_exit_only_mass": exit_mass, "unresolved_mass": unresolved_mass, "unresolved_count": len(unresolved), "minmargin": float(minmargin) if len(occupied) else None, "top_unresolved_margins": unresolved[:20], "actual_flags": flags, "q": q, "r": r, "pension": pension, "scale": scale, "renter_hR_max_active": renter_limit_active, "population_array_sha256_before": g_hash_before, "population_array_sha256_after": array_sha256(g), "scaled_population_sha256": array_sha256(g * scale)}


def _receipt_scale(receipt: Path) -> float:
    obj = json.loads(receipt.read_text())
    found: list[float] = []
    def walk(x: Any) -> None:
        if isinstance(x, dict):
            for k, v in x.items():
                if k in {"population_scale", "original_scale", "scale"} and isinstance(v, (int, float)): found.append(float(v))
                walk(v)
        elif isinstance(x, list):
            for v in x: walk(v)
    walk(obj)
    exact = [x for x in found if x == EXPECTED_SCALE]
    if not exact: raise ValueError("authenticated closure receipt does not contain the required original scale")
    return exact[0]


def authenticate_cited_sources(source_root: Path) -> dict[str, str]:
    """Check the three cited frozen files against their frozen source manifest."""
    manifest = source_root / "code/model/experiments/birth_count_choice/source_manifest.json"
    document = json.loads(manifest.read_text())
    entries = document.get("files", [])
    amended = {x.get("file"): x.get("experiment_sha256") for x in document.get("engine_modified_files", [])}
    wanted = {"shared.py", "household.py", "kernels.py"}
    hashes: dict[str, str] = {}
    for entry in entries:
        copied = str(entry.get("copied", ""))
        name = Path(copied).name
        if name in wanted:
            actual = sha256(source_root / copied)
            expected = amended.get(copied, entry.get("sha256"))
            if actual != expected:
                raise ValueError(f"frozen {name} fails its authenticated source hash")
            hashes[name] = actual
    if set(hashes) != wanted:
        raise ValueError("frozen source manifest lacks a cited household/kernels/shared hash")
    return hashes


def main() -> None:
    ap = argparse.ArgumentParser(description="Inert same-tenure current-budget stay-witness screen")
    ap.add_argument("--packet", type=Path, required=True); ap.add_argument("--closure", type=Path, required=True)
    ap.add_argument("--q", type=float, required=True); ap.add_argument("--r", type=float, required=True); ap.add_argument("--pension", type=float, required=True); ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--frozen-source-root", type=Path, required=True, help="execution_smoke_v5/frozen/source")
    a = ap.parse_args()
    if sha256(a.packet) != PACKET_SHA256: raise ValueError("packet SHA-256 does not match the named inert reference packet")
    if sha256(a.closure) != CLOSURE_SHA256: raise ValueError("closure SHA-256 does not match the authenticated original receipt")
    scale = _receipt_scale(a.closure)
    source_hashes = authenticate_cited_sources(a.frozen_source_root)
    packet = load_packet(a.packet)
    out = screen(packet, q=a.q, r=a.r, pension=a.pension, scale=scale)
    out.update({"packet_sha256": sha256(a.packet), "closure_sha256": sha256(a.closure), "authenticated_frozen_source_sha256": source_hashes, "packet_arrays_sha256": {k: array_sha256(packet[k]) for k in ("b_grid", "stationary_g_pre")}})
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__": main()
