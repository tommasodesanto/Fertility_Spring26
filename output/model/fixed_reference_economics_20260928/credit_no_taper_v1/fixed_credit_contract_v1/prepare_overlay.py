#!/usr/bin/env python3
"""Create the fixed unsecured-credit overlay from the pinned source bundle.

This is deliberately a copier-plus-small-diff, not an importer of the model.
It rejects a changed source bundle and never overwrites an output directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


LABEL = "2007 stationary reference — block0506, September 28 verified export"
HASHES = {
    "parameters.py": "66f86697c2c58ca3864305bf13dd2be71a008905b2beb573f1a4ebafabef5464",
    "solver.py": "b637a655a9344b63f4461ee0fa4796c04bd98188477c4e6ace2c48ae0fc8aec1",
    "kernels.py": "639c9a21797dbc9f2a0e9a891f283c115353c2edfcb89c959a7fe9f32b86ca27",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_once(text: str, old: str, new: str, name: str) -> str:
    if text.count(old) != 1:
        raise RuntimeError(f"{name}: expected one matching edit anchor, found {text.count(old)}")
    return text.replace(old, new)


def repo_root(here: Path) -> Path:
    for parent in (here, *here.parents):
        if (parent / ".git").exists() and (parent / "code/model").exists():
            return parent
    raise RuntimeError("could not locate repository root")


def patch_parameters(text: str) -> str:
    text = replace_once(text, '    "use_tenure_kernel", "w_fixed",\n', '    "use_tenure_kernel", "w_fixed",\n    "unsecured_credit_limit",\n', "dynamic override registration")
    text = replace_once(text, "    P.lambda_d = 0.0\n", "    P.lambda_d = 0.0\n    # None preserves the legacy rollover/taper rule.  A scalar D activates\n    # the isolated renter-only floor b' >= -D.\n    P.unsecured_credit_limit = None\n", "parameters default")
    anchor = "    if \"eta_supply\" in od:\n"
    insert = "    raw_credit = getattr(P, \"unsecured_credit_limit\", None)\n    if raw_credit is not None:\n        if not np.isscalar(raw_credit):\n            raise ValueError(\"unsecured_credit_limit must be None or a finite non-negative scalar\")\n        try:\n            credit = float(raw_credit)\n        except (TypeError, ValueError) as exc:\n            raise ValueError(\"unsecured_credit_limit must be None or a finite non-negative scalar\") from exc\n        if not np.isfinite(credit) or credit < 0.0:\n            raise ValueError(\"unsecured_credit_limit must be None or a finite non-negative scalar\")\n        P.unsecured_credit_limit = credit\n    if raw_credit is not None and bool(getattr(P, \"native_solvency_credit\", False)):\n        raise ValueError(\"unsecured_credit_limit cannot be combined with native_solvency_credit\")\n"
    text = replace_once(text, anchor, insert + anchor, "parameters validation")
    cap = "    lam = float(getattr(P, \"lambda_d\", 0.0))\n"
    checked = "    raw_credit = getattr(P, \"unsecured_credit_limit\", None)\n    if raw_credit is not None:\n        if not np.isscalar(raw_credit):\n            raise ValueError(\"unsecured_credit_limit must be None or a finite non-negative scalar\")\n        try:\n            credit = float(raw_credit)\n        except (TypeError, ValueError) as exc:\n            raise ValueError(\"unsecured_credit_limit must be None or a finite non-negative scalar\") from exc\n        if not np.isfinite(credit) or credit < 0.0:\n            raise ValueError(\"unsecured_credit_limit must be None or a finite non-negative scalar\")\n        P.unsecured_credit_limit = credit\n    lam = float(getattr(P, \"lambda_d\", 0.0))\n"
    return replace_once(text, cap, checked, "direct rebuild validation")


def patch_kernels(text: str) -> str:
    text = replace_once(text, "    transaction_support=False,\n):\n    # Discrete tenure-choice", "    transaction_support=False,\n    require_owner_sale_solvency=False,\n):\n    # Discrete tenure-choice", "deterministic tenure signature")
    text = replace_once(text, "                            v0 = _interp_with_clip(b_grid, Vd[:, 0, id_, nn, cs], ba, strict_interpolated_support, transaction_support)\n                        if v0 > best_v:", "                            v0 = _interp_with_clip(b_grid, Vd[:, 0, id_, nn, cs], ba, strict_interpolated_support, transaction_support)\n                            if require_owner_sale_solvency and bg_b + sp < 0.0:\n                                v0 = NEG_INF\n                        if v0 > best_v:", "deterministic raw sale gate")
    text = replace_once(text, "    transaction_support=False,\n):\n    Nb, nt, I, npar, ncs = Vd.shape", "    transaction_support=False,\n    require_owner_sale_solvency=False,\n):\n    Nb, nt, I, npar, ncs = Vd.shape", "logit tenure signature")
    text = replace_once(text, "                            v0 = _interp_with_clip(b_grid, Vd[:, 0, id_, nn, cs], ba, False, transaction_support)\n                        vals[0] = v0", "                            v0 = _interp_with_clip(b_grid, Vd[:, 0, id_, nn, cs], ba, False, transaction_support)\n                            if require_owner_sale_solvency and bg_b + sp < 0.0:\n                                v0 = NEG_INF\n                        vals[0] = v0", "logit raw sale gate")
    text = replace_once(text, "    natural_floor_v=None,\n):\n", "    natural_floor_v=None,\n    fixed_renter_floor=-np.inf,\n):\n", "renter kernel signature")
    text = replace_once(text, "            unsecured_floor = rollover_floor if rollover_floor < line_floor else line_floor\n            if natural_floor_v is not None:", "            unsecured_floor = rollover_floor if rollover_floor < line_floor else line_floor\n            if np.isfinite(fixed_renter_floor):\n                unsecured_floor = fixed_renter_floor\n            if natural_floor_v is not None:", "renter kernel floor")
    return text


def patch_solver(text: str) -> str:
    old = "def renter_borrowing_floor(P: SimpleNamespace, b: Any, j: int) -> np.ndarray:\n    \"\"\"Renter floor; all renter debt is unsecured.\"\"\"\n\n    return debt_rule_at_age(P, b, j)\n"
    new = "def fixed_unsecured_credit_active(P: SimpleNamespace) -> bool:\n    return getattr(P, \"unsecured_credit_limit\", None) is not None\n\n\ndef renter_borrowing_floor(P: SimpleNamespace, b: Any, j: int) -> np.ndarray:\n    \"\"\"Renter floor; scalar credit is separate from the legacy rollover rule.\"\"\"\n\n    credit = getattr(P, \"unsecured_credit_limit\", None)\n    if credit is None:\n        return debt_rule_at_age(P, b, j)\n    floor = -float(credit)\n    death_possible = int(j) == int(P.J) - 1 or (\n        bool(getattr(P, \"use_age_survival\", False))\n        and float(np.asarray(P.survival_probs)[int(j)]) < 1.0\n    )\n    # Estates remain non-negative: a death branch makes the effective floor\n    # max(-D, 0), independently of the age-taper arrays.\n    if death_possible:\n        floor = max(floor, 0.0)\n    return np.zeros_like(np.asarray(b, dtype=float)) + floor\n"
    text = replace_once(text, old, new, "solver renter floor")
    core_anchor = "    stored_bp: np.ndarray | None,\n    eval_mode: bool,\n):\n    if fecundity_active(P)"
    core_insert = "    stored_bp: np.ndarray | None,\n    eval_mode: bool,\n):\n    if fixed_unsecured_credit_active(P):\n        raise ValueError(\"unsecured_credit_limit requires solve_bellman_full_markov_income; legacy/factored Bellman routes are unsupported\")\n    if fecundity_active(P)"
    text = replace_once(text, core_anchor, core_insert, "core scalar guard")
    text = replace_once(text, "    natural_credit = bool(getattr(P, \"native_solvency_credit\", False))\n", "    natural_credit = bool(getattr(P, \"native_solvency_credit\", False))\n    fixed_credit = fixed_unsecured_credit_active(P)\n    if fixed_credit and natural_credit:\n        raise ValueError(\"unsecured_credit_limit cannot be combined with native_solvency_credit\")\n", "saving conflict")
    old_call = "                bool(getattr(P, \"native_exact_allocation_output\", False)),\n                natural_floor,\n            )"
    new_call = "                bool(getattr(P, \"native_exact_allocation_output\", False)),\n                natural_floor,\n                float(renter_floor[0]) if fixed_credit else -np.inf,\n            )"
    text = replace_once(text, old_call, new_call, "full renter fixed floor call")
    text = replace_once(text, "    transaction_support = bool(getattr(P, \"native_purchase_income\", False))\n", "    transaction_support = bool(getattr(P, \"native_purchase_income\", False))\n    require_owner_sale_solvency = fixed_unsecured_credit_active(P)\n", "tenure sale gate flag")
    text = replace_once(text, "Vd_stay, transaction_support\n        )", "Vd_stay, transaction_support, require_owner_sale_solvency\n        )", "logit native tenure call")
    text = replace_once(text, "Vd_stay, False, transaction_support\n        )", "Vd_stay, False, transaction_support, require_owner_sale_solvency\n        )", "deterministic native tenure call")
    old_fallback = "                    Vopt[:, :, :, 0] = interp_on_grid(b_grid, Vd[:, 0, id_, :, :], ba)\n                for tn in range(1, nt):"
    new_fallback = "                    Vopt[:, :, :, 0] = interp_on_grid(b_grid, Vd[:, 0, id_, :, :], ba)\n                    if require_owner_sale_solvency:\n                        Vopt[b_grid + sp < 0.0, :, :, 0] = -1e10\n                for tn in range(1, nt):"
    text = replace_once(text, old_fallback, new_fallback, "primary Python sale gate")
    # The legacy Markov path has separate calls and a separate Python fallback.
    text = replace_once(text, "Vd, b_grid, heq, hcost, dp_choice, bmo, SD.birth_dp, birth_entry_grant, tenure_choice_kappa, Vd\n            )", "Vd, b_grid, heq, hcost, dp_choice, bmo, SD.birth_dp, birth_entry_grant, tenure_choice_kappa, Vd, False, fixed_unsecured_credit_active(P)\n            )", "markov logit gate")
    text = replace_once(text, "Vd, b_grid, heq, hcost, dp_choice, bmo, SD.birth_dp, birth_entry_grant, Vd\n            )", "Vd, b_grid, heq, hcost, dp_choice, bmo, SD.birth_dp, birth_entry_grant, Vd, False, False, fixed_unsecured_credit_active(P)\n            )", "markov deterministic gate")
    # In the Markov fallback, apply the same raw pre-clipping sale test.
    marker = "                    else:\n                        ba = np.clip(b_grid + sp, b_grid[0], b_grid[-1])\n                        Vopt[:, :, :, 0] = interp_on_grid(b_grid, Vd[:, 0, id_, :, :], ba)\n"
    replacement = marker + "                    if fixed_unsecured_credit_active(P) and to > 0:\n                        Vopt[b_grid + sp < 0.0, :, :, 0] = -1e10\n"
    text = replace_once(text, marker, replacement, "markov Python sale gate")
    return text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    destination = args.destination.resolve()
    if destination.exists():
        raise SystemExit(f"refusing to overwrite existing destination: {destination}")
    root = repo_root(Path(__file__).resolve())
    source = root / "code/model/intergen_eqscale_seq_optimized"
    before: dict[str, str] = {}
    texts: dict[str, str] = {}
    for name, expected in HASHES.items():
        path = source / name
        actual = sha(path)
        if actual != expected:
            raise SystemExit(f"hash mismatch for {path}: {actual} != {expected}")
        before[name] = actual
        texts[name] = path.read_text()
    texts["parameters.py"] = patch_parameters(texts["parameters.py"])
    texts["kernels.py"] = patch_kernels(texts["kernels.py"])
    texts["solver.py"] = patch_solver(texts["solver.py"])
    destination.mkdir(parents=True)
    for name, text in texts.items():
        (destination / name).write_text(text)
    after = {name: sha(destination / name) for name in HASHES}
    (destination / "manifest_before.json").write_text(json.dumps({"reference": LABEL, "sha256": before}, indent=2) + "\n")
    (destination / "manifest_after.json").write_text(json.dumps({"reference": LABEL, "sha256": after}, indent=2) + "\n")


if __name__ == "__main__":
    main()
