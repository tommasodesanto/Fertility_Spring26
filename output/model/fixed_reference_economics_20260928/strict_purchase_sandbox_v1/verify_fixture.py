"""Zero-solve source and purchase-eligibility fixture for this experiment."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT / "code/model/experiments/strict_purchase_sandbox/source"
ORIGINS = {
    "refactor_lab": ROOT / "code/model/refactor_lab",
    "small_credit_lab": ROOT / "output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed/source/small_credit_lab",
}
OLD = """            if purchase_income:
                income_for_purchase = np.array([
                    income_at_state(P, i, j, float(z_value)) for i in range(I)
                ], dtype=float).reshape(I, 1, 1, 1) / Rg
                dp_choice = ctx.dp_arr - income_for_purchase
                bmo_purchase = np.maximum(ctx.bmo - income_for_purchase, b_grid[0])
"""
NEW = """            if purchase_income:
                # Experimental origination: only wealth held at period start is eligible.
                bmo_purchase = np.maximum(ctx.bmo, b_grid[0])
"""


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    manifest = json.loads((HERE / "manifest.json").read_text())
    assert sha(HERE / "incumbent.json") == manifest["incumbent_sha256"]
    for pair, hashes in manifest["source_pairs"].items():
        original, sandbox = (ROOT / rel for rel in pair.split("|"))
        assert sha(original) == hashes["original"], original
        assert sha(sandbox) == hashes["sandbox"], sandbox
    for name, origin in ORIGINS.items():
        copy = SOURCE / name
        files = {p.relative_to(origin) for p in origin.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc" and p.name != ".DS_Store"}
        copied = {p.relative_to(copy) for p in copy.rglob("*") if p.is_file()}
        assert copied == files, name
        for rel in files:
            before = (origin / rel).read_bytes()
            after = (copy / rel).read_bytes()
            if rel == Path("engine/household.py"):
                assert before.count(OLD.encode()) == 1
                assert after == before.replace(OLD.encode(), NEW.encode())
            else:
                assert after == before, (name, rel)
            if rel.suffix == ".py":
                compile(after, str(copy / rel), "exec")
    # Illustrative state: S=2, Q=10, financed share phi=.8.
    # Income y=100 cannot rescue b=-1; with b=0 the same purchase is eligible.
    S, Q, phi = 2.0, 10.0, 0.8
    def eligible(b: float, y: float) -> bool:
        threshold = (1.0 - phi) * Q
        transaction_balance = b + S - Q
        return b + S >= threshold and transaction_balance >= -phi * Q
    assert not eligible(-1.0, 100.0)
    assert eligible(0.0, 0.0) and eligible(0.0, 100.0)
    assert eligible(2.0, 0.0) and eligible(2.0, 100.0)
    # Forward wealth and budget stay b+S-Q and b'=R(b+S-Q)+y-c-K.
    b, R, y, c, K = 2.0, 1.03, 4.0, 1.0, 0.5
    assert b + S - Q == -6.0
    assert R * (b + S - Q) + y - c - K == -3.6799999999999997
    print("PASS: exact one-block edit in both copies; source pins and equation fixture")


if __name__ == "__main__":
    main()
