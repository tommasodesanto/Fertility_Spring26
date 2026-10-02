#!/usr/bin/env python3
"""ONE intentional algorithmic change: indexed exhaustive saving (experimental).

Writes a separate engine stage to --out (never edits engine/ in place):
every engine file copied unchanged except kernels.py, where the top-level
`exhaustive_saving_scalar` definition is replaced by the body of
`saving.exhaustive_saving_indexed` renamed to `exhaustive_saving_scalar` (so
all engine call sites are untouched), preceded by its four helpers
(`_rank`, `_interp_ranked`, `_renter_value`, `_owner_value`) copied verbatim
from saving.py. Formulas, candidate set, order and strict `>` tie rule are
those of saving.py, already tested bit-identical against the oracle.

Checks (failure stops; no fallback):
  * the engine's original exhaustive_saving_scalar equals saving.py's verbatim copy
  * helper names do not collide with engine kernels definitions
  * every other top-level kernels.py segment and every other file is byte-identical
Outputs indexed_saving.diff (one-function unified diff) and transform_receipt.json.
Remains experimental until a fresh 113-path fixed-price certificate passes.

    python3 make_indexed_stage.py --engine engine --saving saving.py --out STAGE/engine
"""
from __future__ import annotations

import argparse
import ast
import difflib
import hashlib
import json
import shutil
from pathlib import Path

TARGET = "exhaustive_saving_scalar"
REPLACEMENT = "exhaustive_saving_indexed"
HELPERS = ("_rank", "_interp_ranked", "_renter_value", "_owner_value")


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def segments(text: str) -> dict:
    lines, out = text.splitlines(keepends=True), {}
    for node in ast.parse(text).body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            start = min([node.lineno] + [d.lineno for d in node.decorator_list])
            out[node.name] = (start, node.end_lineno, "".join(lines[start - 1:node.end_lineno]))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", type=Path, required=True)
    ap.add_argument("--saving", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    kernels_text = (a.engine / "kernels.py").read_text()
    saving_text = a.saving.read_text()
    ker, sav = segments(kernels_text), segments(saving_text)
    if ker[TARGET][2] != sav[TARGET][2]:
        raise SystemExit("engine exhaustive_saving_scalar differs from saving.py verbatim copy; stop")
    clash = [h for h in HELPERS if h in ker]
    if clash:
        raise SystemExit(f"helper names collide with engine kernels: {clash}")
    new_fn = sav[REPLACEMENT][2]
    if new_fn.count(f"def {REPLACEMENT}(") != 1:
        raise SystemExit("unexpected replacement definition text")
    new_fn = new_fn.replace(f"def {REPLACEMENT}(", f"def {TARGET}(", 1)
    block = "".join(sav[h][2] + "\n\n" for h in HELPERS) + new_fn
    start, end, old_seg = ker[TARGET]
    lines = kernels_text.splitlines(keepends=True)
    new_text = "".join(lines[:start - 1]) + block + "".join(lines[end:])
    # Verify: all other top-level kernels definitions byte-identical and present once.
    after = segments(new_text)
    changed = sorted(n for n in ker if n != TARGET and after.get(n, (0, 0, None))[2] != ker[n][2])
    if changed:
        raise SystemExit(f"unintended kernels changes: {changed}")
    if sorted(set(after) - set(ker)) != sorted(HELPERS):
        raise SystemExit("unexpected new definitions")
    a.out.mkdir(parents=True, exist_ok=False)
    files = {}
    for f in sorted(a.engine.iterdir()):
        if f.is_file() and f.name != "kernels.py" and f.suffix in (".py", ".json", ""):
            shutil.copy2(f, a.out / f.name)
            files[f.name] = dict(sha256=sha(f.read_text()), unchanged=True)
    (a.out / "kernels.py").write_text(new_text)
    files["kernels.py"] = dict(sha256_before=sha(kernels_text), sha256_after=sha(new_text), unchanged=False)
    diff = "".join(difflib.unified_diff(kernels_text.splitlines(keepends=True), new_text.splitlines(keepends=True),
                                        "engine/kernels.py", "indexed/kernels.py"))
    (a.out / "indexed_saving.diff").write_text(diff)
    (a.out / "transform_receipt.json").write_text(json.dumps(dict(
        status="experimental_until_113_path_fixed_price_certificate",
        change=f"kernels.{TARGET} body := saving.{REPLACEMENT} (renamed) + helpers {list(HELPERS)}",
        saving_sha256=sha(saving_text), replaced_segment_sha256=sha(old_seg), new_block_sha256=sha(block),
        other_kernels_definitions_unchanged=len(ker) - 1, files=files), indent=1) + "\n")
    print(json.dumps(dict(diff_lines=diff.count("\n"), unchanged_definitions=len(ker) - 1,
                          unchanged_files=sum(v["unchanged"] for v in files.values()))))


if __name__ == "__main__":
    main()
