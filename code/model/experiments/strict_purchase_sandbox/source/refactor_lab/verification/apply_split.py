#!/usr/bin/env python3
"""Prepare the shared/household/distribution/equilibrium split of engine/solver.py.

READY, NOT APPLIED: run only after lead_replay_status.json records a passing
unchanged-source fixed-price replay. Writes a complete engine copy to --out
(never edits engine/ in place):

  shared.py, household.py, distribution.py, equilibrium.py
      top-level definitions from engine/solver.py in their original order,
      each body byte-identical; every module repeats solver.py's own import
      block and imports lower-layer names explicitly (plan has 0 upward refs)
  solver.py
      facade re-exporting every definition by name, so `from . import solver
      as model` callers (joint_nested, paygo, run.py) are unchanged
  split_receipt.json
      per-definition source sha256 and destination module; checked on write

    python3 apply_split.py --engine engine --plan split_plan.json --out DIR
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import shutil
from pathlib import Path

ORDER = ("shared", "household", "distribution", "equilibrium")


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", type=Path, required=True)
    ap.add_argument("--plan", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    plan = json.loads(a.plan.read_text())
    if plan["upward_references"]:
        raise SystemExit("plan has upward references")
    text = (a.engine / "solver.py").read_text()
    lines = text.splitlines(keepends=True)
    tree = ast.parse(text)
    stage = {name: s for s, names in plan["assignment"].items() for name in names}
    imports, nodes = [], []
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom, ast.Try)):
            start = min([node.lineno] + [d.lineno for d in getattr(node, "decorator_list", [])])
            imports.append("".join(lines[start - 1:node.end_lineno]))
        elif isinstance(node, (ast.FunctionDef, ast.ClassDef, ast.Assign, ast.AnnAssign)):
            nodes.append(node)
        elif not (isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant)):
            raise SystemExit(f"unexpected top-level statement at line {node.lineno}")
    def names_of(node):
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            return [node.name]
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        return [n.id for t in targets for n in ast.walk(t) if isinstance(n, ast.Name)]
    owner, segments, receipt = {}, {s: [] for s in ORDER}, []
    for node in nodes:
        names = names_of(node)
        dest = {stage[n] for n in names}
        if len(dest) != 1:
            raise SystemExit(f"definition split across modules: {names}")
        dest = dest.pop()
        start = min([node.lineno] + [d.lineno for d in getattr(node, "decorator_list", [])])
        seg = "".join(lines[start - 1:node.end_lineno])
        segments[dest].append((names, seg, node))
        receipt.append(dict(names=names, module=dest, source_lines=[start, node.end_lineno], sha256=sha(seg)))
        for n in names:
            owner[n] = dest
    if set(owner) != set(stage):
        raise SystemExit("plan and solver definitions differ")
    a.out.mkdir(parents=True, exist_ok=False)
    for f in a.engine.iterdir():
        if f.suffix in (".py", ".json") and f.name != "solver.py":
            shutil.copy2(f, a.out / f.name)
    header = "".join(imports)
    for level, mod in enumerate(ORDER):
        used = set()
        for _, _, node in segments[mod]:
            used |= {n.id for n in ast.walk(node) if isinstance(n, ast.Name) and n.id in owner}
        lower = {}
        for n in sorted(used):
            src = owner[n]
            if ORDER.index(src) > level:
                raise SystemExit(f"upward reference {mod} -> {src}:{n}")
            if src != mod:
                lower.setdefault(src, []).append(n)
        body = [f'"""{mod.capitalize()} stage of the extracted solver (bodies byte-identical; see split_receipt.json)."""\n',
                header]
        for src in ORDER:
            if src in lower:
                body.append(f"from .{src} import ({', '.join(lower[src])})\n")
        body += ["\n\n" + seg for _, seg, _ in segments[mod]]
        (a.out / f"{mod}.py").write_text("".join(body).rstrip() + "\n")
    facade = ['"""Solver facade: re-exports every definition from the four stage modules."""\n']
    for mod in ORDER:
        names = [n for names, _, _ in segments[mod] for n in names]
        facade.append(f"from .{mod} import ({', '.join(names)})  # noqa: F401\n")
    (a.out / "solver.py").write_text("".join(facade))
    # Check: every definition present byte-identically, in original order within its module.
    for mod in ORDER:
        out_text = (a.out / f"{mod}.py").read_text()
        pos = -1
        for names, seg, _ in segments[mod]:
            i = out_text.find(seg)
            if i <= pos:
                raise SystemExit(f"segment missing or reordered: {names}")
            pos = i
        ast.parse(out_text)
    (a.out / "split_receipt.json").write_text(json.dumps(dict(
        source=str(a.engine / "solver.py"), source_sha256=sha(text), definitions=receipt,
        modules={m: dict(sha256=sha((a.out / f"{m}.py").read_text()),
                         lines=(a.out / f"{m}.py").read_text().count("\n")) for m in ORDER + ("solver",)}),
        indent=1) + "\n")
    print(json.dumps({m: (a.out / f"{m}.py").read_text().count("\n") for m in ORDER + ("solver",)}))


if __name__ == "__main__":
    main()
