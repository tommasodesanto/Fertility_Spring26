#!/usr/bin/env python3
"""Mechanical extraction of the executed stationary path (stdlib only; Torch).

Reads the pinned, read-only source, computes the closure of top-level
definitions reachable from ENTRY, and writes `engine/<module>.py` files that
contain only reachable definitions. Kept definitions are copied as exact
source lines. The only edits are import statements:
  * name lists filtered to names that exist in the pruned target module;
  * absolute imports of in-scope modules become package-relative.
Every edit is listed in the receipt. Dynamic lookups (getattr on a module
alias, globals(), sys.modules) are flagged for review, never guessed.

    python3 materialize.py --source-root ROOT --out DIR
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import sys
from pathlib import Path

PACKAGE = "intergen_eqscale_seq_optimized"
SEARCH = ("code/model/" + PACKAGE, "code/model/tools")
ENTRY = [
    ("solver", "solve_markov_income_equilibrium"),
    ("solver", "solve_markov_income_at_prices"),
    ("solver", "precompute_shared"),
    ("solver", "make_grid"),
    ("solver", "configure_current_household_contract"),
    ("solver", "get_fecundity_by_age"),
    ("solver", "InfeasibleThetaError"),
    ("solver", "DEAD_MASS_TOL"),
    ("solver", "DEAD_VALUE_CUTOFF"),
    ("e5f_stationary_paygo", "bind_initial_balanced_pension"),
    ("e5f_stationary_paygo", "certify_initial_pension"),
    ("e5f_stationary_paygo", "solve_balanced_initial_equilibrium"),
    ("diagnostics", "write_diagnostics"),
]


def sha(text: bytes) -> str:
    return hashlib.sha256(text).hexdigest()


class Module:
    def __init__(self, name: str, path: Path):
        self.name, self.path = name, path
        self.raw = path.read_bytes()
        self.text = self.raw.decode()
        self.lines = self.text.splitlines(keepends=True)
        self.tree = ast.parse(self.text)
        self.defs: dict[str, ast.stmt] = {}
        self.aliases: dict[str, tuple[str, str | None]] = {}  # local -> (module, attr|None)
        for node in self.tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                self.defs[node.name] = node
            elif isinstance(node, (ast.Assign, ast.AnnAssign)):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                for t in targets:
                    for n in ast.walk(t):
                        if isinstance(n, ast.Name):
                            self.defs[n.id] = node

    def span(self, node: ast.stmt) -> tuple[int, int]:
        start = min([node.lineno] + [d.lineno for d in getattr(node, "decorator_list", [])])
        return start, node.end_lineno

    def segment(self, node: ast.stmt) -> str:
        a, b = self.span(node)
        return "".join(self.lines[a - 1:b])


def resolve(importer: str, node: ast.stmt, known: set[str]):
    """Yield (local, module, attr) for in-scope imports."""
    if isinstance(node, ast.Import):
        for a in node.names:
            base = a.name.split(".")[-1] if a.name.startswith(PACKAGE + ".") else a.name
            if base in known:
                yield a.asname or a.name, base, None
    elif isinstance(node, ast.ImportFrom):
        mod = node.module or ""
        if node.level >= 1 or mod == PACKAGE or mod.startswith(PACKAGE + "."):
            mod = mod.split(".")[-1] if mod and mod != PACKAGE else ""
        if mod == "":
            for a in node.names:
                if a.name in known:
                    yield a.asname or a.name, a.name, None
        elif mod in known:
            for a in node.names:
                yield a.asname or a.name, mod, a.name


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--overlay", type=Path, help="reviewed replacement modules (same names) applied over the package")
    args = ap.parse_args()
    root = args.source_root.resolve()
    paths = {}
    for rel in SEARCH:
        for p in sorted((root / rel).glob("*.py")):
            paths.setdefault(p.stem, p)
    overlay = {}
    if args.overlay:
        for p in sorted(args.overlay.resolve().glob("*.py")):
            if p.stem not in paths:
                raise SystemExit("overlay module has no base: " + p.name)
            overlay[p.stem] = dict(base=str(paths[p.stem]), base_sha256=sha(paths[p.stem].read_bytes()),
                                   overlay=str(p), overlay_sha256=sha(p.read_bytes()))
            paths[p.stem] = p
    modules: dict[str, Module] = {}

    def get(name: str) -> Module:
        if name not in modules:
            modules[name] = Module(name, paths[name])
        return modules[name]

    known = set(paths)
    need: dict[str, set[str]] = {}
    flags: list[dict] = []
    queue = list(ENTRY)
    seen_imports: set[tuple[str, str]] = set()
    while queue:
        mname, name = queue.pop()
        m = get(mname)
        if name in need.setdefault(mname, set()):
            continue
        # module-level import alias: follow it
        if name not in m.defs:
            hit = False
            for node in m.tree.body:
                for local, tmod, attr in resolve(mname, node, known):
                    if local == name:
                        hit = True
                        seen_imports.add((mname, name))
                        if attr is not None:
                            queue.append((tmod, attr))
            if not hit:
                flags.append(dict(kind="unresolved_entry", module=mname, name=name))
            continue
        need[mname].add(name)
        node = m.defs[name]
        # function-local imports
        local_alias: dict[str, tuple[str, str | None]] = {}
        for sub in ast.walk(node):
            if isinstance(sub, (ast.Import, ast.ImportFrom)):
                for local, tmod, attr in resolve(mname, sub, known):
                    local_alias[local] = (tmod, attr)
                    if attr is not None:
                        queue.append((tmod, attr))
        module_alias = {}
        for top in m.tree.body:
            for local, tmod, attr in resolve(mname, top, known):
                module_alias[local] = (tmod, attr)
        module_alias.update(local_alias)
        for sub in ast.walk(node):
            if isinstance(sub, ast.Name):
                if sub.id in m.defs and sub.id != name:
                    queue.append((mname, sub.id))
                elif sub.id in module_alias:
                    tmod, attr = module_alias[sub.id]
                    if attr is not None:
                        queue.append((tmod, attr))
                    else:
                        get(tmod)
                        need.setdefault(tmod, set())
                if sub.id in ("globals", "vars", "__import__"):
                    flags.append(dict(kind="dynamic_" + sub.id, module=mname, definition=name, line=sub.lineno))
            elif isinstance(sub, ast.Attribute) and isinstance(sub.value, ast.Name):
                al = module_alias.get(sub.value.id)
                if al and al[1] is None:
                    queue.append((al[0], sub.attr))
            elif isinstance(sub, ast.Call) and isinstance(sub.func, ast.Name) and sub.func.id == "getattr" \
                    and sub.args and isinstance(sub.args[0], ast.Name):
                al = module_alias.get(sub.args[0].id)
                if al and al[1] is None:
                    flags.append(dict(kind="getattr_on_module", module=mname, definition=name, line=sub.lineno,
                                      target=al[0], source=ast.get_source_segment(m.text, sub)))
            if isinstance(sub, ast.Attribute) and sub.attr in ("modules",) and isinstance(sub.value, ast.Name) \
                    and sub.value.id == "sys":
                flags.append(dict(kind="sys_modules", module=mname, definition=name, line=sub.lineno))
    # write
    out = args.out
    out.mkdir(parents=True, exist_ok=False)
    receipt = dict(source_root=str(root), entry=ENTRY, modules={}, flags=flags, overlay=overlay)
    for mname in sorted(need):
        m = get(mname)
        keep = need[mname]
        pieces, edits = [], []
        shown = m.path.relative_to(root) if m.path.is_relative_to(root) else m.path.name + " (reviewed overlay)"
        header = (f'"""Extracted from {shown} (sha256 {sha(m.raw)}).\n\n'
                  f"Mechanical copy by refactor_lab/materialize.py: only reachable top-level\n"
                  f"definitions, bodies byte-identical; import edits listed in the receipt.\n\"\"\"\n")
        pieces.append(header)
        kept_defs = []
        for node in m.tree.body:
            text = m.segment(node)
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                new = rewrite_import(mname, node, text, known, need, modules)
                if new != text:
                    edits.append(dict(line=node.lineno, old=text, new=new))
                if new:
                    pieces.append(new)
            elif isinstance(node, ast.Try) and all(isinstance(s, (ast.Import, ast.ImportFrom, ast.Assign,
                                                   ast.FunctionDef, ast.Expr, ast.Pass)) for s in node.body):
                pieces.append(text)  # optional-dependency guards (e.g. numba)
            elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Assign, ast.AnnAssign)):
                names = [n for n, d in m.defs.items() if d is node]
                if any(n in keep for n in names):
                    a0, _ = m.span(node)
                    seg = m.lines[a0 - 1:node.end_lineno]
                    nested = [s for s in ast.walk(node) if isinstance(s, (ast.Import, ast.ImportFrom))]
                    for s in sorted(nested, key=lambda s: -s.lineno):
                        old = "".join(m.lines[s.lineno - 1:s.end_lineno])
                        new = rewrite_import(mname, s, old, known, need, modules)
                        if new != old:
                            edits.append(dict(line=s.lineno, nested_in=names[0], old=old, new=new))
                            seg[s.lineno - a0:s.end_lineno - a0 + 1] = [new or old[:len(old) - len(old.lstrip())] + "pass\n"]
                    text = "".join(seg)
                    pieces.append("\n\n" + text)
                    kept_defs.append(dict(names=names, lines=list(m.span(node)), sha256=sha(text.encode())))
            elif isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) and node is m.tree.body[0]:
                continue  # original docstring replaced by provenance header
            elif isinstance(node, ast.If) and "__main__" in text.split("\n", 1)[0]:
                continue
            else:
                flags.append(dict(kind="toplevel_statement_dropped", module=mname, line=node.lineno,
                                  text=text[:200]))
        body = "".join(pieces).rstrip() + "\n"
        (out / f"{mname}.py").write_text(body)
        dropped = sorted(set(m.defs) - {n for d in kept_defs for n in d["names"]})
        source = str(m.path.relative_to(root)) if m.path.is_relative_to(root) else str(m.path)
        receipt["modules"][mname] = dict(source=source, source_sha256=sha(m.raw),
            output_sha256=sha(body.encode()), kept=kept_defs, dropped=dropped, import_edits=edits,
            source_lines=len(m.lines), output_lines=body.count("\n"))
    (out / "__init__.py").write_text('"""Mechanically extracted engine; see materialize_receipt.json."""\n')
    (out / "materialize_receipt.json").write_text(json.dumps(receipt, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: (v["source_lines"], v["output_lines"], len(v["kept"])) for k, v in receipt["modules"].items()}))
    print("flags", len(flags))


def rewrite_import(mname, node, text, known, need, modules) -> str:
    """Filter imported names to kept definitions; make in-scope imports relative."""
    indent = text[: len(text) - len(text.lstrip())]
    if isinstance(node, ast.Import):
        keep = []
        for a in node.names:
            base = a.name.split(".")[-1] if a.name.startswith(PACKAGE + ".") else a.name
            if base in known:
                if base in need:
                    keep.append(f"from . import {base}" + (f" as {a.asname}" if a.asname else
                                (f" as {a.name}" if a.name != base else "")))
            else:
                keep.append(f"import {a.name}" + (f" as {a.asname}" if a.asname else ""))
        if len(keep) == 1 and keep[0] == text.strip():
            return text
        return "".join(indent + k + "\n" for k in keep)
    mod = node.module or ""
    internal = node.level >= 1 or mod == PACKAGE or mod.startswith(PACKAGE + ".")
    base = (mod.split(".")[-1] if mod and mod != PACKAGE else "") if internal else mod
    if not internal and base not in known:
        return text
    if base == "":
        names = [a for a in node.names if a.name in need]
        clause = "from . import "
    else:
        if base not in need:
            return ""
        exists = set(modules[base].defs) if base in modules else set()
        names = [a for a in node.names if a.name in need[base] or
                 (a.name not in exists and a.name != "*")]
        clause = f"from .{base} import "
    if not names:
        return ""
    if node.level >= 1 and len(names) == len(node.names):
        return text  # already package-relative and nothing filtered
    joined = ", ".join(a.name + (f" as {a.asname}" if a.asname else "") for a in names)
    return f"{indent}{clause}({joined})\n"


if __name__ == "__main__":
    main()
