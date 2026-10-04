#!/usr/bin/env python3
"""Create a self-contained, hash-preserving stationary-GE review bundle.

The exported copy keeps frozen observer sources byte-for-byte. Portable
adaptations resolve historical paths, locate explorer assets and preserve the
bundled cache when publishing the first fresh result. They do not alter
equations, inputs, tolerances or any SHA-256 gate.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import shutil
import stat
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/review_bundle_20261003"
STAGE = OUT / "Fertility_Model_Review_20261003"
ZIP = OUT / "Fertility_Model_Review_20261003.zip"
CASE = ROOT / "output/model/local_solution/cases/20261003T175652812716Z_b1c72f13"
OLD_ROOT = "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
COPIED: set[str] = set()
MODIFIED: set[str] = set()


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def ignored(_dir, names):
    return {n for n in names if n in {"__pycache__", ".DS_Store", ".pytest_cache"}
            or n.endswith((".pyc", ".pyo", ".nbc", ".nbi"))}


def copy_path(relative: str) -> None:
    source, target = ROOT / relative, STAGE / relative
    if not source.exists():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        shutil.copytree(source, target, ignore=ignored, symlinks=False)
    else:
        shutil.copy2(source, target)
    copied = target.rglob("*") if target.is_dir() else [target]
    COPIED.update(p.relative_to(STAGE).as_posix() for p in copied if p.is_file())


def copy_absolute_pin(path: str) -> None:
    """Copy one strict observer pin, rejecting anything outside this checkout."""
    source = Path(path).resolve()
    try:
        relative = source.relative_to(ROOT)
    except ValueError as exc:
        raise RuntimeError(f"Pinned dependency is outside repository: {source}") from exc
    copy_path(relative.as_posix())


def copy_frozen_dependency_closure() -> None:
    """Copy only files authenticated by the stationary reporting observer.

    The observer verifies the identification contract, its listed source pins,
    and the production ancestry contract.  It does not execute calibration
    histories, so copying their directories would make this review bundle both
    misleading and unnecessarily large.
    """
    manifest_path = ROOT / "output/model/fertility_identification_20260928/fixed_reference_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    copy_path(manifest_path.relative_to(ROOT).as_posix())
    for name in ("contract", "objective", "source_manifest", "source_contract", "native_ancestry_contract"):
        copy_absolute_pin(manifest[name]["path"])
    contract = json.loads(Path(manifest["contract"]["path"]).read_text())
    for pin in contract["files"].values():
        copy_absolute_pin(pin["path"])
    # The observer also authenticates its full historical source inventory.
    # These are source files, including historical third-party Python sources;
    # they are evidence for its hash gate, not an executable environment.
    inventory = json.loads(Path(manifest["source_manifest"]["path"]).read_text())
    recovery = ROOT / "code/model/experiments/birth_count_choice/frozen_sources"
    recovered = json.loads((recovery / "mapping.json").read_text())["mapping"]
    overlay = ROOT / "output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/local_runtime/frozen_sources"
    for relative, expected in inventory["files"].items():
        source = ROOT / relative
        if not source.is_file() or digest(source) != expected:
            record = recovered.get(str(source)) or recovered.get(str(source.resolve()))
            candidate = recovery / record["recovered_path"] if record else overlay / source.name
            if not candidate.is_file() or digest(candidate) != expected:
                raise RuntimeError(f"No exact historical source for {relative}")
            source = candidate
        target = STAGE / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        COPIED.add(relative)


def copy_nested_routing_closure() -> dict[str, str]:
    """Retain exact JSON ancestry pins; map only authenticated routing records."""
    recovery = ROOT / "code/model/experiments/birth_count_choice/frozen_sources"
    recovered = json.loads((recovery / "mapping.json").read_text())["mapping"]
    pending = list(STAGE.rglob("*.json")); seen = set(); routing = {}
    def retain(path, expected):
        source = Path(path)
        relative = source.relative_to(ROOT).as_posix()
        target = STAGE / relative
        if not target.is_file() or digest(target) != expected:
            if not source.is_file() or digest(source) != expected:
                record = recovered.get(str(source)) or recovered.get(str(source.resolve()))
                source = recovery / record["recovered_path"] if record else source
            if not source.is_file() or digest(source) != expected:
                raise RuntimeError("Missing exact nested ancestry pin: " + path)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target); COPIED.add(relative)
        if target.suffix == ".json": pending.append(target)
    def pins(value):
        if isinstance(value, dict):
            if isinstance(value.get("path"), str) and value["path"].startswith(OLD_ROOT + "/") and "sha256" in value:
                retain(value["path"], value["sha256"])
            for child in value.values(): pins(child)
        elif isinstance(value, list):
            for child in value: pins(child)
    while pending:
        path = pending.pop()
        if path in seen: continue
        seen.add(path); value = json.loads(path.read_text())
        if not isinstance(value, dict): continue
        is_routing = "target_rows" not in value and any(k in value for k in
            ("source_root", "reference_root", "runtime_tools", "native_ancestry_contract"))
        if is_routing:
            for field in ("files", "base_contract", "objective", "source_manifest",
                          "native_ancestry_contract", "source_contract"):
                if field in value:
                    pins(value[field])
        if "source_root" in value and str(value["source_root"]).startswith(OLD_ROOT + "/") and "source_manifest" in value:
            inventory_path = STAGE / Path(value["source_manifest"]["path"]).relative_to(ROOT)
            inventory = json.loads(inventory_path.read_text())
            for relative, expected in inventory["files"].items():
                retain(str(Path(value["source_root"]) / relative), expected)
        # Objectives are fingerprinted as parsed objects, so never localize them.
        if is_routing:
            routing[digest(path)] = path.relative_to(STAGE).as_posix()
    # The reporting ancestor authenticates this fixed packet without running its
    # historical search. Retain exactly the files named by its lock/inventory.
    reference = ROOT / "tmp/e5f_overnight_local_20260927/portable/nightpair_20260925_v1"
    lock = json.loads((reference / "inputs/launch_lock.json").read_text())
    for name in ("inputs/launch_lock.json", "inputs/objective.json", "inputs/proposal_bank.json",
                 "inputs/source_manifest.json", "ancestor_commute.py"):
        source = reference / name
        retain(str(source), digest(source))
    for name, expected in lock["runtime_file_sha256"].items():
        retain(str(reference / name), expected)
    inventory_path = reference / "inputs/source_manifest.json"
    inventory = json.loads(inventory_path.read_text())
    for name, expected in inventory["files"].items():
        retain(str(reference / "source" / name), expected)
    routing[digest(inventory_path)] = inventory_path.relative_to(ROOT).as_posix()
    # Native setup reads the selected ancestry case before binding observers.
    # Retain its raw receipts/tables/checkpoint rather than an executable search.
    identification = json.loads((ROOT / "output/model/fertility_identification_20260928/contract_v1/contract.json").read_text())
    case = Path(identification["reference_case"])
    for name in ("receipt.json", "initial_state.pkl.gz", "target_fit.csv", "parameters.csv"):
        source = case / name
        retain(str(source), digest(source))
    tax = reference.parent / "paygo_tax_comparison_20260924/run_paygo_two_rate.py"
    retain(str(tax), lock["tax_driver_sha256"])
    return routing


def copy_source_import_closure() -> None:
    """Retain local transitive imports without copying run histories.

    Frozen pins authenticate entry points, not every module those entry points
    import. Traverse Python imports (including literal import_module calls),
    preserving the source bytes. Include the native package's Python sources
    because its configuration also loads helpers by filename at runtime.
    """
    native = ROOT / "code/model/intergen_eqscale_seq_optimized"
    for path in sorted(native.rglob("*.py")):
        if "__pycache__" not in path.parts:
            relative = path.relative_to(ROOT).as_posix()
            if relative not in COPIED:
                copy_path(relative)
    pending = list(STAGE.rglob("*.py"))
    visited: set[str] = set()
    roots = (ROOT / "code/model/tools", ROOT / "code/model")
    while pending:
        staged = pending.pop()
        relative = staged.relative_to(STAGE).as_posix()
        if relative in visited:
            continue
        visited.add(relative)
        source = ROOT / relative
        tree = ast.parse(staged.read_text(), filename=relative)
        modules = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                modules.extend((alias.name, 0) for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                modules.append((node.module or "", node.level))
                modules.extend(((node.module + "." if node.module else "") + alias.name,
                                node.level) for alias in node.names if alias.name != "*")
            elif (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                  and node.func.attr == "import_module" and node.args
                  and isinstance(node.args[0], ast.Constant)
                  and isinstance(node.args[0].value, str)):
                modules.append((node.args[0].value, 0))
        for module, level in modules:
            if level:
                bases = (source.parent.joinpath(*([".."] * (level - 1))),)
            else:
                bases = (source.parent, *roots)
            for base in bases:
                candidate = base.joinpath(*module.split(".")) if module else base
                candidates = (candidate.with_suffix(".py"), candidate / "__init__.py")
                found = next((p for p in candidates if p.is_file()), None)
                if found is None:
                    continue
                # Normalize lexical ../ components, retaining package symlink names.
                found = Path(os.path.abspath(found))
                rel = found.relative_to(ROOT).as_posix()
                if rel not in COPIED:
                    copy_path(rel)
                    pending.append(STAGE / rel)
                # Importing a submodule needs its package initializers as well.
                parent = found.parent
                while parent != ROOT and parent.is_relative_to(ROOT):
                    init = parent / "__init__.py"
                    init_rel = init.relative_to(ROOT).as_posix()
                    if init.is_file() and init_rel not in COPIED:
                        copy_path(init_rel)
                        pending.append(STAGE / init_rel)
                    parent = parent.parent
                break


def write(path: Path, text: str, executable=False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    if executable:
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


LOCALIZER = '''"""Path-only compatibility adapter for the frozen observer.

This module deliberately keeps frozen file bytes and every observer SHA-256
check intact.  It localizes paths only in parsed JSON records, after those
records have been read from their original byte-preserved files.
"""
from pathlib import Path
import hashlib
import json
import importlib.machinery

OLD_ROOT = %r

def localize(value, root):
    if isinstance(value, str) and (value == str(root) or value.startswith(str(root) + "/")):
        return value
    if value == OLD_ROOT:
        return str(root)
    if isinstance(value, str) and value.startswith(OLD_ROOT + "/"):
        return str(Path(root) / value.removeprefix(OLD_ROOT + "/"))
    if isinstance(value, list): return [localize(x, root) for x in value]
    if isinstance(value, dict):
        return {k: v if k == "actual_serialized_parameters" else localize(v, root)
                for k, v in value.items()}
    return value

def install(frozen_module, root):
    root = Path(root).resolve()
    # This allowlist binds localization to exact routing-record bytes. All file
    # SHA checks still read the untouched files; objectives are absent here.
    allowed = json.loads((root / "code/model/review_bundle_routing.json").read_text())
    original_loads = json.loads
    def loads_local(raw, *args, **kwargs):
        value = original_loads(raw, *args, **kwargs)
        data = raw.encode() if isinstance(raw, str) else bytes(raw)
        return localize(value, root) if hashlib.sha256(data).hexdigest() in allowed else value
    json.loads = loads_local
    # Frozen Python ancestors also construct historical absolute Path globals.
    # Relocate those globals after loading their unchanged, authenticated source.
    original_exec = importlib.machinery.SourceFileLoader.exec_module
    def exec_local(loader, module):
        original_exec(loader, module)
        if Path(getattr(module, '__file__', '')).is_relative_to(root):
            for name, value in tuple(vars(module).items()):
                if not name.startswith('__') and isinstance(value, Path):
                    replacement = localize(str(value), root)
                    if replacement != str(value):
                        setattr(module, name, Path(replacement))
    importlib.machinery.SourceFileLoader.exec_module = exec_local
    original_read = frozen_module.read
    # The target objective is itself fingerprinted, including textual fields.
    # Localize only routing manifests; leave objective content byte-equivalent.
    routing = ('fixed_reference_manifest.json',
               'fertility_identification_20260928/contract_v1/contract.json')
    def read_local(path):
        value = original_read(path)
        return localize(value, root) if str(path).endswith(routing) else value
    frozen_module.read = read_local
    frozen_module.ROOT = root
    frozen_module.MANIFEST = root / "output/model/fertility_identification_20260928/fixed_reference_manifest.json"
    # The source file remains byte-identical and is hash-checked by driver.py.
    return frozen_module
''' % OLD_ROOT


def patch_driver() -> None:
    """Make the unpinned facade install the documented path-only adapter."""
    path = STAGE / "output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed/driver.py"
    text = path.read_text()
    needle = "    sys.modules[spec.name] = fp\n    spec.loader.exec_module(fp)\n"
    replacement = needle + (
        "    from review_bundle_localization import install as _install_localized_frozen_observer\n"
        "    fp = _install_localized_frozen_observer(fp, context['reference_root'])\n")
    if text.count(needle) != 1:
        raise RuntimeError("Unexpected authenticated driver shape; refusing adaptation")
    path.write_text(text.replace(needle, replacement))
    COPIED.discard(path.relative_to(STAGE).as_posix())
    MODIFIED.add(path.relative_to(STAGE).as_posix())


def patch_explorer() -> None:
    """Teach the copied viewer to resolve its saved-data paths from its config."""
    path = STAGE / "code/model/tools/economics_explorer.py"
    text = path.read_text()
    needle = "    config = json.loads(config_bytes)\n    cases = {x['id']: SavedCase(x,config['common']) for x in config['cases']}\n"
    replacement = """    config = json.loads(config_bytes)
    # Bundle configurations use paths relative to this file, so extraction to
    # any directory preserves the cached viewer without changing model data.
    for item in config['cases']:
        item['arrays'] = str((config_path.parent / item['arrays']).resolve())
    if 'report_root' in config:
        config['report_root'] = str((config_path.parent / config['report_root']).resolve())
    cases = {x['id']: SavedCase(x,config['common']) for x in config['cases']}
"""
    if text.count(needle) != 1:
        raise RuntimeError("Unexpected explorer shape; refusing portability adaptation")
    path.write_text(text.replace(needle, replacement))
    COPIED.discard(path.relative_to(STAGE).as_posix())
    MODIFIED.add(path.relative_to(STAGE).as_posix())


def patch_cache_publication() -> None:
    """Retain the unpacked example before replacing latest with a fresh pointer."""
    path = STAGE / "code/model/production/storage.py"
    text = path.read_text()
    needle = "    os.replace(temporary, latest)\n"
    replacement = """    if latest.is_dir() and not latest.is_symlink():
        retained = root / 'cases' / 'bundled_reference'
        if retained.exists():
            raise FileExistsError('Refusing to replace retained bundled reference')
        latest.rename(retained)
    os.replace(temporary, latest)
"""
    if text.count(needle) != 1:
        raise RuntimeError("Unexpected cache publication shape; refusing adaptation")
    path.write_text(text.replace(needle, replacement))
    COPIED.discard(path.relative_to(STAGE).as_posix())
    MODIFIED.add(path.relative_to(STAGE).as_posix())


def scripts() -> None:
    shell = '''#!/bin/sh
set -eu
ROOT="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLCONFIGDIR="$ROOT/.matplotlib"
exec "${PYTHON:-python3}" "$ROOT/code/model/%s" "$@"
'''
    write(STAGE / "run_baseline.sh", shell % "run_model.py", True)
    write(STAGE / "plot_policies.sh", shell % "plot_model_policies.py", True)
    write(STAGE / "plot_aggregates.sh", shell % "plot_model_aggregates.py", True)
    explorer = shell % "tools/economics_explorer.py"
    explorer = explorer.replace('"$ROOT/code/model/tools/economics_explorer.py" "$@"',
        '"$ROOT/code/model/tools/economics_explorer.py" --config "$ROOT/output/model/local_solution/latest/explorer_cases.json" "$@"')
    write(STAGE / "start_explorer.sh", explorer, True)


def readme() -> str:
    return '''# Fertility model review bundle

**Review candidate:** code and cached-result inspection are included. A fresh
equilibrium solve after extraction has not passed independent portability
verification. Do not treat this package as a certified standalone reproduction.

This is a portable review copy of the October 3, 2026 stationary Python model.
It contains code, authenticated derived inputs, and a saved baseline result. It
does **not** contain survey microdata, credentials, Git history, calibration
search histories, a Python environment, or a transition model.

## Quick start

Use Python 3.13.15. Create a clean environment, install the pinned runtime,
and run from this directory:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
./run_baseline.sh                 # one fresh stationary general equilibrium
./plot_policies.sh                # saved-result plots only; no solve
./plot_aggregates.sh              # saved-result plots only; no solve
./start_explorer.sh --port 8765   # open the printed localhost URL
```

On Windows, use `py -3.13 code\\model\\run_model.py`,
`py -3.13 code\\model\\plot_model_policies.py`,
`py -3.13 code\\model\\plot_model_aggregates.py`, and
`py -3.13 code\\model\\tools\\economics_explorer.py --config
output\\model\\local_solution\\latest\\explorer_cases.json --port 8765`.
The `.sh` files only set paths and one-core environment variables. Set
`PYTHON=/path/to/python` before a script to use another interpreter.

## Change a parameter

Edit `code/model/parameters/toy_params.py`, then run
`./run_baseline.sh --params toy_params.py`. This writes a separate cache under
`output/model/experiments/toy_params/`; it does not overwrite the saved or
working baseline. `best_params.py` is the current post-interest, soft-credit
chain-13 working input, not a claim of global optimality. The parameter-file
format accepts data only and validates supported values before solving.

## What a run means

`run_baseline.sh` evaluates one fixed parameter vector. Internally it solves
for the stationary housing price that passes the birth-renewal equilibrium
gate, then reports the full 14-row target table. It is not a calibration:
calibration would repeatedly change ten parameters across many such solves.
The default fixed-`H0` closure reports the implied population scale; it is not
a dated transition or a new empirical estimate.

## Saved material and navigation

`output/model/local_solution/latest/` is a cached, validated baseline: 17
standard diagnostics, 8 policy plots, 7 aggregate plots, target and parameter
tables, and explorer arrays. These cached files are for inspection; a fresh
run first preserves them as `cases/bundled_reference/`, then writes a new case.
`CODE_MAP.md` lists all included components.

## Frozen observer and portability

The reporting layer uses a frozen empirical/accounting observer. Its original
source and input files are included with SHA-256 verification. The historical
source names an absolute repository path; `review_bundle_localization.py`
changes only resolved paths after that source hash has been verified. The copied
explorer resolves relative asset paths. On the first fresh run, the copied cache
publisher retains the bundled example under `cases/bundled_reference` before
updating `latest`. These adaptations change no equations, numerical tolerances,
inputs, or acceptance gates.
`SOURCE_MANIFEST.sha256` records every bundled file and whether it was copied
byte-for-byte or is a generated or modified portability file. The strict source
inventory includes historical third-party Python source files needed by the
observer's hash checks; these are not a usable Python environment.

## Estate/birth-menu note

The baseline here has the historical gross estate/default configuration. The
isolated Estate-A experiment is **not** a default and is not runnable from this
small review bundle. In that separate experiment, net bequest is
`b' + (1 - selling_cost) * price * housing`, with selling cost 0.06 applied in
both utility and death-flow accounting, and no additional interest on `b'`.
It retained the chain-13 inputs and used birth-menu caps one and three. Its
result summary is included at `references/estate_a_RESULTS.md` for review only.

## Scientific limits

The bundle is intended to reproduce the local stationary workflow at the supplied
inputs; the fresh-solve portability check remains incomplete. It
does not certify calibration convergence, grid adequacy, publication readiness,
external validity, or a dynamic transition. Read `references/PROVENANCE.md`
before interpreting the cached result.
'''


def make_manifest() -> None:
    rows = []
    for path in sorted(p for p in STAGE.rglob("*") if p.is_file()):
        rel = path.relative_to(STAGE).as_posix()
        if rel in COPIED:
            kind = "byte-for-byte copy"
        elif rel == "code/model/review_bundle_localization.py":
            kind = "generated path-only localization adapter"
        elif rel.endswith("explorer_cases.json"):
            kind = "generated portable cache configuration"
        elif rel in MODIFIED:
            kind = "modified copied file: documented path/cache portability adaptation"
        else:
            kind = "generated bundle file"
        rows.append(f"{digest(path)}  {kind}  {rel}")
    (STAGE / "SOURCE_MANIFEST.sha256").write_text("\n".join(rows) + "\n")


def build(*, replace=False) -> dict:
    COPIED.clear()
    MODIFIED.clear()
    if replace and STAGE.exists():
        shutil.rmtree(STAGE)
    if replace and ZIP.exists():
        ZIP.unlink()
    if STAGE.exists() or ZIP.exists():
        raise FileExistsError(f"Refusing to replace existing review bundle: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    # Production code and the minimal inspection stack.
    for rel in ("code/model/production", "code/model/parameters", "code/model/run_model.py",
                "code/model/plot_model_policies.py", "code/model/plot_model_aggregates.py",
                "code/model/tools/economics_explorer.py", "code/model/tools/economics_explorer.html",
                "code/model/tools/model_policy_tools.py",
                "output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed",
                "output/model/publication_refactor_20260929/local_export_v1/inputs",
                "output/model/fixed_reference_economics_20260928/sources/fixed_price_v1",
                "output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/local_runtime/bootstrap.py",
                "output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/local_runtime/frozen_sources",
                "output/model/fertility_identification_20260928/resume_v1/selected_export/primary",
                "tmp/e5f_overnight_local_20260927/portable/tools_v4",
                "code/model/experiments/purchase_timing_sandbox",
                "output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json",
                "output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/native_postcheck/completed.json"):
        copy_path(rel)
    copy_frozen_dependency_closure()
    routing = copy_nested_routing_closure()
    copy_source_import_closure()
    # Cached baseline is copied as a real directory, never as a symlink.
    shutil.copytree(CASE, STAGE / "output/model/local_solution/latest", ignore=ignored, symlinks=False)
    COPIED.update(p.relative_to(STAGE).as_posix() for p in (STAGE / "output/model/local_solution/latest").rglob("*") if p.is_file())
    config = STAGE / "output/model/local_solution/latest/explorer_cases.json"
    config_data = json.loads(config.read_text())
    config_data["cases"][0]["arrays"] = "explorer_arrays.npz"
    config_data["report_root"] = "."
    config.write_text(json.dumps(config_data, indent=2) + "\n")
    COPIED.discard(config.relative_to(STAGE).as_posix())
    patch_driver()
    patch_explorer()
    patch_cache_publication()
    write(STAGE / "code/model/review_bundle_routing.json", json.dumps(routing, indent=2) + "\n")
    write(STAGE / "code/model/review_bundle_localization.py", LOCALIZER)
    write(STAGE / "requirements.txt", "# Verified with Python 3.13.15\nnumpy==2.2.6\nscipy==1.15.3\nnumba==0.61.2\nmatplotlib==3.10.3\n")
    write(STAGE / "README.md", readme())
    write(STAGE / "CODE_MAP.md", "# Code map\n\n- `code/model/run_model.py`: one stationary GE.\n- `code/model/parameters/`: editable baseline and toy inputs.\n- `code/model/production/`: canonical solver, equilibrium closure, reporting, storage.\n- `plot_*`: cached plot writers.\n- `tools/economics_explorer.py`: local saved-case browser.\n- `output/.../arms/indexed/`: authenticated observer facade; one documented path-only localization line is added to its copied `driver.py`.\n- `output/model/local_solution/latest/`: cached baseline outputs.\n")
    estate = ROOT / "output/model/experiments/birth_count_choice/estate_a_v1/RESULTS.md"
    if estate.is_file():
        copy_path("output/model/experiments/birth_count_choice/estate_a_v1/RESULTS.md")
        (STAGE / "references").mkdir(exist_ok=True)
        shutil.move(STAGE / "output/model/experiments/birth_count_choice/estate_a_v1/RESULTS.md", STAGE / "references/estate_a_RESULTS.md")
        COPIED.add("references/estate_a_RESULTS.md")
    write(STAGE / "references/PROVENANCE.md", "# Provenance\n\nBaseline: post-interest soft-credit chain 13, verified October 3, 2026. The complete cached target and parameter tables are in `../output/model/local_solution/latest/`. This review bundle is isolated from active calibration and cluster work.\n")
    scripts(); make_manifest()
    with zipfile.ZipFile(ZIP, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for path in sorted(p for p in STAGE.rglob("*") if p.is_file()):
            z.write(path, path.relative_to(OUT).as_posix())
    return {"stage": str(STAGE), "archive": str(ZIP), "archive_sha256": digest(ZIP),
            "archive_bytes": ZIP.stat().st_size, "files": sum(1 for p in STAGE.rglob("*") if p.is_file())}


def finalize_verified() -> dict:
    """Publish verified documentation without rebuilding any executable bytes."""
    receipt_path = OUT / "verification.json"
    receipt = json.loads(receipt_path.read_text())
    tested_sha = digest(ZIP)
    if receipt.get("status") != "passed" or receipt.get("archive_sha256") != tested_sha:
        raise RuntimeError("Finalization requires a passed receipt for this exact ZIP")
    environment = json.loads((OUT / "test_environment.json").read_text())
    if (environment.get("created_from_declared_requirements_only") is not True or
            environment.get("requirements_sha256") != digest(STAGE / "requirements.txt") or
            environment.get("python") != receipt.get("test_python")):
        raise RuntimeError("Verified environment/requirements identity differs")
    manifest_path = STAGE / "SOURCE_MANIFEST.sha256"
    original_manifest = manifest_path.read_bytes()
    rows = {}
    for line in original_manifest.decode().splitlines():
        expected, kind, relative = line.split("  ", 2)
        if relative in rows: raise RuntimeError("Duplicate manifest path: " + relative)
        rows[relative] = (expected, kind)
    stage_files = {p.relative_to(STAGE).as_posix() for p in STAGE.rglob("*") if p.is_file()}
    if stage_files != set(rows) | {"SOURCE_MANIFEST.sha256"}:
        raise RuntimeError("Stage contains unmanifested or missing files")
    for relative, (expected, _kind) in rows.items():
        if digest(STAGE / relative) != expected:
            raise RuntimeError("Stage differs from tested source manifest: " + relative)
    with zipfile.ZipFile(ZIP) as archive:
        names = [item.filename for item in archive.infolist()]
        wanted = {STAGE.name + "/" + name for name in stage_files}
        if len(names) != len(set(names)) or set(names) != wanted:
            raise RuntimeError("ZIP file inventory differs from stage")
        for relative in stage_files:
            h = hashlib.sha256()
            with archive.open(STAGE.name + "/" + relative) as stream:
                for block in iter(lambda: stream.read(1 << 20), b""): h.update(block)
            if h.hexdigest() != digest(STAGE / relative):
                raise RuntimeError("ZIP bytes differ from stage: " + relative)
    text = (STAGE / "README.md").read_text()
    candidate = "**Review candidate:** code and cached-result inspection are included. A fresh\nequilibrium solve after extraction has not passed independent portability\nverification. Do not treat this package as a certified standalone reproduction."
    verified = "**Verified current working baseline:** independent clean extraction passed cached\nplots, the explorer, and a fresh stationary equilibrium with exact agreement on\nall 14 target rows, 31 parameter rows, and stored arrays. Original-project file\naccess was blocked. See `VERIFICATION.json`. This is the post-interest, soft-credit\nchain-13 working baseline, not a certified paper calibration or global optimum."
    if text.count(candidate) != 1:
        raise RuntimeError("Unexpected candidate README; refusing finalization")
    text = text.replace(candidate, verified)
    text = text.replace(
        "The bundle is intended to reproduce the local stationary workflow at the supplied\ninputs; the fresh-solve portability check remains incomplete. It\n",
        "The bundle reproduces the tested stationary workflow at the supplied inputs. It\n")
    start = text.index("On Windows, use ")
    end = text.index("## Change a parameter", start)
    text = text[:start] + "The fresh solve was tested on macOS with Python 3.13.15 using one core\nand an environment installed only from `requirements.txt`. Windows execution\nis untested; the implementation includes POSIX dependencies. Set\n`PYTHON=/path/to/python` before a script to use another interpreter.\n\n" + text[end:]
    public = {key: receipt[key] for key in
        ("status", "full_target_rows_exact", "full_parameter_rows_exact", "stored_arrays_exact",
         "fresh_ge_original_project_access_blocked", "cached_reference_preserved",
         "production_defaults_unchanged", "transition_not_included", "cached_plotters", "explorer",
         "fresh_figures", "clean_dependency_environment_verified")
        if key in receipt}
    public.update(tested_archive_sha256=tested_sha, baseline="October 3, 2026 post-interest soft-credit chain 13",
                  requirements_sha256=environment["requirements_sha256"],
                  python="3.13.15", platform_tested="macOS", cores=1,
                  environment="Installed only from requirements.txt", windows_execution="untested")
    changed = ["README.md", "SOURCE_MANIFEST.sha256", "VERIFICATION.json"]
    replacements = {"README.md": text.encode(),
                    "VERIFICATION.json": (json.dumps(public, indent=2) + "\n").encode()}
    updated = dict(rows)
    for relative, data in replacements.items():
        updated[relative] = (hashlib.sha256(data).hexdigest(), "generated verified bundle documentation")
    replacements["SOURCE_MANIFEST.sha256"] = ("\n".join(
        f"{expected}  {kind}  {relative}" for relative, (expected, kind) in sorted(updated.items())) + "\n").encode()
    temporary = ZIP.with_name(ZIP.name + ".finalizing")
    try:
        with zipfile.ZipFile(ZIP) as source, zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as dest:
            for item in source.infolist():
                relative = item.filename.removeprefix(STAGE.name + "/")
                if relative in replacements:
                    dest.writestr(item, replacements[relative])
                else:
                    with source.open(item) as incoming, dest.open(item, "w") as outgoing:
                        shutil.copyfileobj(incoming, outgoing, 1 << 20)
            dest.writestr(STAGE.name + "/VERIFICATION.json", replacements["VERIFICATION.json"])
        # No input/source has been recopied. Assert the exact documentation delta.
        with zipfile.ZipFile(ZIP) as before, zipfile.ZipFile(temporary) as after:
            differences = []
            for name in set(before.namelist()) | set(after.namelist()):
                if name not in before.namelist() or name not in after.namelist():
                    differences.append(name.removeprefix(STAGE.name + "/")); continue
                if before.getinfo(name).CRC != after.getinfo(name).CRC or before.getinfo(name).file_size != after.getinfo(name).file_size:
                    differences.append(name.removeprefix(STAGE.name + "/"))
            if sorted(differences) != changed:
                raise RuntimeError("Unexpected finalization delta: " + repr(differences))
        with zipfile.ZipFile(temporary) as final_archive:
            for relative, (expected, _kind) in updated.items():
                h = hashlib.sha256()
                with final_archive.open(STAGE.name + "/" + relative) as stream:
                    for block in iter(lambda: stream.read(1 << 20), b""): h.update(block)
                if h.hexdigest() != expected:
                    raise RuntimeError("Final ZIP byte identity differs: " + relative)
            if final_archive.read(STAGE.name + "/SOURCE_MANIFEST.sha256") != replacements["SOURCE_MANIFEST.sha256"]:
                raise RuntimeError("Final ZIP manifest differs")
        final_sha = digest(temporary)
        for relative, data in replacements.items():
            target = STAGE / relative
            pending = target.with_name(target.name + ".finalizing")
            pending.write_bytes(data); os.replace(pending, target)
        os.replace(temporary, ZIP)
        receipt.update(tested_archive_sha256=tested_sha, archive_sha256=final_sha,
                       documentation_only_changed=changed)
        pending_receipt = receipt_path.with_name(receipt_path.name + ".finalizing")
        pending_receipt.write_text(json.dumps(receipt, indent=2) + "\n")
        os.replace(pending_receipt, receipt_path)
        return {"archive_sha256": final_sha, "tested_archive_sha256": tested_sha,
                "documentation_only_changed": changed, "archive_bytes": ZIP.stat().st_size}
    finally:
        if temporary.exists(): temporary.unlink()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", action="store_true", help="create the bundle")
    parser.add_argument("--replace", action="store_true", help="replace this exporter's incomplete stage only")
    parser.add_argument("--finalize-verified", action="store_true", help="finalize documentation after exact-ZIP verification")
    args = parser.parse_args()
    if args.build == args.finalize_verified: parser.error("pass exactly one of --build or --finalize-verified")
    if args.finalize_verified and args.replace: parser.error("--replace applies only to --build")
    print(json.dumps(finalize_verified() if args.finalize_verified else build(replace=args.replace), indent=2))

if __name__ == "__main__": main()
