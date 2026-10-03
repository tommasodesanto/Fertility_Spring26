#!/usr/bin/env python3
"""Create a self-contained, hash-preserving stationary-GE review bundle.

The exported copy keeps the frozen observer source byte-for-byte.  Its only
portable adaptation is a small wrapper in ``review_bundle_localization.py``:
after the observer's own source-hash check has passed, it resolves the
observer's historical absolute paths beneath the unpacked bundle root.  It
does not alter numerical code, inputs, tolerances, or any SHA-256 gate.
"""
from __future__ import annotations

import argparse
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

OLD_ROOT = %r

def localize(value, root):
    if isinstance(value, str) and value.startswith(OLD_ROOT + "/"):
        return str(Path(root) / value.removeprefix(OLD_ROOT + "/"))
    if isinstance(value, list): return [localize(x, root) for x in value]
    if isinstance(value, dict): return {k: localize(v, root) for k, v in value.items()}
    return value

def install(frozen_module, root):
    root = Path(root).resolve()
    original_read = frozen_module.read
    def read_local(path): return localize(original_read(path), root)
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


def scripts() -> None:
    shell = '''#!/bin/sh
set -eu
ROOT="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
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

On Windows, run the displayed Python commands directly (the four `.sh` files
only set paths and one-core environment variables). Set `PYTHON=/path/to/python`
before a script to use another interpreter.

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
run writes a new case. `CODE_MAP.md` lists all included components.

## Frozen observer and portability

The reporting layer uses a frozen empirical/accounting observer. Its original
source and input files are included with SHA-256 verification. The historical
source names an absolute repository path; `review_bundle_localization.py`
changes only resolved paths after that source hash has been verified. It does
not edit observer equations, numerical tolerances, inputs, or acceptance gates.
`SOURCE_MANIFEST.sha256` records every bundled file and whether it was copied
byte-for-byte or is the one documented localization facade.

## Estate/birth-menu note

The baseline here has the historical gross estate/default configuration. The
isolated Estate-A experiment is **not** a default and is not runnable from this
small review bundle. In that separate experiment, net bequest is
`b' + (1 - selling_cost) * price * housing`, with selling cost 0.06 applied in
both utility and death-flow accounting, and no additional interest on `b'`.
It retained the chain-13 inputs and used birth-menu caps one and three. Its
result summary is included at `references/estate_a_RESULTS.md` for review only.

## Scientific limits

The bundle reproduces the local stationary workflow at the supplied inputs. It
does not certify calibration convergence, grid adequacy, publication readiness,
external validity, or a dynamic transition. Read `references/PROVENANCE.md`
before interpreting the cached result.
'''


def make_manifest() -> None:
    rows = []
    for path in sorted(p for p in STAGE.rglob("*") if p.is_file()):
        rel = path.relative_to(STAGE).as_posix()
        kind = "localization facade" if rel == "code/model/review_bundle_localization.py" else "copied source/input/cache"
        rows.append(f"{digest(path)}  {kind}  {rel}")
    (STAGE / "SOURCE_MANIFEST.sha256").write_text("\n".join(rows) + "\n")


def build(*, replace=False) -> dict:
    if replace and STAGE.exists():
        # The stage is exporter-owned and is rebuilt atomically as a complete bundle.
        shutil.rmtree(STAGE)
    if replace and ZIP.exists():
        ZIP.unlink()
    if STAGE.exists() or ZIP.exists():
        raise FileExistsError(f"Refusing to replace existing review bundle: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    # Production code and the minimal inspection stack.
    for rel in ("code/model/production", "code/model/parameters", "code/model/run_model.py",
                "code/model/plot_model_policies.py", "code/model/plot_model_aggregates.py",
                "code/model/tools/economics_explorer.py", "code/model/tools/model_policy_tools.py",
                "output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed",
                "output/model/publication_refactor_20260929/local_export_v1/inputs",
                "output/model/fixed_reference_economics_20260928/sources/fixed_price_v1",
                "output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/local_runtime/bootstrap.py",
                "output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/local_runtime/frozen_sources",
                "output/model/fertility_identification_20260928/fixed_reference_manifest.json",
                "output/model/fertility_identification_20260928/contract_v1",
                "output/model/fertility_identification_20260928/resume_v1/selected_export/primary",
                "output/model/overnight_calibration_20260928/contract_v1",
                "tmp/e5f_overnight_local_20260927/portable/tools_v4",
                "tmp/e5f_overnight_local_20260927/portable/night_launch_v4/primary_continuation",
                "code/model/experiments/purchase_timing_sandbox",
                "output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json",
                "output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/native_postcheck/completed.json"):
        copy_path(rel)
    # Cached baseline is copied as a real directory, never as a symlink.
    shutil.copytree(CASE, STAGE / "output/model/local_solution/latest", ignore=ignored, symlinks=False)
    config = STAGE / "output/model/local_solution/latest/explorer_cases.json"
    config.write_text(config.read_text().replace(OLD_ROOT, str(STAGE)))
    patch_driver()
    write(STAGE / "code/model/review_bundle_localization.py", LOCALIZER)
    write(STAGE / "requirements.txt", "numpy==2.2.6\nscipy==1.15.3\nnumba==0.61.2\nmatplotlib==3.10.3\n")
    write(STAGE / "README.md", readme())
    write(STAGE / "CODE_MAP.md", "# Code map\n\n- `code/model/run_model.py`: one stationary GE.\n- `code/model/parameters/`: editable baseline and toy inputs.\n- `code/model/production/`: canonical solver, equilibrium closure, reporting, storage.\n- `plot_*`: cached plot writers.\n- `tools/economics_explorer.py`: local saved-case browser.\n- `output/.../arms/indexed/`: authenticated observer facade; one documented path-only localization line is added to its copied `driver.py`.\n- `output/model/local_solution/latest/`: cached baseline outputs.\n")
    estate = ROOT / "output/model/experiments/birth_count_choice/estate_a_v1/RESULTS.md"
    if estate.is_file():
        copy_path("output/model/experiments/birth_count_choice/estate_a_v1/RESULTS.md")
        (STAGE / "references").mkdir(exist_ok=True)
        shutil.move(STAGE / "output/model/experiments/birth_count_choice/estate_a_v1/RESULTS.md", STAGE / "references/estate_a_RESULTS.md")
    write(STAGE / "references/PROVENANCE.md", "# Provenance\n\nBaseline: post-interest soft-credit chain 13, verified October 3, 2026. The complete cached target and parameter tables are in `../output/model/local_solution/latest/`. This review bundle is isolated from active calibration and cluster work.\n")
    scripts(); make_manifest()
    with zipfile.ZipFile(ZIP, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for path in sorted(p for p in STAGE.rglob("*") if p.is_file()):
            z.write(path, path.relative_to(OUT).as_posix())
    return {"stage": str(STAGE), "archive": str(ZIP), "archive_sha256": digest(ZIP),
            "archive_bytes": ZIP.stat().st_size, "files": sum(1 for p in STAGE.rglob("*") if p.is_file())}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", action="store_true", help="create the bundle")
    parser.add_argument("--replace", action="store_true", help="replace this exporter's incomplete stage only")
    args = parser.parse_args()
    if not args.build: parser.error("pass --build")
    print(json.dumps(build(replace=args.replace), indent=2))

if __name__ == "__main__": main()
