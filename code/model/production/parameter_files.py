"""Load explicit, editable data files without importing or running the model."""
from __future__ import annotations

import ast
import copy
import hashlib
import json
import math
from pathlib import Path


MODEL_ROOT = Path(__file__).resolve().parents[1]
PARAMETER_ROOT = MODEL_ROOT / 'parameters'
PROJECT_ROOT = MODEL_ROOT.parents[1]
BEST_FILE = PARAMETER_ROOT / 'best_params.py'
_REQUIRED = {'PARAMETERS', 'EXTERNAL_INPUTS', 'NATIVE_OVERRIDES',
             'PRICE_GUESS', 'CLOSURE', 'BUDGET_SECONDS'}
_ALLOWED = _REQUIRED | {'PROVENANCE'}
# Allow arithmetic and literal repetition, but no imports, calls, attributes,
# comprehensions, assignments to objects, or executable control flow.
_EXPR = (ast.Constant, ast.Dict, ast.List, ast.Tuple, ast.Set, ast.BinOp,
         ast.UnaryOp, ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.USub,
         ast.UAdd, ast.Load)


def resolve_parameter_file(path='best_params.py'):
    path = Path(path).expanduser()
    if not path.is_absolute():
        path = PARAMETER_ROOT / path
    path = path.resolve(strict=True)
    if not path.is_file() or path.suffix != '.py':
        raise ValueError('parameter file must be an existing .py file')
    return path


def load_parameter_file(path='best_params.py'):
    """Freshly read a data-only Python file, validate controls, and return copies.

    Relative names always refer to code/model/parameters, independent of cwd.
    Validation constructs inputs only; it never imports an equilibrium solver.
    """
    from .inputs import DEFAULT_PARAMETERS, load_inputs
    source = resolve_parameter_file(path)
    raw = source.read_bytes()
    tree = ast.parse(raw, filename=str(source))
    names = set()
    for index, statement in enumerate(tree.body):
        if index == 0 and isinstance(statement, ast.Expr) and isinstance(statement.value, ast.Constant) and isinstance(statement.value.value, str):
            continue
        if not isinstance(statement, ast.Assign) or len(statement.targets) != 1 or not isinstance(statement.targets[0], ast.Name):
            raise ValueError('parameter files accept data assignments only')
        name = statement.targets[0].id
        if name not in _ALLOWED:
            raise ValueError('unknown parameter-file control: ' + name)
        if name in names:
            raise ValueError('duplicate parameter-file control: ' + name)
        names.add(name)
        if any(not isinstance(node, _EXPR) for node in ast.walk(statement.value)):
            raise ValueError('data-only expression required for ' + name)
    missing = _REQUIRED - names
    if missing:
        raise ValueError('missing parameter-file controls: ' + ', '.join(sorted(missing)))
    namespace = {'__builtins__': {}}
    exec(compile(tree, str(source), 'exec'), namespace)
    for name in ('PARAMETERS', 'EXTERNAL_INPUTS', 'NATIVE_OVERRIDES'):
        if not isinstance(namespace[name], dict) or any(not isinstance(k, str) for k in namespace[name]):
            raise ValueError(name + ' must be a dictionary with string keys')
    if set(namespace['PARAMETERS']) != set(DEFAULT_PARAMETERS):
        raise ValueError('PARAMETERS must contain exactly the ten supported coordinates')
    for name in ('PRICE_GUESS', 'BUDGET_SECONDS'):
        value = namespace[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(name + ' must be positive and finite')
    if namespace['CLOSURE'] not in ('fixed_h0', 'population_one'):
        raise ValueError('unsupported CLOSURE')
    provenance = namespace.get('PROVENANCE', {})
    if not isinstance(provenance, dict):
        raise ValueError('PROVENANCE must be a dictionary')
    # Native known-field, shape, physical-domain and conflicting-input checks.
    load_inputs(namespace['PARAMETERS'], namespace['EXTERNAL_INPUTS'], namespace['NATIVE_OVERRIDES'])
    data = {name.lower(): copy.deepcopy(namespace[name]) for name in _REQUIRED}
    data.update(provenance=copy.deepcopy(provenance), config_source=str(source),
                config_sha256=hashlib.sha256(raw).hexdigest(), config_text=raw.decode('utf-8'))
    return data


def output_root_for(path='best_params.py'):
    """Isolate personal files and reject competing files with the same stem.

    Read-only routing: ownership is registered by claim_output_root only when
    a run is requested. Existing roots must belong to this exact source path.
    """
    source = resolve_parameter_file(path)
    if source == BEST_FILE.resolve():
        return PROJECT_ROOT / 'output/model/local_solution'
    stem = source.stem
    if stem == 'best_params':
        # Separate exported winners with identical filenames from distinct runs.
        stem += '__' + source.parent.name + '__' + hashlib.sha256(str(source).encode()).hexdigest()[:8]
    root = PROJECT_ROOT / 'output/model/experiments' / stem
    if root.resolve() == (PROJECT_ROOT / 'output/model/local_solution').resolve():
        raise ValueError('experimental output collides with production')
    marker = root / '.parameter_file_source.json'
    if root.exists():
        if not marker.is_file() or json.loads(marker.read_text()).get('source') != str(source):
            raise ValueError('experimental output name is already owned: ' + str(root))
    return root


def claim_output_root(path='best_params.py'):
    """Register an experiment source immediately before its actual run."""
    source = resolve_parameter_file(path)
    root = output_root_for(source)
    if source != BEST_FILE.resolve() and not root.exists():
        root.mkdir(parents=True, exist_ok=False)
        (root / '.parameter_file_source.json').write_text(json.dumps({'source': str(source)}, indent=2) + '\n')
    return root


def describe_saved_case(case_directory, parameter_file):
    """Identify the cached case and disclose edits since its successful solve."""
    import hashlib
    import json
    case = Path(case_directory).resolve()
    note = f"Saved case: {case}"
    contract = case / "input_contract.json"
    if contract.is_file():
        metadata = json.loads(contract.read_text()).get("parameter_file") or {}
        if metadata.get("config_sha256"):
            source = resolve_parameter_file(parameter_file)
            current = hashlib.sha256(source.read_bytes()).hexdigest()
            if current != metadata["config_sha256"]:
                note += "\nParameter file has changed since this solve; showing the last successful saved run."
    return note
