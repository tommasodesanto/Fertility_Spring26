"""Minimal flat-mapping YAML reader for sandbox/specs/*.yaml.

PyYAML is not installed in code/model/.venv, and every spec here is a flat
mapping of scalar overrides (no lists, no nesting), so a full YAML parser is
not worth adding as a dependency. This reads exactly that subset:

    # comment
    key: value          # int / float / true / false / null / bare string
    overrides:
      nested_key: value
"""
from __future__ import annotations

from pathlib import Path
from typing import Any


def _coerce(token: str) -> Any:
    token = token.strip()
    if token in ("", "~", "null", "None"):
        return None
    if token.lower() == "true":
        return True
    if token.lower() == "false":
        return False
    if token.startswith(('"', "'")) and token.endswith(('"', "'")):
        return token[1:-1]
    try:
        return int(token)
    except ValueError:
        pass
    try:
        return float(token)
    except ValueError:
        pass
    return token


def load_spec(path: Path) -> dict[str, Any]:
    spec: dict[str, Any] = {"overrides": {}}
    current_indent_is_overrides = False
    if not path.exists():
        raise FileNotFoundError(f"Spec not found: {path}")
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        indent = len(raw_line) - len(raw_line.lstrip(" "))
        key, _, value = line.strip().partition(":")
        key = key.strip()
        value = value.strip()
        if indent == 0 and key == "overrides" and not value:
            current_indent_is_overrides = True
            continue
        if indent == 0:
            current_indent_is_overrides = False
            spec[key] = _coerce(value)
        elif current_indent_is_overrides:
            spec["overrides"][key] = _coerce(value)
        else:
            raise ValueError(f"Unsupported spec indentation in {path}: {raw_line!r}")
    return spec
