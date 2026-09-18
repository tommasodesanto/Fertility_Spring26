"""Minimal YAML-subset reader for sandbox/specs/*.yaml.

PyYAML is not installed in code/model/.venv, so this parses exactly the
subset the specs use: a flat mapping whose values are scalars
(int / float / true / false / null / bare or quoted string) or flow-style
lists of scalars, including nested lists (e.g. the Pi_z matrices and the
z_grid / permanent-level arrays in earnings_e6b_*.yaml and
earnings_bgm_*.yaml, and child_earnings_penalty):

    # comment
    key: value
    overrides:
      nested_key: value
      some_list: [0, 0.2, 0.2, 0.2]
      some_matrix: [[1.0, 0.0], [0.0, 1.0]]
"""
from __future__ import annotations

from pathlib import Path
from typing import Any


def _split_top_level(body: str) -> list[str]:
    """Split a flow-list body on top-level commas (bracket- and quote-aware)."""
    parts: list[str] = []
    depth = 0
    quote: str | None = None
    current: list[str] = []
    for char in body:
        if quote is not None:
            current.append(char)
            if char == quote:
                quote = None
        elif char in ("'", '"'):
            quote = char
            current.append(char)
        elif char == "[":
            depth += 1
            current.append(char)
        elif char == "]":
            depth -= 1
            if depth < 0:
                raise ValueError(f"Unbalanced brackets in spec value: [{body}]")
            current.append(char)
        elif char == "," and depth == 0:
            parts.append("".join(current))
            current = []
        else:
            current.append(char)
    if quote is not None:
        raise ValueError(f"Unterminated quote in spec value: [{body}]")
    if depth != 0:
        raise ValueError(f"Unbalanced brackets in spec value: [{body}]")
    parts.append("".join(current))
    return parts


def _coerce_scalar(token: str) -> Any:
    token = token.strip()
    if token in ("", "~", "null", "None"):
        return None
    if token.lower() == "true":
        return True
    if token.lower() == "false":
        return False
    if len(token) >= 2 and token.startswith(("'", '"')) and token.endswith(token[0]):
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


def _coerce(token: str) -> Any:
    """Parse one scalar or flow-list value (nested lists supported)."""
    token = token.strip()
    if token.startswith("["):
        if not token.endswith("]"):
            raise ValueError(f"Unbalanced brackets in spec value: {token!r}")
        body = token[1:-1].strip()
        if body == "":
            return []
        return [_coerce(part) for part in _split_top_level(token[1:-1])]
    return _coerce_scalar(token)


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
