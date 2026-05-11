"""Assemble docs/tour.html from the template, inlining JSON, code excerpts, and widget JS.

Usage:
    python tools/assemble_tour.py
"""

from __future__ import annotations

import html
import json
import re
from pathlib import Path
from typing import Iterable


HERE = Path(__file__).resolve().parent
PORT_ROOT = HERE.parent
TEMPLATE = HERE / "tour_template.html"
WIDGETS = HERE / "tour_widgets.js"
DATA = PORT_ROOT / "docs" / "tour_data.json"
OUTPUT = PORT_ROOT / "docs" / "tour.html"

KERNELS_PY = PORT_ROOT / "dt_cp_model" / "kernels.py"
SOLVER_PY = PORT_ROOT / "dt_cp_model" / "solver.py"

DATA_PLACEHOLDER = re.compile(
    r'<script id="tour-data" type="application/json">.*?</script>',
    re.DOTALL,
)


# Each panel: (anchor, title, math LaTeX, code extraction spec).
# Code spec is one of:
#   ("def", file, function_name, max_lines)
#   ("range", file, start_pattern, end_pattern, max_lines)
DUAL_PANELS = [
    {
        "anchor": "dual-utility",
        "n": "3.1",
        "title": "Period utility",
        "math": r"""u(c, h, n) = \frac{\big[(c - \bar c_0 - \bar c_n n)^\alpha (h - \bar h_0 - \bar h_n n - \bar h_{\mathrm{jump}}\mathbf{1}[n\geq 1])^{1-\alpha}\big]^{1-\sigma} - 1}{1-\sigma}""",
        "code_spec": ("def", KERNELS_PY, "eval_renter_scalar", 36),
    },
    {
        "anchor": "dual-renter-bellman",
        "n": "3.2",
        "title": "Renter Bellman block",
        "math": r"""V_d^{\mathrm{rent}}(b, i, n, \mathrm{cs}, j) = \max_{b', h \leq h_R^{\max}} u(c, h, n) + \beta\, \mathbb{E}_{\mathrm{cs}'}\big[V(b', 0, i, n, \mathrm{cs}', j+1)\big]""",
        "code_spec": ("def", KERNELS_PY, "full_renter_block_kernel", 36),
    },
    {
        "anchor": "dual-owner-bellman",
        "n": "3.3",
        "title": "Owner Bellman block",
        "math": r"""V_d^{\mathrm{own}}(b, t, i, n, \mathrm{cs}, j) = \max_{b'} u\big(c, \chi H_{\mathrm{own}}[t-1], n\big) + \beta\, \mathbb{E}_{\mathrm{cs}'}\big[V(b', t, i, n, \mathrm{cs}', j+1)\big]""",
        "code_spec": ("def", KERNELS_PY, "full_owner_block_kernel", 36),
    },
    {
        "anchor": "dual-tenure-max",
        "n": "3.4",
        "title": "Tenure max + feasibility",
        "math": r"""V_H(b, t_o, i, n, \mathrm{cs}) = \max_{t_n \in \mathcal{F}(b, t_o, i)} V_d(b^*(b, t_o, t_n, i), t_n, i, n, \mathrm{cs})""",
        "code_spec": ("def", KERNELS_PY, "tenure_choice_kernel", 36),
    },
    {
        "anchor": "dual-loc-logsumexp",
        "n": "3.5",
        "title": "Location logsumexp",
        "math": r"""V_I(b, t_o, i, \cdot) = \kappa_{\mathrm{loc}} \log \sum_{i'} \exp\Big(\tfrac{V_H(b_{i'}, t_o, i', \cdot) + E_{i'} - mc \cdot \mathbf{1}[i' \neq i]}{\kappa_{\mathrm{loc}}}\Big)""",
        "code_spec": ("def", KERNELS_PY, "location_logit_kernel", 40),
    },
    {
        "anchor": "dual-fert-logit",
        "n": "3.6",
        "title": "Fertility logit",
        "math": r"""\Pr(n' \mid b, i, j) = \frac{\exp(V_I(b, 0, i, n', \mathrm{cs}'(n'), j)/\kappa_{\mathrm{fert}})}{\sum_{n''}\exp(V_I(b, 0, i, n'', \mathrm{cs}'(n''), j)/\kappa_{\mathrm{fert}})}""",
        "code_spec": ("range", SOLVER_PY, r"^\s+if in_fert:$", r"^\s+return V,", 16),
    },
    {
        "anchor": "dual-forward",
        "n": "3.7",
        "title": "Forward distribution kernel",
        "math": r"""g_{j+1}(b', t_n, i', n, \mathrm{cs}') = \sum_{(b, t_o, i, \mathrm{cs})} g_j(b, t_o, i, n, \mathrm{cs}) \cdot \Pi^{\mathrm{loc}}_{i \to i'}(b, t_o, n, \mathrm{cs}) \cdot \mathbf{1}[t_n = t^*(b, t_o, i', n, \mathrm{cs})] \cdot W_{b'}(b, \cdot) \cdot \Pi^{\mathrm{cs}}_{\mathrm{cs} \to \mathrm{cs}'}""",
        "code_spec": ("def", KERNELS_PY, "forward_distribution_fast_kernel", 50),
    },
    {
        "anchor": "dual-ge",
        "n": "3.8",
        "title": "GE price update",
        "math": r"""p^{(\text{new})}_i = \bar r_i \left(\frac{H^d_i}{H_{0,i}}\right)^{1/\xi_i} \!\!\!/\, \mathrm{ucr}, \qquad p^{(it+1)}_i = p^{(it)}_i + \lambda_i^p \cdot (p^{(\text{new})}_i - p^{(it)}_i)""",
        "code_spec": ("range", SOLVER_PY, r"p_target = np.zeros\(P.I\)", r"err_p = float", 24),
    },
]


def _read_text(p: Path) -> str:
    return p.read_text()


def _extract_def(src_text: str, fn_name: str, max_lines: int = 40) -> tuple[str, int, int]:
    """Extract a def block starting at `def fn_name(` and ending at the next def or eof."""
    lines = src_text.splitlines()
    pat_start = re.compile(r"^(\s*)def\s+" + re.escape(fn_name) + r"\b")
    start = None
    for i, line in enumerate(lines):
        if pat_start.match(line):
            start = i
            break
    if start is None:
        raise RuntimeError(f"function {fn_name} not found")
    # walk forward until next top-level def or class, or max_lines
    end = min(start + max_lines, len(lines))
    for i in range(start + 1, end):
        if re.match(r"^def\s|^class\s|^@\w+", lines[i]):
            end = i
            break
    return "\n".join(lines[start:end]).rstrip(), start + 1, end


def _extract_range(src_text: str, start_pat: str, end_pat: str, max_lines: int = 30) -> tuple[str, int, int]:
    lines = src_text.splitlines()
    sp = re.compile(start_pat)
    ep = re.compile(end_pat)
    start = None
    for i, line in enumerate(lines):
        if sp.search(line):
            start = i
            break
    if start is None:
        raise RuntimeError(f"start pattern {start_pat!r} not found")
    end = min(start + max_lines, len(lines))
    for i in range(start + 1, end):
        if ep.search(lines[i]):
            end = i
            break
    return "\n".join(lines[start:end]).rstrip(), start + 1, end


def extract_code(spec) -> tuple[str, str, int, int]:
    """Returns (code, file_label, start_line, end_line)."""
    kind = spec[0]
    if kind == "def":
        _, path, fn, max_lines = spec
        text = _read_text(path)
        code, s, e = _extract_def(text, fn, max_lines=max_lines)
        return code, path.name + f":{s}-{e}", s, e
    elif kind == "range":
        _, path, start_pat, end_pat, max_lines = spec
        text = _read_text(path)
        code, s, e = _extract_range(text, start_pat, end_pat, max_lines=max_lines)
        return code, path.name + f":{s}-{e}", s, e
    else:
        raise ValueError(f"unknown spec kind {kind}")


def render_dual_panels() -> str:
    """Build the HTML for the eight chapter-3 dual-view panels."""
    parts: list[str] = []
    for p in DUAL_PANELS:
        try:
            code, label, _, _ = extract_code(p["code_spec"])
        except Exception as exc:
            code = f"# extraction failed: {exc}"
            label = "extraction failed"
        code_escaped = html.escape(code)
        parts.append(f"""
  <h3 id="{p['anchor']}">{p['n']} {html.escape(p['title'])}</h3>
  <div class="dual">
    <div class="dual-math">
      <span class="label">math</span>
      $$\\displaystyle {p['math']}$$
    </div>
    <div class="dual-code">
      <span class="label">code · {html.escape(label)}</span>
      <pre><code class="language-python">{code_escaped}</code></pre>
    </div>
  </div>""")
    return "\n".join(parts)


def assemble() -> Path:
    template = _read_text(TEMPLATE)
    data_text = _read_text(DATA)
    widgets_js = _read_text(WIDGETS)

    # Build chapter 3 dual-view panels
    dual_html = render_dual_panels()
    template = template.replace("<!-- DUAL_VIEW_PANELS_GO_HERE -->", dual_html)

    # Inject widgets JS as an inline script just before </body>
    widgets_block = (
        '<script id="tour-widgets">\n' +
        widgets_js +
        '\n</script>'
    )
    template = template.replace("<!-- WIDGETS_GO_HERE -->", widgets_block)

    # Inject the JSON data
    safe = data_text.replace("</", r"<\/")
    if not DATA_PLACEHOLDER.search(template):
        raise RuntimeError("Data placeholder not found in template")
    template = DATA_PLACEHOLDER.sub(
        '<script id="tour-data" type="application/json">' + safe + '</script>',
        template,
        count=1,
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(template)
    size_kb = OUTPUT.stat().st_size / 1024
    print(f"Wrote {OUTPUT} ({size_kb:.1f} KB)")
    return OUTPUT


if __name__ == "__main__":
    assemble()
