#!/usr/bin/env python3
"""Render the frozen audit document; never import or solve the economic model.

Requires pandoc, XeLaTeX, and the macOS Palatino/Helvetica/Menlo fonts.
Run from any directory with Python 3. Intermediate files stay in tmp/pdfs/.
"""
from pathlib import Path
import argparse
import re
import subprocess

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "docs/model/e5f_independent_quantitative_audit.md"
TEMP = ROOT / "tmp/pdfs/e5f_independent_audit"
OUTPUT = ROOT / "output/pdf/e5f_independent_quantitative_audit.pdf"

HEADER = r"""
\usepackage{hyperref,xurl,etoolbox,fancyhdr,titlesec,needspace}
\definecolor{auditblue}{HTML}{173E55}
\definecolor{auditgray}{HTML}{53616B}
\AtBeginDocument{\hypersetup{linkcolor=auditblue,urlcolor=auditblue}}
\setlength{\emergencystretch}{4em}
\setlength{\parskip}{5.5pt}
\setlength{\parindent}{0pt}
\setlength{\tabcolsep}{4pt}
\renewcommand{\arraystretch}{1.15}
\AtBeginEnvironment{longtable}{\footnotesize\setlength{\parskip}{3pt}}
\titleformat{\section}{\Large\sffamily\bfseries\color{auditblue}}{}{0pt}{}
\titleformat{\subsection}{\large\sffamily\bfseries\color{auditblue}}{}{0pt}{}
\titleformat{\subsubsection}{\normalsize\sffamily\bfseries\color{auditblue}}{}{0pt}{}
\titlespacing*{\section}{0pt}{12pt}{10pt}
\titlespacing*{\subsection}{0pt}{12pt}{8pt}
\titlespacing*{\subsubsection}{0pt}{10pt}{6pt}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\footnotesize\sffamily\color{auditgray}FERTILITY / INDEPENDENT QUANTITATIVE AUDIT}
\fancyhead[R]{\footnotesize\sffamily\color{auditgray}4 SEPTEMBER 2026}
\fancyfoot[L]{\footnotesize\sffamily\color{auditgray}Discussion copy | production calibration unchanged}
\fancyfoot[R]{\footnotesize\sffamily\thepage}
\renewcommand{\headrulewidth}{0.3pt}
\setlength{\headheight}{14pt}
\fancypagestyle{plain}{\fancyhf{}\fancyfoot[R]{\footnotesize\sffamily\thepage}\renewcommand{\headrulewidth}{0pt}}
\setcounter{tocdepth}{2}
\urlstyle{same}
\widowpenalty=10000
\clubpenalty=10000
"""

FILTER = r"""
function Header(h)
  local title = pandoc.utils.stringify(h.content)
  if h.level == 2 then
    return {pandoc.RawBlock('latex', '\\Needspace{16\\baselineskip}\\phantomsection'), h}
  end
  return {pandoc.RawBlock('latex', '\\Needspace{8\\baselineskip}\\phantomsection'), h}
end
function Table(t)
  local n = #t.colspecs
  local widths
  local label = pandoc.utils.stringify(t)
  if n == 2 then widths = {0.29,0.71}
  elseif n == 3 then
    if label:match('Authoritative source') then widths = {0.19,0.43,0.38}
    elseif label:match('Impact births per household') then widths = {0.50,0.25,0.25}
    else widths = {0.24,0.37,0.39} end
  elseif n == 4 then
    if label:match('Diagnostic measurement') then widths = {0.44,0.14,0.15,0.27}
    elseif label:match('Impact births per household') then widths = {0.39,0.20,0.20,0.21}
    elseif label:match('Baseline quantity') then widths = {0.40,0.18,0.18,0.24}
    else widths = {0.23,0.24,0.26,0.27} end
  elseif n == 5 then
    if label:match('Panel 003') then widths = {0.29,0.14,0.14,0.26,0.17}
    elseif label:match('Wealth nodes') then widths = {0.32,0.14,0.18,0.18,0.18}
    elseif label:match('Rental limit') then widths = {0.32,0.14,0.18,0.18,0.18}
    else widths = {0.15,0.20,0.25,0.24,0.16} end
  elseif n == 6 then
    if label:match('Estimate') then widths = {0.20,0.10,0.12,0.15,0.20,0.23}
    else widths = {0.32,0.135,0.135,0.14,0.135,0.135} end
  elseif n == 8 then widths = {0.23,0.11,0.09,0.09,0.12,0.12,0.12,0.12}
  end
  if widths then
    for i=1,n do t.colspecs[i] = {t.colspecs[i][1], widths[i]} end
  end
  return t
end
function Code(c)
  if #c.text > 38 then
    return pandoc.RawInline('latex', '\\path{' .. c.text .. '}')
  end
end
"""


def prepare_markdown(text, source_index=True):
    """Preserve link evidence in a readable, portable source index."""
    # The installed TeX format and newer babel package are incompatible; the
    # document is English and needs no language-specific babel conventions.
    text = re.sub(r"^lang:.*\n", "", text, flags=re.MULTILINE)
    text = text.replace("→", " to ")
    text = text.replace("Births/household 2023", "Births per household 2023")
    text = text.replace("Rooms/household 2023", "Rooms per household 2023")
    if not source_index:
        return text
    sources = {}
    def link(match):
        label, target = match.groups()
        if target.startswith("#"):
            return match.group(0)
        if target not in sources:
            sources[target] = (len(sources) + 1, label)
        number = sources[target][0]
        return f"{label} [S{number:03d}](#source-{number:03d})"
    text = re.sub(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)", link, text)
    # Pandoc turns the protected manuscript marker into literal text as intended.
    text += "\n\n\\clearpage\n\n# Source index\n\n"
    text += ("References reproduce the paths and line numbers used during the audit. "
             "Local source files may subsequently change. Run folders and saved receipts "
             "identify the reviewed artifacts; web links identify the cited primary sources.\n\n")
    for target, (number, label) in sources.items():
        text += ("```{=latex}\n\\Needspace{75pt}\\phantomsection"
                 f"\\label{{source-{number:03d}}}\n```\n\n"
                 f"**S{number:03d}. {label}**\n\n")
        text += "```{=latex}\n\\begingroup\\footnotesize\\url{" + target + "}\\par\\endgroup\n```\n\n"
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--date", default="5 September 2026")
    parser.add_argument("--heading", default="FERTILITY / INDEPENDENT QUANTITATIVE AUDIT")
    parser.add_argument("--no-source-index", action="store_true")
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    temporary = TEMP if output == OUTPUT else ROOT / "tmp/pdfs" / output.stem
    temporary.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    header = HEADER.replace("4 SEPTEMBER 2026", args.date.upper()).replace(
        "FERTILITY / INDEPENDENT QUANTITATIVE AUDIT", args.heading)
    (temporary / "header.tex").write_text(header)
    (temporary / "layout.lua").write_text(FILTER)
    (temporary / "document.md").write_text(prepare_markdown(source.read_text(), not args.no_source_index))
    command = [
        "pandoc", str(temporary / "document.md"),
        "--from=markdown+tex_math_single_backslash", "--standalone",
        "--pdf-engine=xelatex", "--lua-filter=" + str(temporary / "layout.lua"),
        "--include-in-header=" + str(temporary / "header.tex"),
        "--resource-path=" + str(ROOT) + ":" + str(source.parent),
        "--extract-media=" + str(temporary / "media"),
        "-V", "documentclass=article", "-V", "fontsize=11pt",
        "-V", "papersize=a4", "-V", "geometry:margin=20mm",
        "-V", "mainfont=Palatino", "-V", "sansfont=Helvetica",
        "-V", "monofont=Menlo", "-V", "linestretch=1.06",
        "-V", "colorlinks=true", "--top-level-division=section",
        "-o", str(temporary / "document.tex"),
    ]
    subprocess.run(command, check=True)
    for _ in range(3):
        subprocess.run(
            ["xelatex", "-interaction=nonstopmode", "-halt-on-error", "document.tex"],
            cwd=temporary, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
    output.write_bytes((temporary / "document.pdf").read_bytes())
    print(output)


if __name__ == "__main__":
    main()
