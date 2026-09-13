#!/usr/bin/env python3
"""Watch the shared September deck and publish successful latexmk builds."""

import argparse
import os
from pathlib import Path
import shlex
import shutil
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "tmp/september_slides_autobuild"
NAME = "september_14_presentation"


def publish():
    source = BUILD / f"{NAME}.pdf"
    with source.open("rb") as stream:
        if stream.read(5) != b"%PDF-":
            raise RuntimeError("Build output is not a PDF; keeping published copies")
    for folder in (ROOT / "latex", ROOT / "output/pdf"):
        destination = folder / source.name
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=folder, suffix=".pdf", delete=False) as stream:
                temporary = Path(stream.name)
            shutil.copyfile(source, temporary)
            temporary.chmod(0o644)
            os.replace(temporary, destination)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    print("Published both September slide PDF copies.", flush=True)


def watch():
    BUILD.mkdir(parents=True, exist_ok=True)
    latexmk = shutil.which("latexmk") or "/Library/TeX/texbin/latexmk"
    command = shlex.join([sys.executable, str(Path(__file__).resolve()), "--publish"])
    perl_string = command.replace("\\", "\\\\").replace("'", "\\'")
    os.chdir(ROOT / "latex")
    os.execv(latexmk, [
        latexmk, "-norc", "-pdf", "-pvc", "-view=none", "-g",
        "-interaction=nonstopmode", "-halt-on-error", "-file-line-error",
        f"-outdir={BUILD}", "-e", f"$success_cmd = '{perl_string}';",
        f"{NAME}.tex",
    ])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--publish", action="store_true", help="latexmk success hook")
    args = parser.parse_args()
    publish() if args.publish else watch()
