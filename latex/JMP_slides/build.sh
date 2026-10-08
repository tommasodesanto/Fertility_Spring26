#!/bin/sh
# Rebuild the JMP slides (main deck and the separate new-results deck) from the production readout.
#   1. Regenerate calibration tables (output/model/production/, via read_results.py) and policy numbers (result folder named in build_prod_tables.py).
#   2. Compile from latex/ (paths in the .tex are relative to it); auxiliary files go to tmp/, the PDF to output/pdf/JMP_slides.pdf.
# Usage: sh latex/JMP_slides/build.sh   (from anywhere)
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
cd "$ROOT" && python3 latex/JMP_slides/build_prod_tables.py > /dev/null
BUILD="$ROOT/tmp/jmp_slides_build"
mkdir -p "$BUILD"
for DECK in JMP_slides JMP_slides_new_results; do
  cd "$ROOT/latex" && latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir="$BUILD" "JMP_slides/$DECK.tex" > "$BUILD/latexmk_$DECK.txt" 2>&1 \
    || { echo "LaTeX failed; see $BUILD/$DECK.log"; exit 1; }
  cp "$BUILD/$DECK.pdf" "$ROOT/output/pdf/$DECK.pdf"
  echo "built $ROOT/output/pdf/$DECK.pdf"
done
