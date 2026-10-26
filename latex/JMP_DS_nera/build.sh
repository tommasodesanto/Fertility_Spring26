#!/bin/zsh
# Build the Nera draft: aux files in build/, final PDF copied next to the tex.
# Plain pdflatex/bibtex sequence: bibtex reports a duplicate entry in ../new_entries.bib (harmless), which stops latexmk.
cd "$(dirname "$0")"
mkdir -p build
for step in 1 bib 2 3; do
  if [[ $step == bib ]]; then bibtex build/JMP_DS_nera >/dev/null 2>&1
  else pdflatex -interaction=nonstopmode -output-directory=build JMP_DS_nera.tex >/dev/null 2>&1; fi
done
cp build/JMP_DS_nera.pdf JMP_DS_nera.pdf
grep -E "undefined|multiply defined|^! " build/JMP_DS_nera.log | sort -u
echo "pages: $(pdfinfo JMP_DS_nera.pdf | awk '/Pages/{print $2}')"
