# Latexmk configuration for this project.
# On this TeX Live setup, pdflatex does not support a separate aux directory
# through latexmk without configuration warnings and inconsistent outputs.
# Keep the build honest: compile in-place, with the PDF in this folder.

$out_dir = '.';
