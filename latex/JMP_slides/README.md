# JMP Slides

The continuing deck is `JMP_slides.tex`, copied unchanged from the September 14 presentation. The reader PDF remains `../../output/pdf/JMP_slides.pdf`.

Compile from the parent `latex/` directory so existing figure paths remain valid:

```sh
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=../tmp/pdfs/jmp_setup/slides JMP_slides/JMP_slides.tex
```

Run twice after restoring the missing legacy asset `../../Outputs/Graphs/own_f_c_y_all.png`. The existing PDF is the unchanged September copy; the latest rebuild was blocked by that asset.
