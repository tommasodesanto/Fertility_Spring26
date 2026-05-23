# Couillard: Build, Baby, Build

Local working copies, not tracked by git:

- `BuildBabyBuild_Couillard_2025.pdf`
- `BuildBabyBuild_Couillard_2025.txt`

Source URL: <https://br.ti.org/pdfs/BuildBabyBuild.pdf>

The text file was extracted locally from the PDF with Ghostscript:

```bash
gs -q -dNOPAUSE -dBATCH -sDEVICE=txtwrite \
  -sOutputFile=docs/literature/couillard/BuildBabyBuild_Couillard_2025.txt \
  docs/literature/couillard/BuildBabyBuild_Couillard_2025.pdf
```

This folder keeps the paper available for local reading without relying on web
search. The PDF and extracted full text are ignored to avoid committing a
third-party paper into the repository.
