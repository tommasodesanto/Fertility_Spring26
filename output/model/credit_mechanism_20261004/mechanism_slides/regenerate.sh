#!/usr/bin/env bash
# Rebuild mechanism figures and slide PDFs from saved diagnostics only.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "$script_dir/../../../.." && pwd -P)"
packet="$script_dir"
build_dir="$packet/build"
figure_dir="$repo_root/latex/JMP_slides/mechanisms/figures"
pdf_dir="$repo_root/output/pdf"
python="$repo_root/code/model/.venv/bin/python"

export NUMBA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export MPLBACKEND=Agg
mkdir -p "$build_dir" "$figure_dir" "$pdf_dir"
export NUMBA_CACHE_DIR="$build_dir/numba_cache"

"$python" "$packet/household/build_household_policies.py"
"$python" "$packet/floor/extract_floor.py"
"$python" "$packet/floor/build_main_figure.py"
"$python" "$packet/responses/build_mechanism_figures.py"

install -m 644 "$packet/household/household_policies.pdf" "$figure_dir/household_policies.pdf"
install -m 644 "$packet/floor/main_floor_state.pdf" "$figure_dir/parenthood_floor.pdf"
install -m 644 "$packet/responses/first_birth_price110_by_beginning_tenure.pdf" "$figure_dir/first_birth_price110_by_beginning_tenure.pdf"
install -m 644 "$packet/responses/credit_phi095_wait_success_gains_debt_states.pdf" "$figure_dir/credit_phi095_wait_success_gains_debt_states.pdf"
install -m 644 "$packet/responses/renter_room_cap_births_ownership.pdf" "$figure_dir/renter_room_cap_births_ownership.pdf"

compile_twice() {
  local source="$1"
  local jobname="$2"
  local pass log
  for pass in 1 2; do
    log="$build_dir/${jobname}.pass${pass}.stdout.log"
    if ! pdflatex -interaction=nonstopmode -halt-on-error -file-line-error \
      -jobname="$jobname" -output-directory="$build_dir" "$source" >"$log" 2>&1; then
      tail -n 60 "$log" >&2
      return 1
    fi
  done
  test -s "$build_dir/$jobname.pdf"
}

cd "$repo_root/latex"
compile_twice "JMP_slides/mechanisms/mechanism_excerpt.tex" "JMP_mechanisms"
compile_twice "JMP_slides/JMP_slides.tex" "JMP_slides"
install -m 644 "$build_dir/JMP_mechanisms.pdf" "$pdf_dir/JMP_mechanisms.pdf"
install -m 644 "$build_dir/JMP_slides.pdf" "$pdf_dir/JMP_slides.pdf"

printf 'Mechanism slides: %s\nFull slides: %s\n' "$pdf_dir/JMP_mechanisms.pdf" "$pdf_dir/JMP_slides.pdf"
