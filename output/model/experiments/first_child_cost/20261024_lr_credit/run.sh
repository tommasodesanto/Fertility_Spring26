#!/bin/bash
export MF_ROOT=/private/tmp/claude-501/-Users-tommasodesanto-Desktop-Projects-Fertility-Fertility-Spring26/8632614c-e3af-47a4-905e-a3434441141e/scratchpad/mf
export MF_UPLOADS=/Users/tommasodesanto/Desktop/Projects/Fertility MF_POINT=newacct
PY=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python
cd "$(dirname "$0")"; mkdir -p logs
for c in "$@"; do set -- $c; $PY lr_credit_cell.py $1 $2 $3 > "logs/${1}_${2}_${3}.log" 2>&1 & done; wait
