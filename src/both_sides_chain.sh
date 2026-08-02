#!/bin/bash
# Reproduce the full both-sides scan over the complete Kreuzer-Skarke
# database (paper 3, Sec. 5) -- the canonical driver for the chain whose
# run is logged in output/both_sides_chain.log.
#
# Per vertex-count file NN = 05..33:
#   1. download data/ks/polytopes-4d-NN-vertices.parquet from the HF mirror
#      (https://huggingface.co/datasets/calabi-yau-data/polytopes-4d) if absent;
#   2. selftest: compare the fast engine against the pure-python reference
#      on 40 random rows of THAT file (ks_sweep.selftest_engines);
#   3. scan with src/both_sides_fast.py --> output/both_sides_vNN.json.
# Files already having a JSON are skipped, so the chain is resumable.
# The 36-vertex polytope (absent from the per-vertex mirror) is handled
# separately by src/missing_polytope.py; src/both_sides_census.py then
# re-derives and asserts every number of paper 3 from the JSONs.
#
# Usage:  ./src/both_sides_chain.sh [procs]     (default: all cores - 2)
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-python3}
PROCS=${1:-$(($(sysctl -n hw.ncpu 2>/dev/null || nproc) - 2))}
BASE="https://huggingface.co/datasets/calabi-yau-data/polytopes-4d/resolve/main"
mkdir -p data/ks output

for NN in 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 \
          26 27 28 29 30 31 32 33; do
    F="polytopes-4d-${NN}-vertices.parquet"
    J="output/both_sides_v${NN}.json"
    # legacy names from the original run cover v05-v09 in combined files
    case "$NN" in
        05|06|07) [ -s output/both_sides_v0507_fast.json ] && J="" ;;
        08|09)    [ -s output/both_sides_v0809_fast.json ] && J="" ;;
    esac
    if [ -z "$J" ] || [ -s "$J" ]; then echo "[skip] v$NN done"; continue; fi
    if [ ! -s "data/ks/$F" ]; then
        echo "[$(date +%T)] downloading v$NN ..."
        curl -fsSL -o "data/ks/$F" "$BASE/$F"
    fi
    echo "[$(date +%T)] selftest v$NN"
    $PY -c "import sys; sys.path.insert(0,'src'); \
from ks_sweep import selftest_engines; \
print('selftest rows ok:', selftest_engines('data/ks/$F', 40))"
    echo "[$(date +%T)] scanning v$NN"
    $PY src/both_sides_fast.py "data/ks/$F" --procs "$PROCS" --json "$J"
done
echo "BOTH-SIDES-CHAIN-ALL-DONE"
