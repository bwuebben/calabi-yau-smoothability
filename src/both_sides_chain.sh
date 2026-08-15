#!/bin/bash
# Reproduce the full both-sides scan over the complete Kreuzer--Skarke
# database (paper 3, Sec. 5).  The driver pins and verifies every input,
# validates an existing result before resuming past it, writes new results
# atomically, and records a complete transcript.
#
# Per vertex-count file NN = 05..33:
#   1. download the parquet from an immutable Hugging Face revision if absent;
#   2. verify byte size and SHA-256 against the committed input manifest;
#   3. if a result exists, verify its schema, row count, aggregates, and every
#      positive hit with the exact reference toolkit before skipping it;
#   4. otherwise cross-check the engines on 40 rows and run the full fast scan.
# The separately stored 36-vertex polytope is input-verified and then checked
# by missing_polytope.py.  Finally both_sides_census.py re-derives all paper-3
# claims from the result artifacts.
#
# This protects provenance and positive results.  It does not constitute the
# independent re-scan of every negative row discussed in paper 3.
#
# Usage: ./src/both_sides_chain.sh [procs]  (default: all cores minus two)
set -euo pipefail
cd "$(dirname "$0")/.."

PY=${PY:-./venv/bin/python}
REVISION="60c0e119a03608418df538191f65da3f43b5b819"
BASE="https://huggingface.co/datasets/calabi-yau-data/polytopes-4d/resolve/$REVISION"
MANIFEST="manifests/ks_polytopes_4d_sha256.tsv"
VERIFIER="src/verify_both_sides_artifact.py"
RUN_LOG=${RUN_LOG:-output/both_sides_chain_transcript.txt}

if [ ! -x "$PY" ]; then
    echo "Python interpreter is not executable: $PY" >&2
    exit 1
fi

CPU_COUNT=$(sysctl -n hw.ncpu 2>/dev/null || nproc)
if [ "$CPU_COUNT" -gt 2 ]; then DEFAULT_PROCS=$((CPU_COUNT - 2)); else DEFAULT_PROCS=1; fi
PROCS=${1:-$DEFAULT_PROCS}
if [ "$PROCS" -lt 1 ]; then
    echo "process count must be positive" >&2
    exit 1
fi

mkdir -p data/ks/.partial output
exec > >(tee -a "$RUN_LOG") 2>&1

echo "=== both-sides chain start: $(date -u '+%Y-%m-%dT%H:%M:%SZ') ==="
echo "dataset revision: $REVISION"
echo "source commit: $(git rev-parse HEAD 2>/dev/null || echo unavailable)"
if ! git diff --quiet -- . 2>/dev/null || ! git diff --cached --quiet -- . 2>/dev/null; then
    echo "source worktree: dirty"
else
    echo "source worktree: clean"
fi
echo "python: $PY"
echo "processes: $PROCS"

ensure_input() {
    local filename=$1
    local target="data/ks/$filename"
    local partial="data/ks/.partial/$filename"
    if [ ! -e "$target" ]; then
        echo "[$(date '+%H:%M:%S')] downloading $filename"
        curl -fL --retry 5 --retry-delay 2 --continue-at - \
             -o "$partial" "$BASE/$filename"
        "$PY" "$VERIFIER" --input "$partial" --manifest "$MANIFEST"
        mv "$partial" "$target"
    else
        "$PY" "$VERIFIER" --input "$target" --manifest "$MANIFEST"
    fi
}

for NN in 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 \
          26 27 28 29 30 31 32 33; do
    F="polytopes-4d-${NN}-vertices.parquet"
    J="output/both_sides_v${NN}.json"
    # Preserve the two combined low-vertex artifacts from the original run.
    case "$NN" in
        05|06|07)
            [ -e output/both_sides_v0507_fast.json ] && \
                J="output/both_sides_v0507_fast.json"
            ;;
        08|09)
            [ -e output/both_sides_v0809_fast.json ] && \
                J="output/both_sides_v0809_fast.json"
            ;;
    esac

    ensure_input "$F"
    if [ -e "$J" ]; then
        echo "[$(date '+%H:%M:%S')] validating existing result for v$NN"
        "$PY" "$VERIFIER" --input "data/ks/$F" --manifest "$MANIFEST" \
              --result "$J"
        echo "[skip] v$NN has a validated complete result"
        continue
    fi

    echo "[$(date '+%H:%M:%S')] cross-checking engines on v$NN"
    "$PY" -c "import sys; sys.path.insert(0, 'src'); \
from ks_sweep import selftest_engines; \
print('selftest rows ok:', selftest_engines('data/ks/$F', 40))"
    echo "[$(date '+%H:%M:%S')] scanning v$NN"
    "$PY" src/both_sides_fast.py "data/ks/$F" --procs "$PROCS" \
          --dataset-revision "$REVISION" --json "$J"
    "$PY" "$VERIFIER" --input "data/ks/$F" --manifest "$MANIFEST" \
          --result "$J"
done

F36="polytopes-4d-36-vertices.parquet"
ensure_input "$F36"
"$PY" src/missing_polytope.py
"$PY" src/both_sides_census.py

echo "=== both-sides chain complete: $(date -u '+%Y-%m-%dT%H:%M:%SZ') ==="
