#!/usr/bin/env bash
set -euo pipefail

ENTRYPOINT_BIN=${ENTRYPOINT_BIN:-/entrypoint.sh}
DATA_DIR=${DATA_DIR:-/opt/test_cifs}
PY_BIN_OVERRIDE=${PY_BIN_OVERRIDE:-}

if [ ! -x "$ENTRYPOINT_BIN" ]; then
  echo "Entrypoint script not found or not executable: $ENTRYPOINT_BIN" >&2
  exit 1
fi

if [ ! -d "$DATA_DIR" ]; then
  echo "Test data directory not found: $DATA_DIR" >&2
  exit 1
fi

if [ -n "$PY_BIN_OVERRIDE" ]; then
  if ! command -v "$PY_BIN_OVERRIDE" >/dev/null 2>&1; then
    echo "Requested Python override not available: $PY_BIN_OVERRIDE" >&2
    exit 1
  fi
  export PY_BIN="$PY_BIN_OVERRIDE"
fi

cleanup() {
  if [ -n "${TMP_WORKDIR:-}" ] && [ -d "$TMP_WORKDIR" ]; then
    rm -rf "$TMP_WORKDIR"
  fi
}
trap cleanup EXIT

TMP_WORKDIR=$(mktemp -d /tmp/entrypoint_cli_test.XXXXXX)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

pushd "$TMP_WORKDIR" >/dev/null

# Smoke test: combined help should render
"$ENTRYPOINT_BIN" --help >/dev/null

echo "Running sanity-only smoke test via entrypoint" >&2
SINGLE_OUTPUT="$TMP_WORKDIR/single-output"
"$ENTRYPOINT_BIN" \
  --run-sanity-only \
  --sanity-total-timeout 90 \
  --output-dir "$SINGLE_OUTPUT" \
  --n-jobs 1 \
  -i "$DATA_DIR/qmof-da78929.cif"

for required in \
  "$SINGLE_OUTPUT/results/qmof-da78929/qmof-da78929.json" \
  "$SINGLE_OUTPUT/results/qmof-da78929/qmof-da78929__base.json" \
  "$SINGLE_OUTPUT/results/qmof-da78929/qmof-da78929.log"; do
  if [ ! -f "$required" ]; then
    echo "Missing expected output artifact: $required" >&2
    exit 1
  fi
done

echo "Running oxichecker-only smoke test via entrypoint" >&2
"$ENTRYPOINT_BIN" \
  --run-oxichecker-only \
  --n-jobs 1 \
  -i "$SINGLE_OUTPUT"

if [ ! -f "$TMP_WORKDIR/oxichecker_results.csv" ]; then
  echo "OxiChecker run did not produce oxichecker_results.csv" >&2
  exit 1
fi

echo "Running postprocess-only smoke test via entrypoint" >&2
"$ENTRYPOINT_BIN" \
  --run-postprocess-only \
  --n-jobs 1 \
  --output-dir "$SINGLE_OUTPUT"

if [ ! -f "$TMP_WORKDIR/sanity_results.csv" ]; then
  echo "Postprocessing run did not produce sanity_results.csv" >&2
  exit 1
fi

popd >/dev/null

echo "Entrypoint CLI smoke tests completed successfully" >&2
