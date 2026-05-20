#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   run_cli_smoke.sh [step]
#
# When called without an argument (or with `all`) the script runs the full
# sequence in one process. When called with a step name it executes that step
# only and exits — the Dockerfile relies on this mode so each step shows up as
# its own line in the BuildKit build progress output.

ENTRYPOINT_BIN=${ENTRYPOINT_BIN:-/entrypoint.sh}
DATA_DIR=${DATA_DIR:-/opt/test_cifs}
FIXTURES_DIR=${FIXTURES_DIR:-/opt/test_fixtures}
FIXTURES_NO_PLATON_DIR=${FIXTURES_NO_PLATON_DIR:-/opt/test_fixtures_no_platon}
PY_BIN=${PY_BIN:-/opt/conda/bin/python3}
PY_BIN_OVERRIDE=${PY_BIN_OVERRIDE:-}
STATE_DIR=${STATE_DIR:-/var/smoke}
SMOKE_NJOBS=${SMOKE_NJOBS:-4}
SMOKE_TOTAL_TIMEOUT=${SMOKE_TOTAL_TIMEOUT:-600}
WITHOUT_PLATON=${WITHOUT_PLATON:-}

# Pick the sanity_results.csv fixture matching the build mode.
if [ -n "$WITHOUT_PLATON" ]; then
  FIXTURE_SANITY_CSV="$FIXTURES_NO_PLATON_DIR/sanity_results.csv"
  FIXTURE_STAGE1_CSV="$FIXTURES_NO_PLATON_DIR/stage1/sanity_results.csv"
else
  FIXTURE_SANITY_CSV="$FIXTURES_DIR/sanity_results.csv"
  FIXTURE_STAGE1_CSV="$FIXTURES_DIR/stage1/sanity_results.csv"
fi
export WITHOUT_PLATON FIXTURE_SANITY_CSV FIXTURE_STAGE1_CSV

SINGLE_OUTPUT="$STATE_DIR/single-output"

if [ ! -x "$ENTRYPOINT_BIN" ]; then
  echo "Entrypoint script not found or not executable: $ENTRYPOINT_BIN" >&2
  exit 1
fi
if [ ! -d "$DATA_DIR" ]; then
  echo "Test data directory not found: $DATA_DIR" >&2
  exit 1
fi
if [ ! -d "$FIXTURES_DIR" ]; then
  echo "Fixtures directory not found: $FIXTURES_DIR" >&2
  exit 1
fi
if [ -n "$PY_BIN_OVERRIDE" ]; then
  if ! command -v "$PY_BIN_OVERRIDE" >/dev/null 2>&1; then
    echo "Requested Python override not available: $PY_BIN_OVERRIDE" >&2
    exit 1
  fi
  PY_BIN="$PY_BIN_OVERRIDE"
fi
export PY_BIN

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1


step_prepare() {
  rm -rf "$STATE_DIR"
  mkdir -p "$STATE_DIR" "$SINGLE_OUTPUT"
  echo "[smoke prepare] state dir ready at $STATE_DIR"
}

step_help() {
  "$ENTRYPOINT_BIN" --help >/dev/null
  echo "[smoke help] entrypoint --help renders OK"
}

step_sanity_only() {
  "$ENTRYPOINT_BIN" \
    --run-sanity-only \
    --sanity-total-timeout "$SMOKE_TOTAL_TIMEOUT" \
    --output-dir "$SINGLE_OUTPUT" \
    --n-jobs "$SMOKE_NJOBS" \
    -i "$DATA_DIR"
  if [ ! -d "$SINGLE_OUTPUT/results" ]; then
    echo "Sanity stage produced no results/ directory" >&2
    exit 1
  fi
  n_dirs=$(find "$SINGLE_OUTPUT/results" -mindepth 1 -maxdepth 1 -type d | wc -l)
  echo "[smoke sanity-only] $n_dirs CIF result dir(s) produced"
}

step_check_stage_jsons() {
  "$PY_BIN" - "$SINGLE_OUTPUT/results" "$FIXTURES_DIR/results" <<'PY'
import json, os, sys
from pathlib import Path

actual_root = Path(sys.argv[1])
fixture_root = Path(sys.argv[2])
without_platon = bool(os.environ.get("WITHOUT_PLATON"))

# libcif's stage JSON also contains a `libcif_out` string with formatted floats
# (cell parameters, density, volume). Only the boolean validity flag is stable
# across runs, so we compare just that field.
LIBCIF_FIELDS = ("libcif_validity",)
# PLATON's stage JSON contains a `platon_data` list with the full PLATON report:
# header with binary build date and check.def version, alert wordings — these
# vary between PLATON builds and check.def updates. Only the boolean validity
# flag is stable across builds, so we compare just that field (like for libcif).
PLATON_FIELDS = ("platon_validity",)
# Stages compared in full (no floats in their schema).
FULL_MATCH_STAGES = ("platon", "oxichecker", "framechecker")
ALL_STAGES = FULL_MATCH_STAGES + ("libcif",)
# In WITHOUT_PLATON builds the runner never produces __platon.json, so we
# exclude the platon stage from the comparison entirely (the shared fixture
# tree still has __platon.json files for the WITH_PLATON build).
if without_platon:
    FULL_MATCH_STAGES = tuple(s for s in FULL_MATCH_STAGES if s != "platon")
    ALL_STAGES = tuple(s for s in ALL_STAGES if s != "platon")


def load(path):
    with open(path) as f:
        return json.load(f)


fixture_cifs = {p.name for p in fixture_root.iterdir() if p.is_dir()}
actual_cifs = {p.name for p in actual_root.iterdir() if p.is_dir()}

errors = []
if fixture_cifs - actual_cifs:
    errors.append(f"actual is missing CIF result dirs: {sorted(fixture_cifs - actual_cifs)}")
if actual_cifs - fixture_cifs:
    errors.append(f"actual has extra CIF result dirs: {sorted(actual_cifs - fixture_cifs)}")

n_files = 0
for cif in sorted(fixture_cifs & actual_cifs):
    f_dir = fixture_root / cif
    a_dir = actual_root / cif
    for stage in ALL_STAGES:
        fname = f"{cif}__{stage}.json"
        f_path = f_dir / fname
        a_path = a_dir / fname
        if f_path.exists() and not a_path.exists():
            errors.append(f"{cif}: actual is missing {fname}")
            continue
        if a_path.exists() and not f_path.exists():
            errors.append(f"{cif}: actual has unexpected {fname}")
            continue
        if not f_path.exists():
            continue
        n_files += 1
        exp = load(f_path)
        got = load(a_path)
        if stage == "libcif":
            for field in LIBCIF_FIELDS:
                if got.get(field) != exp.get(field):
                    errors.append(
                        f"FIXTURE MISMATCH [{cif} / libcif / {field}]\n"
                        f"  expected: {exp.get(field)!r}\n"
                        f"  actual:   {got.get(field)!r}"
                    )
        elif stage == "platon":
            for field in PLATON_FIELDS:
                if got.get(field) != exp.get(field):
                    errors.append(
                        f"FIXTURE MISMATCH [{cif} / platon / {field}]\n"
                        f"  expected: {exp.get(field)!r}\n"
                        f"  actual:   {got.get(field)!r}"
                    )
        else:
            if got != exp:
                errors.append(
                    f"FIXTURE MISMATCH [{cif} / {stage}]\n"
                    f"  expected: {exp!r}\n"
                    f"  actual:   {got!r}"
                )

if without_platon:
    leftover = sorted(p.name for p in actual_root.glob("*/*__platon.json"))
    if leftover:
        errors.append(
            f"WITHOUT_PLATON build produced unexpected __platon.json files: {leftover[:5]}"
            + (f" ... and {len(leftover) - 5} more" if len(leftover) > 5 else "")
        )

if errors:
    for err in errors[:20]:
        sys.stderr.write(err + "\n")
    if len(errors) > 20:
        sys.stderr.write(f"... and {len(errors) - 20} more error(s)\n")
    raise SystemExit(1)

mode = "WITHOUT_PLATON" if without_platon else "with PLATON"
print(
    f"[smoke check-stage-jsons] all stage-JSON fixtures match ({mode}; "
    f"{n_files} files across {len(fixture_cifs & actual_cifs)} CIF dirs)"
)
PY
}

step_oxichecker_only() {
  cd "$STATE_DIR"
  "$ENTRYPOINT_BIN" \
    --run-oxichecker-only \
    --n-jobs "$SMOKE_NJOBS" \
    -i "$SINGLE_OUTPUT"
  if [ ! -f "$STATE_DIR/oxichecker_results.csv" ]; then
    echo "OxiChecker run did not produce oxichecker_results.csv" >&2
    exit 1
  fi
  echo "[smoke oxichecker-only] oxichecker_results.csv produced"
}

step_check_oxichecker_csv() {
  "$PY_BIN" - "$STATE_DIR/oxichecker_results.csv" "$FIXTURES_DIR/oxichecker_results.csv" <<'PY'
import csv, sys


def read_indexed(path, key="cif"):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    seen = {}
    for r in rows:
        if r[key] in seen:
            raise SystemExit(f"duplicate {key}={r[key]!r} in {path}")
        seen[r[key]] = r
    return seen


actual_path, expected_path = sys.argv[1], sys.argv[2]
got = read_indexed(actual_path)
exp = read_indexed(expected_path)

missing = set(exp) - set(got)
extra = set(got) - set(exp)
if missing:
    raise SystemExit(f"oxichecker_results.csv missing rows for: {sorted(missing)[:10]}")
if extra:
    raise SystemExit(f"oxichecker_results.csv unexpected rows for: {sorted(extra)[:10]}")

errors = []
for cif in sorted(exp):
    if got[cif] != exp[cif]:
        diff_cols = [k for k in exp[cif] if got[cif].get(k) != exp[cif].get(k)]
        errors.append(
            f"oxichecker_results.csv mismatch for {cif!r} on cols {diff_cols}:\n"
            f"  expected: { {k: exp[cif][k] for k in diff_cols} }\n"
            f"  actual:   { {k: got[cif].get(k) for k in diff_cols} }"
        )
if errors:
    for e in errors[:5]:
        sys.stderr.write(e + "\n")
    if len(errors) > 5:
        sys.stderr.write(f"... and {len(errors)-5} more row mismatch(es)\n")
    raise SystemExit(1)

print(f"[smoke check-oxichecker-csv] CSV matches fixture ({len(got)} row(s))")
PY
}

step_postprocess_only() {
  "$ENTRYPOINT_BIN" \
    --run-postprocess-only \
    --n-jobs "$SMOKE_NJOBS" \
    --output-dir "$SINGLE_OUTPUT"
  if [ ! -f "$SINGLE_OUTPUT/sanity_results.csv" ]; then
    echo "Postprocessing run did not produce sanity_results.csv" >&2
    exit 1
  fi
  echo "[smoke postprocess-only] sanity_results.csv produced"
}

step_check_sanity_csv() {
  "$PY_BIN" - "$SINGLE_OUTPUT/sanity_results.csv" "$FIXTURE_SANITY_CSV" <<'PY'
import csv, sys

# Text/boolean/integer/hash columns only. Floating-point columns
# (`density`, `volume`, `reduced_formula`) are intentionally skipped because
# their last-bit values are not bit-stable across runs.
TEXT_COLS = [
    "cif",
    "Basic Validity",
    "LibCif Validity",
    "LibCif_Warning",
    "PLATON Validity",
    "OxiChecker Validity",
    "Sanity",
    "content_hash",
    "formula",
    "group_str",
    "group_id",
    "structure_hash_strict",
    "structure_hash",
    "is_graph_constructed",
    "graph_dim",
    "HAS_OMS",
    "DECORATED_GRAPH_HASH",
    "UNDECORATED_GRAPH_HASH",
    "DECORATED_SCAFFOLD_HASH",
    "UNDECORATED_SCAFFOLD_HASH",
    "has atomic overlaps",
    "has overcoordinated c",
    "has overcoordinated h",
    "has overcoordinated n",
    "has lone molecule",
    "has bad rare earth",
    "has bad alkali alkaline",
    "has bad terminal oxo",
]


def read_indexed(path, key="cif"):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    seen = {}
    for r in rows:
        if r[key] in seen:
            raise SystemExit(f"duplicate {key}={r[key]!r} in {path}")
        seen[r[key]] = r
    return seen


actual_path, expected_path = sys.argv[1], sys.argv[2]
got = read_indexed(actual_path)
exp = read_indexed(expected_path)

missing = set(exp) - set(got)
extra = set(got) - set(exp)
if missing:
    raise SystemExit(f"sanity_results.csv missing rows for: {sorted(missing)[:10]}")
if extra:
    raise SystemExit(f"sanity_results.csv unexpected rows for: {sorted(extra)[:10]}")

errors = []
checked_cols = []
for col in TEXT_COLS:
    if col not in next(iter(exp.values()), {}):
        # Column absent in fixture (e.g. PLATON Validity in WITHOUT_PLATON
        # builds, or the extra graph-hash cols in stage1-only fixtures) —
        # skip it instead of treating it as a fixture error.
        continue
    checked_cols.append(col)

for cif in sorted(exp):
    g, e = got[cif], exp[cif]
    for col in checked_cols:
        if col not in g:
            raise SystemExit(f"actual output is missing column {col!r}")
        if g[col] != e[col]:
            errors.append(
                f"{cif} column {col!r}: expected={e[col]!r}, actual={g[col]!r}"
            )
if errors:
    for e in errors[:10]:
        sys.stderr.write(e + "\n")
    if len(errors) > 10:
        sys.stderr.write(f"... and {len(errors)-10} more cell mismatch(es)\n")
    raise SystemExit(1)

print(
    f"[smoke check-sanity-csv] sanity_results.csv text columns match "
    f"({len(got)} row(s), {len(checked_cols)} columns checked)"
)
PY
}


step_check_help_content() {
  out=$("$ENTRYPOINT_BIN" --help 2>&1)
  for section in \
      "Run mode options:" \
      "Sanity runner options (prefix --sanity-):" \
      "OxiChecker options (prefix --oxichecker-):" \
      "Postprocessing options (prefix --post-):" \
      "Combined options:"; do
    if ! grep -qF "$section" <<<"$out"; then
      echo "Missing section in --help output: $section" >&2
      exit 1
    fi
  done
  if [ -n "$WITHOUT_PLATON" ]; then
    if grep -qi "platon" <<<"$out"; then
      echo "WITHOUT_PLATON build but --help still mentions PLATON:" >&2
      grep -i "platon" <<<"$out" >&2
      exit 1
    fi
    echo "[smoke check-help-content] all 5 sections present; no PLATON mentions (WITHOUT_PLATON)"
  else
    echo "[smoke check-help-content] all 5 expected sections present in --help"
  fi
}

step_negative_unknown_flag() {
  set +e
  out=$("$ENTRYPOINT_BIN" --not-a-real-flag 2>&1)
  rc=$?
  set -e
  if [ "$rc" -ne 2 ]; then
    echo "Expected exit code 2 for unknown flag, got $rc" >&2
    echo "Output: $out" >&2
    exit 1
  fi
  if ! grep -qi "Unknown option" <<<"$out"; then
    echo "Expected stderr to contain 'Unknown option', got: $out" >&2
    exit 1
  fi
  echo "[smoke negative-unknown-flag] exits 2 with 'Unknown option'"
}

step_negative_missing_input() {
  set +e
  out=$("$ENTRYPOINT_BIN" --run-sanity-only 2>&1)
  rc=$?
  set -e
  if [ "$rc" -ne 2 ]; then
    echo "Expected exit code 2 for missing --input, got $rc" >&2
    echo "Output: $out" >&2
    exit 1
  fi
  if ! grep -qi "input is required" <<<"$out"; then
    echo "Expected stderr to mention that --input is required, got: $out" >&2
    exit 1
  fi
  echo "[smoke negative-missing-input] exits 2 with 'input is required'"
}

step_stage1_only() {
  s1_dir="$STATE_DIR/stage1"
  rm -rf "$s1_dir"
  "$ENTRYPOINT_BIN" \
    --run-stage1-only \
    --output-dir "$s1_dir" \
    --n-jobs 1 \
    -i "$DATA_DIR/qmof-da78929.cif"

  diff -q \
    "$s1_dir/results/qmof-da78929/qmof-da78929__libcif.json" \
    "$FIXTURES_DIR/results/qmof-da78929/qmof-da78929__libcif.json" >/dev/null \
    || { echo "stage1-only: libcif.json differs from main fixture" >&2; exit 1; }

  "$PY_BIN" - "$s1_dir/sanity_results.csv" "$FIXTURE_STAGE1_CSV" <<'PY'
import csv, sys

# Stage1-only CSV has a slimmer schema (no graph hashes / HAS_OMS columns).
# Floating-point columns are excluded as in the full-pipeline check.
TEXT_COLS = [
    "cif",
    "Basic Validity",
    "LibCif Validity",
    "LibCif_Warning",
    "PLATON Validity",
    "OxiChecker Validity",
    "Sanity",
    "content_hash",
    "formula",
    "group_str",
    "group_id",
    "structure_hash_strict",
    "structure_hash",
    "is_graph_constructed",
    "graph_dim",
    "has atomic overlaps",
    "has overcoordinated c",
    "has overcoordinated h",
    "has overcoordinated n",
    "has lone molecule",
    "has bad rare earth",
    "has bad alkali alkaline",
    "has bad terminal oxo",
]


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


got = read(sys.argv[1])
exp = read(sys.argv[2])

if len(got) != 1 or len(exp) != 1:
    raise SystemExit(
        f"stage1-only: expected exactly 1 row each, got {len(got)} actual / {len(exp)} expected"
    )
g, e = got[0], exp[0]
errors = []
checked = []
for col in TEXT_COLS:
    if col not in e:
        # PLATON Validity is absent in the WITHOUT_PLATON stage1 fixture; skip.
        continue
    if col not in g:
        raise SystemExit(f"stage1 actual output is missing column {col!r}")
    checked.append(col)
    if g[col] != e[col]:
        errors.append(f"  column {col!r}: expected={e[col]!r}, actual={g[col]!r}")
if errors:
    sys.stderr.write("stage1-only sanity_results.csv mismatches:\n")
    for line in errors:
        sys.stderr.write(line + "\n")
    raise SystemExit(1)
print(
    f"[smoke stage1-only] sanity_results.csv text columns match stage1 fixture "
    f"({len(checked)} columns checked)"
)
PY
}

step_single_file_input() {
  sf_dir="$STATE_DIR/single-file"
  rm -rf "$sf_dir"
  "$ENTRYPOINT_BIN" \
    --run-sanity-only \
    --sanity-total-timeout "$SMOKE_TOTAL_TIMEOUT" \
    --output-dir "$sf_dir" \
    --n-jobs 1 \
    -i "$DATA_DIR/qmof-da78929.cif"

  if [ ! -d "$sf_dir/results/qmof-da78929" ] ; then
    echo "single-file-input: missing result dir for the single CIF" >&2
    exit 1
  fi

  "$PY_BIN" - \
      "$sf_dir/results/qmof-da78929" \
      "$FIXTURES_DIR/results/qmof-da78929" <<'PY'
import json, os, sys
from pathlib import Path

actual = Path(sys.argv[1])
fixture = Path(sys.argv[2])
without_platon = bool(os.environ.get("WITHOUT_PLATON"))

LIBCIF_FIELDS = ("libcif_validity",)
PLATON_FIELDS = ("platon_validity",)
FULL_MATCH_STAGES = ("oxichecker", "framechecker") if without_platon else ("platon", "oxichecker", "framechecker")

errors = []
n = 0
for stage in FULL_MATCH_STAGES + ("libcif",):
    fname = f"qmof-da78929__{stage}.json"
    a = actual / fname
    e = fixture / fname
    if not e.exists():
        continue
    if not a.exists():
        errors.append(f"single-file: missing {fname}")
        continue
    with open(a) as f: got = json.load(f)
    with open(e) as f: exp = json.load(f)
    n += 1
    if stage == "libcif":
        for field in LIBCIF_FIELDS:
            if got.get(field) != exp.get(field):
                errors.append(f"single-file libcif.{field} mismatch")
    elif stage == "platon":
        for field in PLATON_FIELDS:
            if got.get(field) != exp.get(field):
                errors.append(f"single-file platon.{field} mismatch")
    elif got != exp:
        errors.append(f"single-file {stage}.json content mismatch")
if errors:
    for err in errors:
        sys.stderr.write(err + "\n")
    raise SystemExit(1)
print(f"[smoke single-file-input] {n} stage-JSONs match the main fixtures")
PY
}

step_recursive() {
  nested="$STATE_DIR/nested-input"
  rm -rf "$nested"
  mkdir -p "$nested/sub1/sub2"
  cp "$DATA_DIR/qmof-da78929.cif" "$nested/sub1/sub2/"

  out_dir="$STATE_DIR/recursive-out"
  rm -rf "$out_dir"

  "$ENTRYPOINT_BIN" \
    --run-stage1-only \
    --sanity-recursive \
    --output-dir "$out_dir" \
    --n-jobs 1 \
    -i "$nested"

  if [ ! -d "$out_dir/results/qmof-da78929" ] ; then
    echo "recursive: nested CIF was not picked up (no result dir)" >&2
    exit 1
  fi
  echo "[smoke recursive] --sanity-recursive descended into sub1/sub2/"
}

step_custom_post_output() {
  custom="$STATE_DIR/custom-result.csv"
  rm -f "$custom"
  "$ENTRYPOINT_BIN" \
    --run-postprocess-only \
    --output-dir "$SINGLE_OUTPUT" \
    --post-output "$custom" \
    --n-jobs 1
  if [ ! -f "$custom" ] ; then
    echo "custom-post-output: CSV not created at $custom" >&2
    exit 1
  fi
  if ! head -n 1 "$custom" | grep -q "^cif," ; then
    echo "custom-post-output: file at $custom does not look like a sanity CSV" >&2
    head -n 1 "$custom" >&2
    exit 1
  fi
  echo "[smoke custom-post-output] CSV written to custom --post-output path"
}

step_no_obabel_run() {
  no_dir="$STATE_DIR/no-obabel"
  rm -rf "$no_dir"
  "$ENTRYPOINT_BIN" \
    --run-sanity-only \
    --sanity-no-obabel-run \
    --sanity-total-timeout "$SMOKE_TOTAL_TIMEOUT" \
    --output-dir "$no_dir" \
    --n-jobs 1 \
    -i "$DATA_DIR/qmof-da78929.cif"

  if compgen -G "$no_dir/babel_cifs/OBABEL_*.cif" >/dev/null ; then
    echo "no-obabel-run: babel_cifs/ should be empty of OBABEL_*.cif but isn't" >&2
    ls -la "$no_dir/babel_cifs/" >&2 || true
    exit 1
  fi
  echo "[smoke no-obabel-run] babel_cifs/ contains no OBABEL_*.cif (preprocessing skipped)"
}

step_timeout_soft() {
  # Test that --sanity-total-timeout actually cuts off per-structure
  # processing. T is set well below the natural per-CIF processing time
  # (~8s on qmof-da78929) so the timeout must fire. We then assert that the
  # observed per-CIF processing time is within the user-tolerated band
  # T*1.2 + 2 seconds — "soft" with respect to the gap between the configured
  # timeout and the actual interrupt moment.
  local T="${SMOKE_TIMEOUT_TEST_T:-1}"
  local bound=$(( T * 12 / 10 + 2 ))   # T*1.2 + 2 (integer floor)
  local to_dir="$STATE_DIR/timeout-soft"
  local logfile="$STATE_DIR/timeout-soft.log"
  rm -rf "$to_dir"

  set +e
  "$ENTRYPOINT_BIN" \
    --run-sanity-only \
    --sanity-total-timeout "$T" \
    --output-dir "$to_dir" \
    --n-jobs 1 \
    -i "$DATA_DIR/qmof-da78929.cif" >"$logfile" 2>&1
  local rc=$?
  set -e

  if [ "$rc" -ne 0 ]; then
    echo "timeout-soft: pipeline exited non-zero ($rc)" >&2
    cat "$logfile" >&2
    exit 1
  fi

  # The pipeline must have hit the per-structure timeout, otherwise the test
  # didn't actually exercise the cut-off (T was too generous for this CIF).
  if ! grep -qE "TIMEOUT|Total execution timeout" \
       "$to_dir/results/qmof-da78929/qmof-da78929.log" 2>/dev/null; then
    echo "timeout-soft: T=${T}s did not trigger the per-structure timeout — pick a smaller T" >&2
    cat "$logfile" >&2
    exit 1
  fi

  # Pull per-CIF processing time from the final tqdm line "N/N [MM:SS<...]"
  # (its first MM:SS is the elapsed). Parsing in Python keeps the regex sane.
  local elapsed
  elapsed=$("$PY_BIN" - "$logfile" <<'PY'
import re, sys
text = open(sys.argv[1]).read()
m = re.findall(r'1/1 \[(\d+):(\d+)<', text)
if not m:
    sys.exit("could not find tqdm progress line")
mm, ss = m[-1]
print(int(mm) * 60 + int(ss))
PY
  )

  if [ "$elapsed" -gt "$bound" ]; then
    echo "timeout-soft: per-CIF processing ${elapsed}s exceeds bound ${bound}s (T=${T}s, T*1.2+2)" >&2
    cat "$logfile" >&2
    exit 1
  fi
  echo "[smoke timeout-soft] timeout fired; per-CIF time=${elapsed}s within bound=${bound}s (T=${T}s, T*1.2+2)"
}

step_idempotency() {
  a_dir="$STATE_DIR/idem-a"
  b_dir="$STATE_DIR/idem-b"
  rm -rf "$a_dir" "$b_dir"
  for d in "$a_dir" "$b_dir"; do
    "$ENTRYPOINT_BIN" \
      --run-stage1-only \
      --output-dir "$d" \
      --n-jobs 1 \
      -i "$DATA_DIR/qmof-da78929.cif"
  done
  diff -q \
    "$a_dir/results/qmof-da78929/qmof-da78929__libcif.json" \
    "$b_dir/results/qmof-da78929/qmof-da78929__libcif.json" >/dev/null \
    || { echo "idempotency: libcif.json differs between consecutive stage1 runs" >&2; exit 1; }
  echo "[smoke idempotency] two consecutive stage1 runs produce identical libcif.json"
}


step="${1:-all}"
case "$step" in
  prepare)               step_prepare ;;
  help)                  step_help ;;
  check-help-content)    step_check_help_content ;;
  negative-unknown-flag) step_negative_unknown_flag ;;
  negative-missing-input) step_negative_missing_input ;;
  sanity-only)           step_sanity_only ;;
  check-stage-jsons)     step_check_stage_jsons ;;
  oxichecker-only)       step_oxichecker_only ;;
  check-oxichecker-csv)  step_check_oxichecker_csv ;;
  postprocess-only)      step_postprocess_only ;;
  check-sanity-csv)      step_check_sanity_csv ;;
  stage1-only)           step_stage1_only ;;
  single-file-input)     step_single_file_input ;;
  recursive)             step_recursive ;;
  custom-post-output)    step_custom_post_output ;;
  no-obabel-run)         step_no_obabel_run ;;
  timeout-soft)          step_timeout_soft ;;
  idempotency)           step_idempotency ;;
  all)
    step_prepare
    step_help
    step_check_help_content
    step_negative_unknown_flag
    step_negative_missing_input
    step_sanity_only
    step_check_stage_jsons
    step_oxichecker_only
    step_check_oxichecker_csv
    step_postprocess_only
    step_check_sanity_csv
    step_stage1_only
    step_single_file_input
    step_recursive
    step_custom_post_output
    step_no_obabel_run
    step_timeout_soft
    step_idempotency
    echo "Entrypoint CLI smoke tests completed successfully"
    ;;
  *)
    echo "Unknown smoke step: $step" >&2
    echo "Valid steps: prepare help check-help-content negative-unknown-flag negative-missing-input sanity-only check-stage-jsons oxichecker-only check-oxichecker-csv postprocess-only check-sanity-csv stage1-only single-file-input recursive custom-post-output no-obabel-run timeout-soft idempotency all" >&2
    exit 2
    ;;
esac
