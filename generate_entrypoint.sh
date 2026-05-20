#!/usr/bin/env bash
set -euo pipefail

SANITY_BIN=${SANITY_BIN:-/code/sanity_runner}
PY_BIN=${PY_BIN:-/opt/conda/bin/python3}
OXI_SCRIPT=${OXI_SCRIPT:-/opt/oxichecker/main.py}
POST_SCRIPT=${POST_SCRIPT:-/code/postprocessing_table_creation.py}
ENTRYPOINT_PATH=${ENTRYPOINT_PATH:-/entrypoint.sh}

SANITY_HELP="$($SANITY_BIN --help 2>&1 || true)"
OXI_HELP="$($PY_BIN -W ignore::SyntaxWarning "$OXI_SCRIPT" --help 2>&1 || true)"
POST_HELP="$($PY_BIN "$POST_SCRIPT" --help 2>&1 || true)"

export SANITY_BIN PY_BIN OXI_SCRIPT POST_SCRIPT ENTRYPOINT_PATH SANITY_HELP OXI_HELP POST_HELP

"$PY_BIN" <<'PY'
import os
import re
import textwrap
from pathlib import Path

SANITY_BIN = os.environ["SANITY_BIN"]
PY_BIN = os.environ["PY_BIN"]
OXI_SCRIPT = os.environ["OXI_SCRIPT"]
POST_SCRIPT = os.environ["POST_SCRIPT"]
ENTRYPOINT_PATH = Path(os.environ["ENTRYPOINT_PATH"])
SANITY_HELP = os.environ.get("SANITY_HELP", "")
OXI_HELP = os.environ.get("OXI_HELP", "")
POST_HELP = os.environ.get("POST_HELP", "")

OPT_LINE_RE = re.compile(r"^\s*(?:-[^,\s]+,\s*)?(--[A-Za-z0-9][\w-]*)")

def token_implies_value(token: str) -> bool:
    token = token.strip()
    if not token:
        return False
    if token.startswith("{") or token.startswith("<"):
        return True
    stripped = token.strip("[](){}<>,.")
    if not stripped:
        return False
    if stripped.upper() == stripped and any(c.isalpha() for c in stripped):
        return True
    if stripped.endswith("..."):
        core = stripped[:-3]
        return core.upper() == core and bool(core)
    return False

def parse_options(help_text: str):
    options = []
    seen = set()
    lines = help_text.splitlines()
    
    for i, line in enumerate(lines):
        lstrip = line.lstrip()
        if not lstrip.startswith(('-', '–')):
            continue
        match = OPT_LINE_RE.match(line)
        if not match:
            continue
        long_opt = match.group(1)
        if long_opt in seen:
            continue
        seen.add(long_opt)
        after = line[line.find(long_opt) + len(long_opt):].strip()
        takes_value = False
        if after:
            token = after.split()[0]
            takes_value = token_implies_value(token)
        
        # Extract full description including continuation lines
        detail = line.strip()
        parts = re.split(r"\s{2,}", detail, maxsplit=1)
        description = parts[1] if len(parts) == 2 else ""
        
        # Look for continuation lines
        j = i + 1
        while j < len(lines):
            next_line = lines[j]
            # If next line starts with option pattern, stop
            if OPT_LINE_RE.match(next_line):
                break
            # If next line is indented and not empty, it's likely a continuation
            if next_line.startswith('  ') and next_line.strip() and not next_line.lstrip().startswith('-'):
                if description:
                    description += " " + next_line.strip()
                else:
                    description = next_line.strip()
                j += 1
            else:
                break
        
        options.append({
            "long": long_opt,
            "takes_value": takes_value,
            "description": description.strip(),
        })
    return options

def parse_positionals(help_text: str):
    positionals = []
    usage_lines = []
    for line in help_text.splitlines():
        if line.lower().startswith("usage"):
            usage_lines.append(line)
        elif usage_lines and (line.startswith(" ") or line.startswith("\t")):
            stripped = line.strip()
            if stripped:
                usage_lines.append(stripped)
            else:
                break
        elif usage_lines:
            break
    usage_text = " ".join(usage_lines)
    for name, dots in re.findall(r"\[([A-Z][A-Z0-9_]+)\](\.{3})?", usage_text):
        if name == "OPTIONS":
            continue
        positionals.append({
            "name": name,
            "variadic": bool(dots)
        })
    return positionals

sanity_options = parse_options(SANITY_HELP)
sanity_positionals = parse_positionals(SANITY_HELP)
oxichecker_options = parse_options(OXI_HELP)
post_options = parse_options(POST_HELP)

# Remove options handled via combined flags generated below
sanity_positionals = [pos for pos in sanity_positionals if pos["name"] != "INPUT_FILE"]
sanity_options = [opt for opt in sanity_options if opt["long"] not in {"--max-workers", "--output-dir"}]
oxichecker_options = [opt for opt in oxichecker_options if opt["long"] not in {"--input", "--num_cores"}]
post_options = [opt for opt in post_options if opt["long"] not in {"--folder", "--n-jobs"}]

if not SANITY_HELP.strip():
    raise SystemExit("Failed to obtain --help from sanity runner")
if not OXI_HELP.strip():
    raise SystemExit("Failed to obtain --help from OxiChecker script")

lines = []
lines.append("#!/usr/bin/env bash")
lines.append("set -euo pipefail")
lines.append("")
lines.append(f"SANITY_BIN=\"{SANITY_BIN}\"")
lines.append(f"PY_BIN=\"{PY_BIN}\"")
lines.append(f"OXI_SCRIPT=\"{OXI_SCRIPT}\"")
lines.append(f"POST_SCRIPT=\"{POST_SCRIPT}\"")
lines.append("")
lines.append("entrypoint_help() { cat <<'EOF'")
lines.append("Usage: entrypoint [OPTIONS]")
lines.append("")
lines.append("Run mode options:")
lines.append("  --run-sanity-only               Run only sanity runner")
lines.append("  --run-oxichecker-only               Run only OxiChecker")
lines.append("  --run-postprocess-only          Run only postprocessing")
lines.append("")
lines.append("Sanity runner options (prefix --sanity-):")
for opt in sanity_options:
    pref = f"--sanity-{opt['long'][2:]}"
    display = pref + (" VALUE" if opt["takes_value"] else "")
    desc = opt["description"]
    if desc:
        lines.append(f"  {display:<35} {desc}")
    else:
        lines.append(f"  {display}")
lines.append("")
lines.append("Combined options:")
lines.append("  -i, --input VALUE               Shared input path for both runners")
lines.append("  --n-jobs VALUE                  Shared worker count (sanity --max-workers, OxiChecker --num_cores)")
lines.append("  --output-dir VALUE              Output directory for sanity runner (default: <inputname>_sanity_results)")
lines.append("")
lines.append("OxiChecker options (prefix --oxichecker-):")
for opt in oxichecker_options:
    pref = f"--oxichecker-{opt['long'][2:]}"
    display = pref + (" VALUE" if opt["takes_value"] else "")
    desc = opt["description"]
    if desc:
        lines.append(f"  {display:<35} {desc}")
    else:
        lines.append(f"  {display}")
lines.append("")
lines.append("Postprocessing options (prefix --post-):")
for opt in post_options:
    pref = f"--post-{opt['long'][2:]}"
    display = pref + (" VALUE" if opt["takes_value"] else "")
    desc = opt["description"]
    if desc:
        lines.append(f"  {display:<35} {desc}")
    else:
        lines.append(f"  {display}")
lines.append("")
lines.append("Use --help or -h to display this message.")
lines.append("EOF")
lines.append("}")
lines.append("")
lines.append("declare -a SANITY_ARGS=()")
lines.append("declare -a SANITY_POSITIONALS=()")
lines.append("declare -a OXICHECKER_ARGS=()")
lines.append("declare -a POST_ARGS=()")
lines.append("oxichecker_input=\"\"")
lines.append("n_jobs_value=\"\"")
lines.append("output_dir_value=\"\"")
lines.append("run_mode=\"all\"  # all, sanity, oxichecker, postprocess")
lines.append("")
lines.append("while [[ $# -gt 0 ]]; do")
lines.append("  case \"$1\" in")
lines.append("    --help|-h)")
lines.append("      entrypoint_help")
lines.append("      exit 0")
lines.append("      ;;")
lines.append("    --run-sanity-only)")
lines.append("      run_mode=\"sanity\"")
lines.append("      shift")
lines.append("      continue")
lines.append("      ;;")
lines.append("    --run-oxichecker-only)")
lines.append("      run_mode=\"oxichecker\"")
lines.append("      shift")
lines.append("      continue")
lines.append("      ;;")
lines.append("    --run-postprocess-only)")
lines.append("      run_mode=\"postprocess\"")
lines.append("      shift")
lines.append("      continue")
lines.append("      ;;")
lines.append("    -i|--input)")
lines.append("      [[ $# -ge 2 ]] || { echo \"Missing value for --input\" >&2; exit 1; }")
lines.append("      value=$2")
lines.append("      SANITY_POSITIONALS+=(\"$value\")")
lines.append("      oxichecker_input=\"$value\"")
lines.append("      shift 2")
lines.append("      continue")
lines.append("      ;;")
lines.append("    --input=*)")
lines.append("      value=${1#*=}")
lines.append("      SANITY_POSITIONALS+=(\"$value\")")
lines.append("      oxichecker_input=\"$value\"")
lines.append("      shift")
lines.append("      continue")
lines.append("      ;;")
lines.append("    --n-jobs)")
lines.append("      [[ $# -ge 2 ]] || { echo \"Missing value for --n-jobs\" >&2; exit 1; }")
lines.append("      n_jobs_value=$2")
lines.append("      shift 2")
lines.append("      continue")
lines.append("      ;;")
lines.append("    --n-jobs=*)")
lines.append("      n_jobs_value=${1#*=}")
lines.append("      shift")
lines.append("      continue")
lines.append("      ;;")
lines.append("    --output-dir)")
lines.append("      [[ $# -ge 2 ]] || { echo \"Missing value for --output-dir\" >&2; exit 1; }")
lines.append("      output_dir_value=$2")
lines.append("      shift 2")
lines.append("      continue")
lines.append("      ;;")
lines.append("    --output-dir=*)")
lines.append("      output_dir_value=${1#*=}")
lines.append("      shift")
lines.append("      continue")
lines.append("      ;;")

for opt in sanity_options:
    long = opt["long"]
    pref = f"--sanity-{long[2:]}"
    if opt["takes_value"]:
        lines.append(f"    {pref})")
        lines.append("      [[ $# -ge 2 ]] || { echo \"Missing value for %s\" >&2; exit 1; }" % pref)
        lines.append("      SANITY_ARGS+=(\"%s\" \"$2\")" % long)
        lines.append("      shift 2")
        lines.append("      continue")
        lines.append("      ;;")
        lines.append(f"    {pref}=*)")
        lines.append("      value=${1#*=}")
        lines.append("      SANITY_ARGS+=(\"%s\" \"$value\")" % long)
        lines.append("      shift")
        lines.append("      continue")
        lines.append("      ;;")
    else:
        lines.append(f"    {pref})")
        lines.append("      SANITY_ARGS+=(\"%s\")" % long)
        lines.append("      shift")
        lines.append("      continue")
        lines.append("      ;;")

for pos in sanity_positionals:
    opt_name = pos['name'].lower().replace('_', '-')
    pref = f"--sanity-{opt_name}"
    lines.append(f"    {pref})")
    lines.append("      [[ $# -ge 2 ]] || { echo \"Missing value for %s\" >&2; exit 1; }" % pref)
    lines.append("      SANITY_POSITIONALS+=(\"$2\")")
    lines.append("      shift 2")
    lines.append("      continue")
    lines.append("      ;;")
    lines.append(f"    {pref}=*)")
    lines.append("      value=${1#*=}")
    lines.append("      SANITY_POSITIONALS+=(\"$value\")")
    lines.append("      shift")
    lines.append("      continue")
    lines.append("      ;;")

for opt in oxichecker_options:
    long = opt["long"]
    pref = f"--oxichecker-{long[2:]}"
    if opt["takes_value"]:
        lines.append(f"    {pref})")
        lines.append("      [[ $# -ge 2 ]] || { echo \"Missing value for %s\" >&2; exit 1; }" % pref)
        lines.append("      OXICHECKER_ARGS+=(\"%s\" \"$2\")" % long)
        lines.append("      shift 2")
        lines.append("      continue")
        lines.append("      ;;")
        lines.append(f"    {pref}=*)")
        lines.append("      value=${1#*=}")
        lines.append("      OXICHECKER_ARGS+=(\"%s\" \"$value\")" % long)
        lines.append("      shift")
        lines.append("      continue")
        lines.append("      ;;")
    else:
        lines.append(f"    {pref})")
        lines.append("      OXICHECKER_ARGS+=(\"%s\")" % long)
        lines.append("      shift")
        lines.append("      continue")
        lines.append("      ;;")

for opt in post_options:
    long = opt["long"]
    pref = f"--post-{long[2:]}"
    if opt["takes_value"]:
        lines.append(f"    {pref})")
        lines.append("      [[ $# -ge 2 ]] || { echo \"Missing value for %s\" >&2; exit 1; }" % pref)
        lines.append("      POST_ARGS+=(\"%s\" \"$2\")" % long)
        lines.append("      shift 2")
        lines.append("      continue")
        lines.append("      ;;")
        lines.append(f"    {pref}=*)")
        lines.append("      value=${1#*=}")
        lines.append("      POST_ARGS+=(\"%s\" \"$value\")" % long)
        lines.append("      shift")
        lines.append("      continue")
        lines.append("      ;;")
    else:
        lines.append(f"    {pref})")
        lines.append("      POST_ARGS+=(\"%s\")" % long)
        lines.append("      shift")
        lines.append("      continue")
        lines.append("      ;;")

lines.append("    --)")
lines.append("      shift")
lines.append("      break")
lines.append("      ;;")
lines.append("    *)")
lines.append("      echo \"Unknown option: $1\" >&2")
lines.append("      exit 2")
lines.append("      ;;")
lines.append("  esac")
lines.append("done")
lines.append("")
lines.append("if [[ $# -gt 0 ]]; then")
lines.append("  echo \"Unexpected arguments: $*\" >&2")
lines.append("  exit 2")
lines.append("fi")
lines.append("")
lines.append('# Validate input requirements')
lines.append('if [[ "$run_mode" != "postprocess" && -z "$oxichecker_input" ]]; then')
lines.append('  echo "Error: --input is required for sanity and oxichecker modes." >&2')
lines.append('  exit 2')
lines.append('fi')
lines.append('')
lines.append('# Set up common arguments')
lines.append('if [[ "$run_mode" == "sanity" || "$run_mode" == "all" ]]; then')
lines.append('  # Only add oxichecker_input if not already in SANITY_POSITIONALS')
lines.append('  input_already_added=false')
lines.append('  for arg in "${SANITY_POSITIONALS[@]}"; do')
lines.append('    if [[ "$arg" == "$oxichecker_input" ]]; then')
lines.append('      input_already_added=true')
lines.append('      break')
lines.append('    fi')
lines.append('  done')
lines.append('  if [[ "$input_already_added" == false && -n "$oxichecker_input" ]]; then')
lines.append('    SANITY_POSITIONALS+=("$oxichecker_input")')
lines.append('  fi')
lines.append('  if [[ -n "$n_jobs_value" ]]; then')
lines.append('    SANITY_ARGS+=("--max-workers" "$n_jobs_value")')
lines.append('  fi')
lines.append('  if [[ -n "$output_dir_value" ]]; then')
lines.append('    SANITY_ARGS+=("--output-dir" "$output_dir_value")')
lines.append('  fi')
lines.append('fi')
lines.append('')
lines.append('if [[ "$run_mode" == "oxichecker" || "$run_mode" == "all" ]]; then')
lines.append('  OXICHECKER_ARGS+=("--input" "$oxichecker_input")')
lines.append('  if [[ -n "$n_jobs_value" ]]; then')
lines.append('    OXICHECKER_ARGS+=("--num_cores" "$n_jobs_value")')
lines.append('  fi')
lines.append('fi')
lines.append('')
lines.append('if [[ "$run_mode" == "postprocess" || "$run_mode" == "all" ]]; then')
lines.append('  # Set up postprocessing arguments')
lines.append('  output_dir=""')
lines.append('  if [[ -n "$output_dir_value" ]]; then')
lines.append('    output_dir="$output_dir_value"')
lines.append('  else')
lines.append('    for i in "${!SANITY_ARGS[@]}"; do')
lines.append('      if [[ "${SANITY_ARGS[$i]}" == "--output-dir" ]]; then')
lines.append('        output_dir="${SANITY_ARGS[$((i+1))]}"')
lines.append('        break')
lines.append('      fi')
lines.append('    done')
lines.append('  fi')
lines.append('  ')
lines.append('  # If still empty, calculate the default output directory like sanity_runner does')
lines.append('  if [[ -z "$output_dir" && -n "$oxichecker_input" ]]; then')
lines.append('    if [[ -f "$oxichecker_input" ]]; then')
lines.append('      # Single file: use stem + _sanity_results')
lines.append('      input_stem=$(basename "$oxichecker_input" .cif)')
lines.append('      output_dir="${input_stem}_sanity_results"')
lines.append('    elif [[ -d "$oxichecker_input" ]]; then')
lines.append('      # Directory: use dirname + _sanity_results')
lines.append('      input_name=$(basename "$oxichecker_input")')
lines.append('      output_dir="${input_name}_sanity_results"')
lines.append('    else')
lines.append('      output_dir="input_sanity_results"')
lines.append('    fi')
lines.append('    # Make output_dir absolute if it\'s relative')
lines.append('    if [[ ! "$output_dir" = /* ]]; then')
lines.append('      output_dir="$(pwd)/$output_dir"')
lines.append('    fi')
lines.append('  fi')
lines.append('  ')
lines.append('  if [[ -n "$output_dir" ]]; then')
lines.append('    POST_ARGS+=("--folder" "$output_dir")')
lines.append('  fi')
lines.append('  ')
lines.append('  if [[ -n "$n_jobs_value" ]]; then')
lines.append('    POST_ARGS+=("--n-jobs" "$n_jobs_value")')
lines.append('  fi')
lines.append('  ')
lines.append('  # Check for OxiChecker output file in the output directory')
lines.append('  if [[ -n "$output_dir" && -f "$output_dir/oxichecker_results.csv" ]]; then')
lines.append('    POST_ARGS+=("--oxichecker-path" "$output_dir/oxichecker_results.csv")')
lines.append('  fi')
lines.append('fi')
lines.append("")
lines.append('# Check for help requests')
lines.append('sanity_help_requested=false')
lines.append('oxichecker_help_requested=false')
lines.append('post_help_requested=false')
lines.append('')
lines.append('for arg in "${SANITY_ARGS[@]}"; do')
lines.append('  if [[ "$arg" == "--help" || "$arg" == "-h" ]]; then')
lines.append('    sanity_help_requested=true')
lines.append('    break')
lines.append('  fi')
lines.append('done')
lines.append('')
lines.append('for arg in "${OXICHECKER_ARGS[@]}"; do')
lines.append('  if [[ "$arg" == "--help" || "$arg" == "-h" ]]; then')
lines.append('    oxichecker_help_requested=true')
lines.append('    break')
lines.append('  fi')
lines.append('done')
lines.append('')
lines.append('for arg in "${POST_ARGS[@]}"; do')
lines.append('  if [[ "$arg" == "--help" || "$arg" == "-h" ]]; then')
lines.append('    post_help_requested=true')
lines.append('    break')
lines.append('  fi')
lines.append('done')
lines.append('')
lines.append('# Handle help requests - show combined help instead of individual help messages')
lines.append('if [[ "$sanity_help_requested" == true || "$oxichecker_help_requested" == true || "$post_help_requested" == true ]]; then')
lines.append('  echo "Individual help messages are not available. Use entrypoint --help for available options."')
lines.append('  exit 0')
lines.append('fi')
lines.append('')
lines.append('# Function to find and use P1 CIF files if available')
lines.append('find_p1_cifs() {')
lines.append('  local input_path="$1"')
lines.append('  local output_dir="$2"  # The sanity runner output directory')
lines.append('  local p1_cifs_found=0')
lines.append('  ')
lines.append('  # Convert to absolute path')
lines.append('  local abs_input_path=$(realpath "$input_path" 2>/dev/null || echo "$input_path")')
lines.append('  ')
lines.append('  # If output_dir is provided, look for P1 files there')
lines.append('  if [[ -n "$output_dir" && -d "$output_dir/p1_cifs" ]]; then')
lines.append('    local p1_dir="$output_dir/p1_cifs"')
lines.append('    local p1_files=("$p1_dir"/P1_*.cif)')
lines.append('    if [[ -f "${p1_files[0]}" ]]; then')
lines.append('      p1_cifs_found=${#p1_files[@]}')
lines.append('      echo "Using P1 CIF files from sanity output: $p1_dir" >&2')
lines.append('      echo "Found $p1_cifs_found P1 CIF files" >&2')
lines.append('      echo "$p1_dir"')
lines.append('      return 0')
lines.append('    fi')
lines.append('  fi')
lines.append('  ')
lines.append('  # Fallback: check if input has P1 files (for manual runs)')
lines.append('  if [[ -f "$abs_input_path" ]]; then')
lines.append('    local basename=$(basename "$abs_input_path" .cif)')
lines.append('    local input_dir=$(dirname "$abs_input_path")')
lines.append('    local p1_file="$input_dir/p1_cifs/P1_${basename}.cif"')
lines.append('    ')
lines.append('    if [[ -f "$p1_file" ]]; then')
lines.append('      echo "Using P1 CIF file: $p1_file" >&2')
lines.append('      echo "$p1_file"')
lines.append('      return 0')
lines.append('    else')
lines.append('      echo "No P1 CIF file found, using original: $abs_input_path" >&2')
lines.append('      echo "$abs_input_path"')
lines.append('      return 0')
lines.append('    fi')
lines.append('  elif [[ -d "$abs_input_path" ]]; then')
lines.append('    local p1_dir="$abs_input_path/p1_cifs"')
lines.append('    ')
lines.append('    if [[ -d "$p1_dir" ]]; then')
lines.append('      local p1_files=("$p1_dir"/P1_*.cif)')
lines.append('      if [[ -f "${p1_files[0]}" ]]; then')
lines.append('        p1_cifs_found=${#p1_files[@]}')
lines.append('        echo "Using P1 CIF files from input directory: $p1_dir" >&2')
lines.append('        echo "Found $p1_cifs_found P1 CIF files" >&2')
lines.append('        echo "$p1_dir"')
lines.append('        return 0')
lines.append('      fi')
lines.append('    fi')
lines.append('    ')
lines.append('    local original_cifs=$(find "$abs_input_path" -maxdepth 1 -name "*.cif" | wc -l)')
lines.append('    echo "No P1 CIF files found, using original directory: $abs_input_path" >&2')
lines.append('    echo "Found $original_cifs original CIF files" >&2')
lines.append('    echo "$abs_input_path"')
lines.append('    return 0')
lines.append('  else')
lines.append('    echo "Input path does not exist: $abs_input_path" >&2')
lines.append('    echo "$abs_input_path"')
lines.append('    return 1')
lines.append('  fi')
lines.append('}')
lines.append('')
lines.append('# Execute based on run mode')
lines.append('case "$run_mode" in')
lines.append('  "sanity")')
lines.append('    echo "=== Running Sanity Runner Only ==="')
lines.append('    if ! "$SANITY_BIN" "${SANITY_ARGS[@]}" "${SANITY_POSITIONALS[@]}"; then')
lines.append('      exit $?')
lines.append('    fi')
lines.append('    ;;')
lines.append('  "oxichecker")')
lines.append('    echo "=== Running OxiChecker Only ==="')
lines.append('    echo "=== OxiChecker CIF Selection ==="')
lines.append('    # For oxichecker-only mode, we don\'t have output_dir from sanity, so pass empty')
lines.append('    oxichecker_input_path=$(find_p1_cifs "$oxichecker_input" "")')
lines.append('    echo "==========================="')
lines.append('    ')
lines.append('    # Update OXICHECKER_ARGS with the selected input path')
lines.append('    for i in "${!OXICHECKER_ARGS[@]}"; do')
lines.append('      if [[ "${OXICHECKER_ARGS[$i]}" == "--input" ]]; then')
lines.append('        OXICHECKER_ARGS[$((i+1))]="$oxichecker_input_path"')
lines.append('        break')
lines.append('      fi')
lines.append('    done')
lines.append('    ')
lines.append('    # Set default output location if not specified')
lines.append('    output_specified=false')
lines.append('    for i in "${!OXICHECKER_ARGS[@]}"; do')
lines.append('      if [[ "${OXICHECKER_ARGS[$i]}" == "--output" ]]; then')
lines.append('        output_specified=true')
lines.append('        break')
lines.append('      fi')
lines.append('    done')
lines.append('    ')
lines.append('    if [[ "$output_specified" == false ]]; then')
lines.append('      OXICHECKER_ARGS+=("--output" "oxichecker_results.csv")')
lines.append('    fi')
lines.append('    ')
lines.append('    exec "$PY_BIN" -W ignore::SyntaxWarning "$OXI_SCRIPT" "${OXICHECKER_ARGS[@]}"')
lines.append('    ;;')
lines.append('  "postprocess")')
lines.append('    echo "=== Running Postprocessing Only ==="')
lines.append('    exec "$PY_BIN" "$POST_SCRIPT" "${POST_ARGS[@]}"')
lines.append('    ;;')
lines.append('  "all")')
lines.append('    echo "=== Running Full Pipeline ==="')
lines.append('    ')
lines.append('    # Step 1: Sanity Runner')
lines.append('    echo "=== Step 1: Sanity Runner ==="')
lines.append('    if ! "$SANITY_BIN" "${SANITY_ARGS[@]}" "${SANITY_POSITIONALS[@]}"; then')
lines.append('      exit $?')
lines.append('    fi')
lines.append('    ')
lines.append('    # Step 2: OxiChecker')
lines.append('    echo "=== Step 2: OxiChecker ==="')
lines.append('    echo "=== OxiChecker CIF Selection ==="')
lines.append('    output_dir=""')
lines.append('    if [[ -n "$output_dir_value" ]]; then')
lines.append('      output_dir="$output_dir_value"')
lines.append('    else')
lines.append('      # Check if --output-dir was explicitly set in SANITY_ARGS')
lines.append('      for i in "${!SANITY_ARGS[@]}"; do')
lines.append('        if [[ "${SANITY_ARGS[$i]}" == "--output-dir" ]]; then')
lines.append('          output_dir="${SANITY_ARGS[$((i+1))]}"')
lines.append('          break')
lines.append('        fi')
lines.append('      done')
lines.append('      ')
lines.append('      # If still empty, calculate the default output directory like sanity_runner does')
lines.append('      if [[ -z "$output_dir" ]]; then')
lines.append('        if [[ -f "$oxichecker_input" ]]; then')
lines.append('          # Single file: use stem + _sanity_results')
lines.append('          input_stem=$(basename "$oxichecker_input" .cif)')
lines.append('          output_dir="${input_stem}_sanity_results"')
lines.append('        elif [[ -d "$oxichecker_input" ]]; then')
lines.append('          # Directory: use dirname + _sanity_results')
lines.append('          input_name=$(basename "$oxichecker_input")')
lines.append('          output_dir="${input_name}_sanity_results"')
lines.append('        else')
lines.append('          output_dir="input_sanity_results"')
lines.append('        fi')
lines.append('        # Make output_dir absolute if it\'s relative')
lines.append('        if [[ ! "$output_dir" = /* ]]; then')
lines.append('          output_dir="$(pwd)/$output_dir"')
lines.append('        fi')
lines.append('      fi')
lines.append('    fi')
lines.append('    ')
lines.append('    # Use find_p1_cifs function with the determined output directory')
lines.append('    oxichecker_input_path=$(find_p1_cifs "$oxichecker_input" "$output_dir")')
lines.append('    echo "==========================="')
lines.append('    ')
lines.append('    # Update OXICHECKER_ARGS with the selected input path')
lines.append('    for i in "${!OXICHECKER_ARGS[@]}"; do')
lines.append('      if [[ "${OXICHECKER_ARGS[$i]}" == "--input" ]]; then')
lines.append('        OXICHECKER_ARGS[$((i+1))]="$oxichecker_input_path"')
lines.append('        break')
lines.append('      fi')
lines.append('    done')
lines.append('    ')
lines.append('    # Ensure OxiChecker output goes to the sanity output directory')
lines.append('    output_specified=false')
lines.append('    for i in "${!OXICHECKER_ARGS[@]}"; do')
lines.append('      if [[ "${OXICHECKER_ARGS[$i]}" == "--output" ]]; then')
lines.append('        output_specified=true')
lines.append('        break')
lines.append('      fi')
lines.append('    done')
lines.append('    ')
lines.append('    if [[ "$output_specified" == false && -n "$output_dir" ]]; then')
lines.append('      OXICHECKER_ARGS+=("--output" "$output_dir/oxichecker_results.csv")')
lines.append('    elif [[ "$output_specified" == false ]]; then')
lines.append('      OXICHECKER_ARGS+=("--output" "oxichecker_results.csv")')
lines.append('    fi')
lines.append('    ')
lines.append('    # Set decomposition directory to be inside sanity output folder')
lines.append('    decompose_dir_specified=false')
lines.append('    for i in "${!OXICHECKER_ARGS[@]}"; do')
lines.append('      if [[ "${OXICHECKER_ARGS[$i]}" == "--decompose-dir" ]]; then')
lines.append('        decompose_dir_specified=true')
lines.append('        break')
lines.append('      fi')
lines.append('    done')
lines.append('    ')
lines.append('    if [[ "$decompose_dir_specified" == false && -n "$output_dir" ]]; then')
lines.append('      OXICHECKER_ARGS+=("--decompose-dir" "$output_dir/decompose_oxichecker_results")')
lines.append('    fi')
lines.append('    ')
lines.append('    if ! "$PY_BIN" -W ignore::SyntaxWarning "$OXI_SCRIPT" "${OXICHECKER_ARGS[@]}"; then')
lines.append('      exit $?')
lines.append('    fi')
lines.append('    ')
lines.append('    # Step 3: Postprocessing with OxiChecker results')
lines.append('    echo "=== Step 3: Postprocessing ==="')
lines.append('    ')
lines.append('    # Add OxiChecker results path if the file was created')
lines.append('    if [[ -n "$output_dir" && -f "$output_dir/oxichecker_results.csv" ]]; then')
lines.append('      # Check if oxichecker-path is already in POST_ARGS')
lines.append('      oxichecker_path_specified=false')
lines.append('      for i in "${!POST_ARGS[@]}"; do')
lines.append('        if [[ "${POST_ARGS[$i]}" == "--oxichecker-path" ]]; then')
lines.append('          oxichecker_path_specified=true')
lines.append('          break')
lines.append('        fi')
lines.append('      done')
lines.append('      ')
lines.append('      if [[ "$oxichecker_path_specified" == false ]]; then')
lines.append('        POST_ARGS+=("--oxichecker-path" "$output_dir/oxichecker_results.csv")')
lines.append('      fi')
lines.append('    fi')
lines.append('    ')
lines.append('    exec "$PY_BIN" "$POST_SCRIPT" "${POST_ARGS[@]}"')
lines.append('    ;;')
lines.append('  *)')
lines.append('    echo "Invalid run mode: $run_mode" >&2')
lines.append('    exit 2')
lines.append('    ;;')
lines.append('esac')

script = "\n".join(lines) + "\n"
ENTRYPOINT_PATH.write_text(script)
ENTRYPOINT_PATH.chmod(0o755)
PY
