#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
"""Compare DUMP_*.xyz outputs of cif_to_building_blocks_v2 and v3."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

DEFAULT_V2 = Path("build/examples/cif_to_building_blocks_v2")
DEFAULT_V3 = Path("build/examples/cif_to_building_blocks_v3")
DEFAULT_DATA_DIR = Path("tests/qmof_cifs_for_v2_v3_comparison")
DEFAULT_PARAMS = (0.5, 5.0, -0.45, "jmol")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run both building-block decomposers and compare their DUMP_*.xyz outputs. "
            "Differences are reported ignoring the ordering of conformations."
        ),
    )
    parser.add_argument(
        "cifs",
        metavar="CIF",
        nargs="*",
        help="Paths to CIF files. Defaults to all .cif files in tests/qmof_cifs_for_v2_v3_comparison.",
    )
    parser.add_argument(
        "--v2",
        type=Path,
        default=DEFAULT_V2,
        help=f"Path to cif_to_building_blocks_v2 (default: {DEFAULT_V2})",
    )
    parser.add_argument(
        "--v3",
        type=Path,
        default=DEFAULT_V3,
        help=f"Path to cif_to_building_blocks_v3 (default: {DEFAULT_V3})",
    )
    parser.add_argument("--r-min", type=float, default=DEFAULT_PARAMS[0])
    parser.add_argument("--r-max", type=float, default=DEFAULT_PARAMS[1])
    parser.add_argument("--tolerance", type=float, default=DEFAULT_PARAMS[2])
    parser.add_argument("--strategy", default=DEFAULT_PARAMS[3])
    parser.add_argument(
        "--keep",
        action="store_true",
        help="Keep per-run outputs under build/building_block_comparisons/",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print stdout/stderr from the executables when runs succeed.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1, use 0 for auto).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Write detailed comparison log to this path.",
    )
    return parser.parse_args()


def resolve_cifs(args: argparse.Namespace) -> List[Path]:
    if args.cifs:
        return [Path(cif).resolve() for cif in args.cifs]
    return sorted((DEFAULT_DATA_DIR).resolve().glob("*.cif"))


@dataclass(frozen=True)
class RunConfig:
    v2: Path
    v3: Path
    r_min: float
    r_max: float
    tolerance: float
    strategy: str
    keep: bool
    verbose: bool
    keep_root: Path


@dataclass(frozen=True)
class CompareTask:
    index: int
    cif: Path


@dataclass
class CompareOutcome:
    index: int
    cif: Path
    status: str  # "match", "difference", "error"
    lines: List[str]


def run_decomposer(exe: Path, cif: Path, out_dir: Path, config: RunConfig) -> subprocess.CompletedProcess[str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(exe.resolve()),
        str(cif),
        str(config.r_min),
        str(config.r_max),
        str(config.tolerance),
        config.strategy,
    ]
    return subprocess.run(
        cmd,
        cwd=out_dir,
        capture_output=True,
        text=True,
        check=False,
    )


def normalise_line(line: str) -> str:
    return " ".join(line.split())


def split_frames(content: str) -> Iterable[str]:
    chunks = re.split(r"\n\s*\n", content.strip())
    for chunk in chunks:
        cleaned = chunk.strip()
        if cleaned:
            yield "\n".join(normalise_line(line) for line in cleaned.splitlines())


@dataclass
class DumpData:
    counter: Counter[str]


def collect_dump_data(out_dir: Path) -> Dict[str, DumpData]:
    dumps: Dict[str, DumpData] = {}
    for path in sorted(out_dir.rglob("DUMP_*.xyz")):
        rel = str(path.relative_to(out_dir))
        text = path.read_text()
        frames = Counter(split_frames(text))
        dumps[rel] = DumpData(counter=frames)
    return dumps


def format_counter(counter: Counter[str]) -> List[str]:
    parts = []
    for frame, count in counter.most_common():
        lines = frame.splitlines()
        comment = lines[1] if len(lines) > 1 else lines[0]
        digest = hashlib.sha256(frame.encode()).hexdigest()[:16]
        parts.append(f"count={count} digest={digest} comment={comment[:60]}")
    return parts


def compare_dump_sets(v2: Dict[str, DumpData], v3: Dict[str, DumpData]) -> Dict[str, Dict[str, Counter[str]]]:
    differences: Dict[str, Dict[str, Counter[str]]] = {}
    keys = sorted(set(v2) | set(v3))
    for key in keys:
        if key not in v2:
            differences[key] = {
                "missing_in_v2": Counter(),
                "missing_in_v3": Counter(v3[key].counter),
            }
        elif key not in v3:
            differences[key] = {
                "missing_in_v2": Counter(v2[key].counter),
                "missing_in_v3": Counter(),
            }
        else:
            extra_v2 = v2[key].counter - v3[key].counter
            extra_v3 = v3[key].counter - v2[key].counter
            if extra_v2 or extra_v3:
                differences[key] = {
                    "missing_in_v2": extra_v3,
                    "missing_in_v3": extra_v2,
                }
    return differences


def render_differences(
    cif: Path,
    differences: Dict[str, Dict[str, Counter[str]]],
    v2_proc: subprocess.CompletedProcess[str],
    v3_proc: subprocess.CompletedProcess[str],
    verbose: bool,
) -> List[str]:
    lines = [f"{cif.name}: differences detected in DUMP files"]
    for dump_name, diff in differences.items():
        lines.append(f"  {dump_name}:")
        only_v3 = diff.get("missing_in_v2", Counter())
        if only_v3:
            lines.append("    frames present only in v3:")
            for line in format_counter(only_v3):
                lines.append(f"      {line}")
        only_v2 = diff.get("missing_in_v3", Counter())
        if only_v2:
            lines.append("    frames present only in v2:")
            for line in format_counter(only_v2):
                lines.append(f"      {line}")
    if verbose:
        lines.append("  v2 stdout:\n" + v2_proc.stdout)
        lines.append("  v3 stdout:\n" + v3_proc.stdout)
        if v2_proc.stderr:
            lines.append("  v2 stderr:\n" + v2_proc.stderr)
        if v3_proc.stderr:
            lines.append("  v3 stderr:\n" + v3_proc.stderr)
    return lines


def _prepare_output_dirs(task: CompareTask, config: RunConfig) -> tuple[Path, Path]:
    if config.keep:
        base = config.keep_root / task.cif.stem
        v2_dir = base / "v2"
        v3_dir = base / "v3"
        if v2_dir.exists():
            shutil.rmtree(v2_dir)
        if v3_dir.exists():
            shutil.rmtree(v3_dir)
        v2_dir.mkdir(parents=True, exist_ok=True)
        v3_dir.mkdir(parents=True, exist_ok=True)
        return v2_dir, v3_dir

    temp_dir = tempfile.mkdtemp(prefix=f"compare_{task.cif.stem}_", dir=str(config.keep_root))
    base = Path(temp_dir)
    v2_dir = base / "v2"
    v3_dir = base / "v3"
    v2_dir.mkdir(parents=True, exist_ok=True)
    v3_dir.mkdir(parents=True, exist_ok=True)
    return v2_dir, v3_dir


def _cleanup_dirs(v2_dir: Path, v3_dir: Path, config: RunConfig) -> None:
    if config.keep:
        return
    shutil.rmtree(v2_dir.parent, ignore_errors=True)


def process_cif(task: CompareTask, config: RunConfig) -> CompareOutcome:
    v2_dir: Optional[Path] = None
    v3_dir: Optional[Path] = None
    v2_proc: Optional[subprocess.CompletedProcess[str]] = None
    v3_proc: Optional[subprocess.CompletedProcess[str]] = None

    try:
        v2_dir, v3_dir = _prepare_output_dirs(task, config)
        v2_proc = run_decomposer(config.v2, task.cif, v2_dir, config)
        v3_proc = run_decomposer(config.v3, task.cif, v3_dir, config)

        if v2_proc.returncode != 0 or v3_proc.returncode != 0:
            lines = [
                (
                    f"{task.cif.name}: execution failure ("
                    f"v2 code {v2_proc.returncode}, v3 code {v3_proc.returncode})"
                )
            ]
            if config.verbose:
                if v2_proc.stdout:
                    lines.append("  v2 stdout:\n" + v2_proc.stdout)
                if v3_proc.stdout:
                    lines.append("  v3 stdout:\n" + v3_proc.stdout)
                if v2_proc.stderr:
                    lines.append("  v2 stderr:\n" + v2_proc.stderr)
                if v3_proc.stderr:
                    lines.append("  v3 stderr:\n" + v3_proc.stderr)
            return CompareOutcome(task.index, task.cif, "error", lines)

        dump_v2 = collect_dump_data(v2_dir)
        dump_v3 = collect_dump_data(v3_dir)

        differences = compare_dump_sets(dump_v2, dump_v3)
        if differences:
            lines = render_differences(task.cif, differences, v2_proc, v3_proc, config.verbose)
            return CompareOutcome(task.index, task.cif, "difference", lines)

        lines = [f"{task.cif.name}: DUMP outputs match"]
        if config.verbose:
            lines.append("  stdout:\n" + v2_proc.stdout)
            if v2_proc.stderr:
                lines.append("  stderr:\n" + v2_proc.stderr)
        return CompareOutcome(task.index, task.cif, "match", lines)

    except Exception as exc:  # noqa: BLE001
        message = f"{task.cif.name}: exception during comparison: {exc}"
        return CompareOutcome(task.index, task.cif, "error", [message])
    finally:
        if v2_dir is not None and v3_dir is not None:
            _cleanup_dirs(v2_dir, v3_dir, config)


def display_progress(done: int, total: int, matched: int, mismatched: int, errors: int) -> None:
    if total == 0:
        return

    bar_length = 30
    filled = int(bar_length * done / total)
    bar = "#" * filled + "-" * (bar_length - filled)

    matched_pct = (matched / total) * 100 if total else 0.0
    mismatched_pct = (mismatched / total) * 100 if total else 0.0

    parts = [
        f"[{bar}] {done}/{total}",
        f"matched: {matched} ({matched_pct:.1f}%)",
        f"mismatched: {mismatched} ({mismatched_pct:.1f}%)",
    ]
    if errors:
        parts.append(f"errors: {errors}")

    message = " | ".join(parts)
    sys.stderr.write("\r" + message.ljust(120))
    sys.stderr.flush()
    if done == total:
        sys.stderr.write("\n")
        sys.stderr.flush()


def main() -> int:
    args = parse_args()
    v2_exe = args.v2.resolve()
    v3_exe = args.v3.resolve()
    if not v2_exe.is_file():
        print(f"Missing v2 executable: {v2_exe}", file=sys.stderr)
        return 2
    if not v3_exe.is_file():
        print(f"Missing v3 executable: {v3_exe}", file=sys.stderr)
        return 2

    cif_paths = resolve_cifs(args)
    if not cif_paths:
        print("No CIF files to test", file=sys.stderr)
        return 2

    keep_root = Path("build/building_block_comparisons").resolve()
    keep_root.mkdir(parents=True, exist_ok=True)

    jobs = args.n_jobs
    if jobs <= 0:
        jobs = os.cpu_count() or 1

    output_path: Optional[Path] = None
    if args.output is not None:
        output_path = args.output.resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)

    config = RunConfig(
        v2=v2_exe,
        v3=v3_exe,
        r_min=args.r_min,
        r_max=args.r_max,
        tolerance=args.tolerance,
        strategy=args.strategy,
        keep=args.keep,
        verbose=args.verbose,
        keep_root=keep_root,
    )

    tasks = [CompareTask(index=i, cif=path) for i, path in enumerate(cif_paths)]
    total = len(tasks)

    matched = 0
    mismatches = 0
    errors = 0
    done = 0
    next_to_emit = 0
    ordered_results: List[Optional[CompareOutcome]] = [None] * total
    collected_lines: Optional[List[str]] = [] if output_path is not None else None

    def emit(outcome: CompareOutcome) -> None:
        nonlocal matched, mismatches, errors, done, next_to_emit
        if outcome.status == "match":
            matched += 1
        elif outcome.status == "difference":
            mismatches += 1
        else:
            errors += 1
        done += 1
        ordered_results[outcome.index] = outcome

        while next_to_emit < total and ordered_results[next_to_emit] is not None:
            current = ordered_results[next_to_emit]
            assert current is not None
            if output_path is None:
                for line in current.lines:
                    print(line)
            else:
                assert collected_lines is not None
                collected_lines.extend(current.lines)
            next_to_emit += 1

        display_progress(done, total, matched, mismatches, errors)

    if jobs == 1:
        for task in tasks:
            outcome = process_cif(task, config)
            emit(outcome)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as executor:
            future_to_task = {executor.submit(process_cif, task, config): task for task in tasks}
            for future in concurrent.futures.as_completed(future_to_task):
                task = future_to_task[future]
                try:
                    outcome = future.result()
                except Exception as exc:  # noqa: BLE001
                    message = f"{task.cif.name}: worker error: {exc}"
                    outcome = CompareOutcome(task.index, task.cif, "error", [message])
                emit(outcome)

    if output_path is not None and collected_lines is not None:
        lines_to_write = collected_lines.copy()
        while lines_to_write and not lines_to_write[-1]:
            lines_to_write.pop()
        text = "\n".join(lines_to_write)
        if text and not text.endswith("\n"):
            text += "\n"
        output_path.write_text(text)

    return 0 if mismatches == 0 and errors == 0 else 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
