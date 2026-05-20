#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
import os
import subprocess
import time
import sys
import resource
from multiprocessing import Pool, Manager, Process, current_process
from pathlib import Path
import click

# Global log queue (to be initialized in each process)
log_queue = None

def init_worker(queue):
    global log_queue
    log_queue = queue

def set_memory_limit(max_mem_gb):
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    limit_bytes = int(max_mem_gb * (1024 ** 3))
    resource.setrlimit(resource.RLIMIT_AS, (limit_bytes, hard))
    
    # Disable core dumps to prevent core.NNNN files
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))

def process_cif(args):
    binary_path, cif_path, params, max_mem_gb, max_time = args
    global log_queue

    cif_path = Path(cif_path)
    cif_name = cif_path.stem
    work_dir = cif_path.parent / cif_name
    work_dir.mkdir(exist_ok=True)

    target_cif = work_dir / cif_path.name
    if cif_path.resolve() != target_cif.resolve():
        target_cif.write_bytes(cif_path.read_bytes())

    cmd = [str(Path(binary_path).resolve()), str(target_cif.name)] + list(params)

    t0 = time.time()
    try:
        set_memory_limit(max_mem_gb)

        subprocess.run(
            cmd,
            cwd=work_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True,
            timeout=max_time
        )
        elapsed = time.time() - t0
        log_queue.put(f"[SUCCESS] {cif_name}. Done in {elapsed:.2f} s")

    except subprocess.TimeoutExpired:
        log_queue.put(f"[FAILED]  {cif_name}: Time limit exceeded")
    except MemoryError:
        log_queue.put(f"[FAILED]  {cif_name}: Memory limit exceeded")
    except subprocess.CalledProcessError as e:
        msg = e.stderr
        if "memory" in msg.lower() or "alloc" in msg.lower() or "killed" in msg.lower():
            log_queue.put(f"[FAILED]  {cif_name}: Memory limit exceeded")
        else:
            error_text = "; ".join(str(msg).split("\n"))
            log_queue.put(f"[FAILED]  {cif_name}: {error_text}")
    except Exception as e:
        error_text = "; ".join(str(e).split("\n"))
        log_queue.put(f"[FAILED]  {cif_name}: {error_text}")

def logger_worker(log_queue, log_file):
    output = open(log_file, 'a', buffering=1) if log_file else sys.stdout
    try:
        while True:
            msg = log_queue.get()
            if msg == "__STOP__":
                break
            print(msg, file=output, flush=True)
    finally:
        if log_file:
            output.close()

@click.command()
@click.option('--binary-path', required=True, type=click.Path(exists=True, dir_okay=False), help="Path to the binary executable.")
@click.option('--source-dir', default='.', type=click.Path(exists=True, file_okay=False), help="Directory with .cif files.")
@click.option('--max-workers', default=64, show_default=True, help="Number of parallel workers.")
@click.option('--params', multiple=True, default=("0.5", "5.0", "-0.45", "jmol"), show_default=True, help="Parameters for the binary.")
@click.option('--max-memory', default=8, show_default=True, type=float, help="Max memory per process (GB).")
@click.option('--max-time', default=15, show_default=True, type=float, help="Max wall time per process (seconds).")
@click.option('--log-file', type=click.Path(dir_okay=False, writable=True), default=None, help="Optional path to log file.")
def main(binary_path, source_dir, max_workers, params, max_memory, max_time, log_file):
    cif_files = list(Path(source_dir).glob("*.cif"))
    print(f"Found {len(cif_files)} CIF files.")
    
    with Manager() as manager:
        queue = manager.Queue()
        log_proc = Process(target=logger_worker, args=(queue, log_file))
        log_proc.start()

        # Send initial log line
        queue.put(f"Found {len(cif_files)} CIF files.")

        args_list = [
            (binary_path, cif, params, max_memory, max_time)
            for cif in cif_files
        ]

        with Pool(processes=max_workers, initializer=init_worker, initargs=(queue,)) as pool:
            pool.map(process_cif, args_list)

        queue.put("__STOP__")
        log_proc.join()

if __name__ == "__main__":
    main()
