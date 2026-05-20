# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
import subprocess as s
from pathlib import Path
import os
from tempfile import TemporaryDirectory
import shutil
import logging
import signal
import threading
from framechecker import FrameChecker

logger = logging.getLogger(__name__)

# Import the global process tracking
from process_tracker import _active_processes, _process_lock

PLATON_BINARY_PATH = '/code/platon'
LIBCIF_BINARY_PATH = "/opt/libcif/bin/cif_check"

def run_platon(filepath, timeout: int = 120) -> tuple[bool, list[str]]:
    filepath = str(filepath)
    logger.info("Running PLATON analyzer")
    try:
        with TemporaryDirectory() as tmpdir:
            basename = os.path.basename(filepath)
            shutil.copy(filepath, os.path.join(tmpdir, basename) )
            run = s.Popen([PLATON_BINARY_PATH, '-u', basename], cwd=tmpdir, stderr=s.PIPE,
                         stdout=s.PIPE, preexec_fn=os.setsid)
            # Register process for global timeout tracking
            with _process_lock:
                _active_processes.append(run)
            try:
                _ = run.communicate(timeout=timeout)
            finally:
                # Remove process from tracking when done
                with _process_lock:
                    if run in _active_processes:
                        _active_processes.remove(run)
            chkfilename = os.path.join(tmpdir, basename.replace(".cif", ".chk"))
            with open(chkfilename) as f:
                lines = [x.strip() for x in f.readlines()]
                last_line_id = [idx for idx,line in enumerate(lines) if "ALERT_Level and ALERT_Type Summary" in line][0]
                lines = lines[:last_line_id]
            logger.info("Stored data from PLATON report file")
        return True, lines
    except s.TimeoutExpired as e:
        logger.error(f"PLATON timeout: PLATON analysis timed out after {timeout} seconds!")
        try:
            # Kill the entire process group
            os.killpg(os.getpgid(run.pid), signal.SIGTERM)
            run.wait(timeout=5)
        except (ProcessLookupError, s.TimeoutExpired):
            try:
                os.killpg(os.getpgid(run.pid), signal.SIGKILL)
            except ProcessLookupError:
                pass
        return False, []
    except Exception as e:
        logger.error(f"PLATON analysis failed: {e}")
        return False, []
    

def run_libcif(filepath, timeout: int = 60):
    filepath = str(filepath)
    logger.info("Running LIBCIF analyzer")
    try:
        run = s.Popen([LIBCIF_BINARY_PATH, filepath, 
                    "0.5", "5.0", "1.2"], 
                    stderr=s.PIPE, stdout=s.PIPE, preexec_fn=os.setsid)
        # Register process for global timeout tracking
        with _process_lock:
            _active_processes.append(run)
        try:
            out, err = [x.decode() for x in run.communicate(timeout=timeout)]
        finally:
            # Remove process from tracking when done
            with _process_lock:
                if run in _active_processes:
                    _active_processes.remove(run)
        logger.info("Finished LIBCIF")
        return {"libcif_validity": True, 'libcif_out': out,  "libcif_error" : err}
    except s.TimeoutExpired as e:
        logger.error(f"LibCIF analysis timed out after {timeout} seconds: {e}")
        try:
            # Kill the entire process group
            os.killpg(os.getpgid(run.pid), signal.SIGTERM)
            run.wait(timeout=5)
        except (ProcessLookupError, s.TimeoutExpired):
            try:
                os.killpg(os.getpgid(run.pid), signal.SIGKILL)
            except ProcessLookupError:
                pass
        return {"libcif_validity": False, 'libcif_out': "",  "libcif_error" : f"Timeout after {timeout}s"}
    except Exception as e:
        logger.error(f"LibCif failed: {e}")
        return {"libcif_validity": False, 'libcif_out': "",  "libcif_error" : str(e)}

def run_framechecker(structure, graph, timeout: int = 180):
    logger.info("Running FrameChecker")
    try:
        checker = FrameChecker(structure, graph)
        tmp_dict = {
            "HAS_OMS" : int(checker.has_oms),
            "HAS_ATOMIC_OVERLAPS" : int(checker.has_atomic_overlaps),
            "HAS_OVERCOORDINATED_C" : int(checker.has_overcoordinated_c),
            "HAS_OVERCOORDINATED_H" : int(checker.has_overcoordinated_h),
            "HAS_OVERCOORDINATED_N" : int(checker.has_overcoordinated_n),
            "HAS_LONE_MOLECULE" : int(checker.has_lone_molecule),
            "HAS_BAD_RARE_EARTH" : int(checker.has_undercoordinated_rare_earth),
            "HAS_BAD_ALKALI_ALKALINE" : int(checker.has_undercoordinated_alkali_alkaline),
            "HAS_BAD_TERMINAL_OXO" : int(checker.has_suspicicious_terminal_oxo),
            "DECORATED_GRAPH_HASH" : checker.graph_hash,
            "UNDECORATED_GRAPH_HASH" : checker.undecorated_graph_hash,
            "DECORATED_SCAFFOLD_HASH" : checker.decorated_scaffold_hash,
            "UNDECORATED_SCAFFOLD_HASH" : checker.undecorated_scaffold_hash }
        logger.info("Finished FrameChecker")
        return {'framecheker_validity' : True, 'framecheker_results' : tmp_dict}
    except Exception as e:
        logger.error(f"FrameChecker Failed: {e}")
        return {'framecheker_validity' : False, 'framecheker_results' : None, 'framechecker_error': str(e)}
        
