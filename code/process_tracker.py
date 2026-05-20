# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.
import threading

# Global tracking for child processes that need to be killed on timeout
_active_processes = []
_process_lock = threading.Lock()