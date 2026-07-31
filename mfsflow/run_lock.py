"""Single-writer lock for a pipeline output directory."""

import json
import os
import socket
from datetime import datetime


class OutputRunLock:
    """Prevent two MfsFlow processes from writing the same output directory."""

    def __init__(self, out_dir):
        self.path = os.path.join(out_dir, ".run.lock")
        self._acquired = False

    @staticmethod
    def _process_alive(pid):
        try:
            os.kill(int(pid), 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        except (TypeError, ValueError, OSError):
            return False
        return True

    def acquire(self):
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        payload = {
            "pid": os.getpid(),
            "host": socket.gethostname(),
            "started_at": datetime.now().isoformat(timespec="seconds"),
        }
        try:
            fd = os.open(self.path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
        except FileExistsError:
            owner = None
            try:
                with open(self.path, "r", encoding="utf-8") as handle:
                    owner = json.load(handle)
            except (OSError, ValueError):
                raise RuntimeError(
                    f"Output directory is locked by an unknown process: {self.path}. "
                    "Remove the lock only after confirming no MfsFlow process is running."
                )

            if (
                not isinstance(owner, dict)
                or not owner.get("host")
                or not owner.get("pid")
            ):
                raise RuntimeError(
                    f"Output directory is locked by an unknown process: {self.path}. "
                    "Remove the lock only after confirming no MfsFlow process is running."
                )
            same_host = owner.get("host") == payload["host"]
            owner_pid = owner.get("pid")
            if same_host and owner_pid and self._process_alive(owner_pid):
                raise RuntimeError(
                    f"Output directory is already being processed by PID {owner_pid}: {self.path}"
                )
            if not same_host:
                raise RuntimeError(
                    f"Output directory has an active or stale lock from host {owner.get('host')}: {self.path}. "
                    "Confirm the other host is finished before removing the lock."
                )

            try:
                os.remove(self.path)
            except OSError as exc:
                raise RuntimeError(f"Could not remove stale output lock {self.path}: {exc}") from exc
            return self.acquire()

        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(payload, handle)
                handle.write("\n")
        except Exception:
            try:
                os.remove(self.path)
            except OSError:
                pass
            raise
        self._acquired = True

    def release(self):
        if not self._acquired:
            return
        try:
            owner = None
            try:
                with open(self.path, "r", encoding="utf-8") as handle:
                    owner = json.load(handle)
            except (OSError, ValueError):
                owner = None
            if isinstance(owner, dict) and owner.get("pid") == os.getpid() and owner.get("host") == socket.gethostname():
                os.remove(self.path)
        finally:
            self._acquired = False

    def __enter__(self):
        self.acquire()
        return self

    def __exit__(self, _exc_type, _exc_value, _traceback):
        self.release()
