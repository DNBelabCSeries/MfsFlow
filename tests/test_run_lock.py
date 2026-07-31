import json
import os
import socket
import tempfile
import unittest

from mfsflow.run_lock import OutputRunLock


class OutputRunLockTests(unittest.TestCase):
    def test_lock_lifecycle(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            lock_path = os.path.join(tmpdir, ".run.lock")
            with OutputRunLock(tmpdir):
                self.assertTrue(os.path.isfile(lock_path))
                with self.assertRaisesRegex(RuntimeError, "already being processed"):
                    OutputRunLock(tmpdir).acquire()
            self.assertFalse(os.path.exists(lock_path))

    def test_malformed_lock_is_not_removed_automatically(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            lock_path = os.path.join(tmpdir, ".run.lock")
            with open(lock_path, "w", encoding="utf-8") as handle:
                handle.write("incomplete")

            with self.assertRaisesRegex(RuntimeError, "unknown process"):
                OutputRunLock(tmpdir).acquire()
            self.assertTrue(os.path.exists(lock_path))

    def test_stale_same_host_lock_can_be_replaced(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            lock_path = os.path.join(tmpdir, ".run.lock")
            with open(lock_path, "w", encoding="utf-8") as handle:
                json.dump({"pid": 999999999, "host": socket.gethostname()}, handle)

            lock = OutputRunLock(tmpdir)
            lock.acquire()
            try:
                with open(lock_path, "r", encoding="utf-8") as handle:
                    owner = json.load(handle)
                self.assertEqual(owner["pid"], os.getpid())
            finally:
                lock.release()


if __name__ == "__main__":
    unittest.main()
