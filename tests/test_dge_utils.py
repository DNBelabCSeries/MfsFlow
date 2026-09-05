import unittest

from mfsflow.scripts.dge_utils import (
    dynamic_chunksize,
    resolve_worker_count,
    workload_order,
)


class DgeUtilsTests(unittest.TestCase):
    def test_worker_count_respects_task_count_and_optional_cap(self):
        self.assertEqual(resolve_worker_count(20, 3), 3)
        self.assertEqual(resolve_worker_count(20, 30, {"max_dge_workers": 6}), 6)

    def test_worker_count_rejects_invalid_cap(self):
        with self.assertRaises(ValueError):
            resolve_worker_count(4, 10, {"max_dge_workers": 0})

    def test_workload_order_is_heaviest_first_and_deterministic(self):
        self.assertEqual(
            workload_order([("B", 2), ("A", 2), ("C", 5)]),
            ["C", "A", "B"],
        )

    def test_dynamic_chunksize_leaves_work_for_stealing(self):
        self.assertEqual(dynamic_chunksize(100, 4), 7)
        self.assertEqual(dynamic_chunksize(2, 8), 1)


if __name__ == "__main__":
    unittest.main()
