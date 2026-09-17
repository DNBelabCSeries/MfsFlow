import os
import tempfile
import unittest

from mfsflow.scripts.dge_utils import (
    close_pass1_store,
    dynamic_chunksize,
    finalize_pass1_store,
    load_pass1_global_counts,
    load_pass1_read_bundle,
    load_pass1_read_counts,
    load_pass1_umi_bundle,
    load_pass1_umi_counts,
    open_pass1_store,
    pass1_barcode_workloads,
    pass1_barcodes,
    resolve_worker_count,
    workload_order,
    store_pass1_result,
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

    def test_pass1_store_round_trips_chunk_results(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            connection, path = open_pass1_store(tmpdir, "sample")
            try:
                store_pass1_result(
                    connection,
                    {
                        "exon": {"BC1": {"G1": 3}},
                        "intron": {"BC1": {"G1": 1}},
                    },
                    {
                        "exon": {"BC1": {"G1": {"AAAA": 2}}},
                        "intron": {"BC1": {"G1": {"AAAT": 1}}},
                    },
                    {"BC1": {"AAAA": 2, "AAAT": 1}},
                )
                finalize_pass1_store(connection)

                self.assertEqual(pass1_barcodes(connection, "read"), ["BC1"])
                self.assertEqual(load_pass1_read_counts(connection, "BC1", "exon"), {"G1": 3})
                self.assertEqual(
                    load_pass1_read_bundle(connection, "BC1"),
                    {"exon": {"G1": 3}, "intron": {"G1": 1}},
                )
                self.assertEqual(
                    dict(load_pass1_umi_counts(connection, "BC1", "exon")["G1"]),
                    {"AAAA": 2},
                )
                umi_bundle, global_counts = load_pass1_umi_bundle(connection, "BC1")
                self.assertEqual(dict(umi_bundle["exon"]["G1"]), {"AAAA": 2})
                self.assertEqual(dict(umi_bundle["intron"]["G1"]), {"AAAT": 1})
                self.assertEqual(dict(global_counts), {"AAAA": 2, "AAAT": 1})
                self.assertEqual(
                    dict(load_pass1_global_counts(connection, "BC1")),
                    {"AAAA": 2, "AAAT": 1},
                )
                self.assertEqual(pass1_barcode_workloads(connection, True), [("BC1", 4)])
            finally:
                close_pass1_store(connection, path)
            self.assertFalse(os.path.exists(path))


if __name__ == "__main__":
    unittest.main()
