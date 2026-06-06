import os
import sys
import unittest

from mfsflow.scripts.umi_utils import cluster_umis


class UmiUtilsTests(unittest.TestCase):
    def test_threshold_zero_keeps_all_umis(self):
        self.assertEqual(
            cluster_umis({"AAAA": 5, "AAAT": 1}, threshold=0),
            {"AAAA": "AAAA", "AAAT": "AAAT"},
        )

    def test_threshold_one_matches_existing_behavior(self):
        self.assertEqual(
            cluster_umis({"AAAA": 10, "AAAT": 1, "AATT": 1, "TTTT": 1}, threshold=1),
            {"AAAA": "AAAA", "AAAT": "AAAA", "AATT": "AATT", "TTTT": "TTTT"},
        )

    def test_threshold_two_adds_two_mismatch_neighbors(self):
        self.assertEqual(
            cluster_umis({"AAAA": 10, "AAAT": 1, "AATT": 1, "TTTT": 1}, threshold=2),
            {"AAAA": "AAAA", "AAAT": "AAAA", "AATT": "AAAA", "TTTT": "TTTT"},
        )


if __name__ == "__main__":
    unittest.main()
