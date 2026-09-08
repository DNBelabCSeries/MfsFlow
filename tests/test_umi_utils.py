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

    def test_threshold_three_is_consistent_across_dispatch_boundary(self):
        parent = "AAAAAAAAAA"
        child = "TTTAAAAAAA"
        small = {parent: 10, child: 1}

        large = dict(small)
        for index in range(98):
            suffix = format(index, "04b")
            large["CCCCCC" + suffix] = 1

        self.assertEqual(cluster_umis(small, threshold=3)[child], parent)
        self.assertEqual(cluster_umis(large, threshold=3)[child], parent)


if __name__ == "__main__":
    unittest.main()
