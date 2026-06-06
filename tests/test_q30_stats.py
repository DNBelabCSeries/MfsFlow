import os
import sys
import tempfile
import unittest

from mfsflow.pipeline_modules import merge_q30_stats


class Q30StatsTests(unittest.TestCase):
    def test_merge_q30_stats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_merge = os.path.join(tmpdir, "tmp_merge")
            out_dir = os.path.join(tmpdir, "XPRESS_PROCESSING")
            os.makedirs(tmp_merge)
            os.makedirs(out_dir)

            with open(os.path.join(tmp_merge, "sample.aa.Q30stats.txt"), "w") as f:
                f.write("metric\ttotal_bases\tq30_bases\n")
                f.write("R1\t10\t8\n")
                f.write("BC\t5\t5\n")
            with open(os.path.join(tmp_merge, "sample.ab.Q30stats.txt"), "w") as f:
                f.write("metric\ttotal_bases\tq30_bases\n")
                f.write("R1\t30\t22\n")
                f.write("BC\t5\t4\n")

            out_file = merge_q30_stats(tmp_merge, "sample", out_dir)
            self.assertTrue(os.path.exists(out_file))

            with open(out_file) as f:
                rows = {line.split("\t")[0]: line.strip().split("\t") for line in f if not line.startswith("metric")}

            self.assertEqual(rows["R1"], ["R1", "40", "30", "0.750000"])
            self.assertEqual(rows["BC"], ["BC", "10", "9", "0.900000"])


if __name__ == "__main__":
    unittest.main()
