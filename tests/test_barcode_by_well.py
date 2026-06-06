import os
import sys
import tempfile
import unittest

try:
    import pandas as pd
except ImportError:
    pd = None

if pd is not None:
    from mfsflow.scripts.barcode_detection import write_barcodes_by_well
else:
    write_barcodes_by_well = None


class BarcodeByWellTests(unittest.TestCase):
    @unittest.skipIf(pd is None, "pandas is not installed")
    def test_write_barcodes_by_well_combines_umi_and_internal_counts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            expect = os.path.join(tmpdir, "expect_id_barcode.tsv")
            with open(expect, "w") as f:
                f.write("wellID\tumi_barcodes\tinternal_barcodes\n")
                f.write("P1A1\tAAAA,CCCC\tGGGG,TTTT\n")

            kept_df = pd.DataFrame({
                "XC": ["AAAA", "CCCC", "GGGG", "TTTT", "NNNN"],
                "n": [10, 5, 3, 2, 99],
            })
            out_file = os.path.join(tmpdir, "by_well.tsv")

            self.assertTrue(write_barcodes_by_well(kept_df, expect, out_file))
            with open(out_file) as f:
                lines = [line.strip().split("\t") for line in f]

            self.assertEqual(lines[0][:4], ["wellID", "umi_reads", "internal_reads", "total_reads"])
            self.assertEqual(lines[1][:4], ["P1A1", "15", "5", "20"])


if __name__ == "__main__":
    unittest.main()
