import os
import tempfile
import unittest

try:
    import pandas as pd
    from mfsflow.scripts.barcode_detection import cell_bc_selection
except ImportError:
    pd = None
    cell_bc_selection = None


class BarcodeSelectionTests(unittest.TestCase):
    @unittest.skipIf(pd is None, "pandas is not installed")
    def test_known_whitelist_does_not_fallback_to_top_barcodes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            whitelist = os.path.join(tmpdir, "expect_barcode.tsv")
            with open(whitelist, "w") as handle:
                handle.write("CCCCCCCCCCCCCCCCCCCC\n")

            counts = pd.DataFrame({
                "XC": ["AAAAAAAAAAAAAAAAAAAA", "GGGGGGGGGGGGGGGGGGGG"],
                "n": [100, 90],
            })
            config = {
                "barcodes": {
                    "nReadsperCell": 1,
                    "automatic": False,
                    "barcode_file": whitelist,
                }
            }

            with self.assertRaisesRegex(ValueError, "None of the annotated barcodes"):
                cell_bc_selection(counts, config)

    @unittest.skipIf(pd is None, "pandas is not installed")
    def test_selection_keeps_low_count_rows_available_for_downstream_binning(self):
        counts = pd.DataFrame({
            "XC": ["AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAT", "GGGGGGGGGGGGGGGGGGGG"],
            "n": [100, 2, 1],
        })
        config = {
            "barcodes": {
                "nReadsperCell": 10,
                "automatic": False,
                "barcode_num": 1,
            }
        }

        selected = cell_bc_selection(counts, config)
        self.assertEqual(
            ["AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAT", "GGGGGGGGGGGGGGGGGGGG"],
            selected["XC"].tolist(),
        )
        self.assertEqual([True, False, False], selected["keep"].tolist())


if __name__ == "__main__":
    unittest.main()
