import os
import tempfile
import unittest

try:
    import pandas as pd
    from mfsflow.scripts.barcode_detection import (
        cell_bc_selection,
        write_unmatched_whitelist_barcodes,
    )
except ImportError:
    pd = None
    cell_bc_selection = None
    write_unmatched_whitelist_barcodes = None


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

    @unittest.skipIf(pd is None, "pandas is not installed")
    def test_unmatched_whitelist_file_contains_top_100_by_reads(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            counts = pd.DataFrame({
                "XC": ["WHITELIST"] + [f"BC{i:03d}" for i in range(105)],
                "n": [1000] + list(range(105)),
            })
            out_file = os.path.join(tmpdir, "sample.unmatched_whitelist_barcodes.tsv")

            written = write_unmatched_whitelist_barcodes(
                counts,
                {"WHITELIST"},
                out_file,
            )
            result = pd.read_csv(out_file, sep="\t")

            self.assertEqual(written, 100)
            self.assertEqual(result.columns.tolist(), ["barcode", "reads"])
            self.assertEqual(result.iloc[0].to_dict(), {"barcode": "BC104", "reads": 104})
            self.assertEqual(result.iloc[-1].to_dict(), {"barcode": "BC005", "reads": 5})
            self.assertNotIn("WHITELIST", set(result["barcode"]))

    @unittest.skipIf(pd is None, "pandas is not installed")
    def test_unmatched_whitelist_file_is_header_only_when_everything_matches(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            counts = pd.DataFrame({"XC": ["AAAA", "CCCC"], "n": [20, 10]})
            out_file = os.path.join(tmpdir, "unmatched.tsv")

            written = write_unmatched_whitelist_barcodes(
                counts,
                {"AAAA", "CCCC"},
                out_file,
            )

            self.assertEqual(written, 0)
            with open(out_file) as handle:
                self.assertEqual(handle.read(), "barcode\treads\n")


if __name__ == "__main__":
    unittest.main()
