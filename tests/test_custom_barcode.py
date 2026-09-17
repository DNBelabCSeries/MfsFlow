import os
import tempfile
import unittest

from mfsflow.bootstrap import create_barcode_tables
from mfsflow.path_layout import config_dir, ensure_layout

try:
    from mfsflow.scripts.barcode_detection import read_whitelist
except ImportError:
    read_whitelist = None


class CustomBarcodeTests(unittest.TestCase):
    def test_custom_table_is_validated_and_normalized(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            custom = os.path.join(tmpdir, "custom.tsv")
            with open(custom, "w", encoding="utf-8") as handle:
                handle.write("wellID\tumi_barcode\tinternal_barcode\n")
                handle.write("a1\tacgt\ttgca\n")

            out_dir = os.path.join(tmpdir, "XPRESS_PROCESSING")
            ensure_layout(out_dir)
            config = {
                "out_dir": out_dir,
                "toolkit_directory": tmpdir,
                "sample": {"sample_type": "custom", "sample_id": None},
                "barcodes": {"barcode_file": custom},
            }

            create_barcode_tables(config)

            with open(os.path.join(config_dir(out_dir), "expect_id_barcode.tsv"), encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "wellID\tumi_barcodes\tinternal_barcodes\nA1\tACGT\tTGCA\n")
            with open(config["barcodes"]["barcode_file"], encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "ACGT\nTGCA\n")

    @unittest.skipIf(read_whitelist is None, "pandas is not installed")
    def test_custom_whitelist_matching_is_case_insensitive(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            whitelist = os.path.join(tmpdir, "whitelist.tsv")
            with open(whitelist, "w", encoding="utf-8") as handle:
                handle.write("acgt\n")
            self.assertEqual(read_whitelist(whitelist), {"ACGT"})

    def test_custom_table_rejects_barcode_assigned_to_multiple_wells(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            custom = os.path.join(tmpdir, "custom.tsv")
            with open(custom, "w", encoding="utf-8") as handle:
                handle.write("A1\tACGT\t\n")
                handle.write("A2\tACGT\t\n")

            out_dir = os.path.join(tmpdir, "XPRESS_PROCESSING")
            ensure_layout(out_dir)
            config = {
                "out_dir": out_dir,
                "toolkit_directory": tmpdir,
                "sample": {"sample_type": "custom", "sample_id": None},
                "barcodes": {"barcode_file": custom},
            }

            with self.assertRaisesRegex(ValueError, "assigned to both"):
                create_barcode_tables(config)


if __name__ == "__main__":
    unittest.main()
