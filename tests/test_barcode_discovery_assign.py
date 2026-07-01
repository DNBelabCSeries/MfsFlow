"""Tests for run_barcode_discovery assign_barcodes behavior."""

import unittest
from unittest import mock

from mfsflow.bootstrap import run_barcode_discovery


_FAKE_SELECTED = [{
    "candidate_type": "manual",
    "candidate_id": "20",
    "matched_reads": 1000,
    "matched_expected_barcodes": 5,
}], []


class BarcodeDiscoveryAssignTests(unittest.TestCase):
    def test_assign_barcodes_false_skips_write_expected_tables(self):
        config = {"sample": {}, "barcodes": {}, "toolkit_directory": "."}
        with mock.patch("mfsflow.bootstrap.discover_barcodes", return_value=_FAKE_SELECTED), \
             mock.patch("mfsflow.bootstrap.write_expected_tables") as write_tables, \
             mock.patch("mfsflow.bootstrap.build_expected_records", return_value=[]):
            run_barcode_discovery(config, "proj", "/tmp/analysis", assign_barcodes=False)

        write_tables.assert_not_called()
        self.assertEqual(config["sample"]["discovered_sample_type"], "manual")
        self.assertEqual(config["sample"]["discovered_sample_ids"], "20")
        self.assertNotIn("barcode_file", config["barcodes"])

    def test_assign_barcodes_true_writes_tables_and_sets_barcode_file(self):
        config = {"sample": {}, "barcodes": {}, "toolkit_directory": "."}
        with mock.patch("mfsflow.bootstrap.discover_barcodes", return_value=_FAKE_SELECTED), \
             mock.patch("mfsflow.bootstrap.write_expected_tables",
                        return_value=("/tmp/pipe.tsv", "/tmp/summary.tsv")) as write_tables, \
             mock.patch("mfsflow.bootstrap.build_expected_records", return_value=[]):
            run_barcode_discovery(config, "proj", "/tmp/analysis")

        write_tables.assert_called_once()
        self.assertEqual(config["barcodes"]["barcode_file"], "/tmp/pipe.tsv")
        self.assertEqual(config["sample"]["discovered_sample_type"], "manual")


if __name__ == "__main__":
    unittest.main()
