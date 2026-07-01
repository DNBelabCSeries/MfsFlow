import os
import sys
import tempfile
import unittest
from unittest import mock

sys.modules.setdefault("yaml", mock.Mock())
from mfsflow.barcode_discovery import discover_barcodes, write_expected_tables


class BarcodeDiscoveryTests(unittest.TestCase):
    def test_discover_selects_best_candidate_and_writes_selected_tables(self):
        records = [
            {
                "candidate_type": "manual",
                "candidate_id": "9",
                "wellID": "MANUAL9",
                "barcode": "AAAAAAAAAAAAAAAAAAAA",
                "barcode_type": "umi",
            },
            {
                "candidate_type": "manual",
                "candidate_id": "9",
                "wellID": "MANUAL9",
                "barcode": "CCCCCCCCCCCCCCCCCCCC",
                "barcode_type": "internal",
            },
            {
                "candidate_type": "manual",
                "candidate_id": "10",
                "wellID": "MANUAL10",
                "barcode": "GGGGGGGGGGGGGGGGGGGG",
                "barcode_type": "umi",
            },
            {
                "candidate_type": "manual",
                "candidate_id": "10",
                "wellID": "MANUAL10",
                "barcode": "TTTTTTTTTTTTTTTTTTTT",
                "barcode_type": "internal",
            },
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            bcstats = os.path.join(tmpdir, "sample.BCstats.txt")
            with open(bcstats, "w") as handle:
                handle.write("AAAAAAAAAAAAAAAAAAAA\t100\n")
                handle.write("CCCCCCCCCCCCCCCCCCCA\t80\n")
                handle.write("GGGGGGGGGGGGGGGGGGGG\t5\n")

            report = os.path.join(tmpdir, "discovery.tsv")
            selected, selected_records = discover_barcodes(bcstats, records, report)
            self.assertEqual([("manual", "9")], [(x["candidate_type"], x["candidate_id"]) for x in selected])
            self.assertEqual(2, len(selected_records))
            self.assertTrue(os.path.exists(report))

            pipe, summary = write_expected_tables(selected_records, tmpdir)
            with open(pipe) as handle:
                barcodes = {line.strip() for line in handle if line.strip()}
            self.assertEqual({"AAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCC"}, barcodes)
            with open(summary) as handle:
                self.assertIn("MANUAL9", handle.read())

    def test_discover_counts_shared_barcodes_once_per_candidate(self):
        records = [
            {
                "candidate_type": "manual",
                "candidate_id": "9",
                "wellID": "MANUAL9",
                "barcode": "ACGTACGTACGTACGTAAAA",
                "barcode_type": "umi",
            },
            {
                "candidate_type": "manual",
                "candidate_id": "9",
                "wellID": "MANUAL9",
                "barcode": "ACGTACGTACGTACGTCCCC",
                "barcode_type": "internal",
            },
            {
                "candidate_type": "auto",
                "candidate_id": "1",
                "wellID": "P1A1",
                "barcode": "ACGTACGTACGTACGTAAAA",
                "barcode_type": "umi",
            },
            {
                "candidate_type": "auto",
                "candidate_id": "1",
                "wellID": "P1A1",
                "barcode": "ACGTACGTACGTACGTCCCC",
                "barcode_type": "internal",
            },
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            bcstats = os.path.join(tmpdir, "sample.BCstats.txt")
            with open(bcstats, "w") as handle:
                handle.write("ACGTACGTACGTACGTAAAA\t100\n")
                handle.write("ACGTACGTACGTACGTCCCC\t80\n")

            report = os.path.join(tmpdir, "discovery.tsv")
            with self.assertRaisesRegex(ValueError, "could not distinguish manual vs auto"):
                discover_barcodes(bcstats, records, report, min_unique_barcodes=1)

            with open(report) as handle:
                report_text = handle.read()
            self.assertIn("manual\t9\t180", report_text)
            self.assertIn("auto\t1\t180", report_text)


if __name__ == "__main__":
    unittest.main()
