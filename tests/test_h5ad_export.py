import json
import os
import sys
import unittest
from unittest import mock

from mfsflow.scripts.h5ad_export import _json_safe, export_h5ad


class H5adExportTests(unittest.TestCase):
    def test_missing_anndata_reports_actionable_error(self):
        real_import = __import__

        def fake_import(name, *args, **kwargs):
            if name == "anndata":
                raise ImportError("missing anndata")
            return real_import(name, *args, **kwargs)

        with mock.patch("builtins.__import__", side_effect=fake_import):
            with self.assertRaisesRegex(RuntimeError, "H5AD export requires anndata"):
                export_h5ad("/tmp/does-not-matter", "project")

    def test_json_safe_converts_list_of_dict_to_json_string(self):
        config = {
            "project": "F350032621",
            "fastq_groups": [
                {"read1": "/data/r1.fq.gz", "read2": "/data/r2.fq.gz", "barcode": "ACGTACGT"},
                {"read1": "/data/r1_2.fq.gz", "read2": "/data/r2_2.fq.gz", "barcode": "TTTTAAAA"},
            ],
            "sample": {"sample_type": "discover"},
        }
        result = _json_safe(config)
        self.assertEqual(result["project"], "F350032621")
        self.assertEqual(result["sample"], {"sample_type": "discover"})
        self.assertIsInstance(result["fastq_groups"], str)
        parsed = json.loads(result["fastq_groups"])
        self.assertEqual(len(parsed), 2)
        self.assertEqual(parsed[0]["barcode"], "ACGTACGT")

    def test_json_safe_preserves_plain_list_and_scalar(self):
        result = _json_safe({"threads": 20, "stages": ["Filtering", "Mapping"]})
        self.assertEqual(result["threads"], 20)
        self.assertEqual(result["stages"], ["Filtering", "Mapping"])


if __name__ == "__main__":
    unittest.main()
