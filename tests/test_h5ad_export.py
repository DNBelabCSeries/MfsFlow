import os
import sys
import unittest
from unittest import mock

from mfsflow.scripts.h5ad_export import export_h5ad


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


if __name__ == "__main__":
    unittest.main()
