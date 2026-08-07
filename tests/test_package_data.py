import os
import pkgutil
import unittest
from pathlib import Path

import mfsflow


REQUIRED_RESOURCES = (
    "yaml/example_config.yaml",
    "yaml/manual_barcode_list.yaml",
    "yaml/auto_barcode_list.yaml",
    "report_assets/template_auto.html",
    "report_assets/template_manual.html",
    "report_assets/plotly-2.26.0.min.js",
    "report_assets/logo.png",
    "software/STAR",
    "software/featureCounts",
    "software/pigz",
    "software/samtools",
    "software/seqkit",
)


class PackageDataTests(unittest.TestCase):
    def test_required_resources_are_available_from_installed_package(self):
        for resource in REQUIRED_RESOURCES:
            with self.subTest(resource=resource):
                self.assertIsNotNone(
                    pkgutil.get_data("mfsflow", resource),
                    f"Missing packaged resource: {resource}",
                )

    def test_packaged_bioinformatics_tools_are_executable(self):
        package_dir = Path(mfsflow.__file__).resolve().parent
        for resource in REQUIRED_RESOURCES:
            if not resource.startswith("software/"):
                continue
            with self.subTest(resource=resource):
                path = package_dir / resource
                self.assertTrue(path.is_file(), f"Missing installed tool: {resource}")
                self.assertTrue(os.access(path, os.X_OK), f"Installed tool is not executable: {resource}")


if __name__ == "__main__":
    unittest.main()
