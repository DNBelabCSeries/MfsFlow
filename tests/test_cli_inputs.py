import os
import sys
import gzip
import json
import tempfile
import unittest
from types import ModuleType, SimpleNamespace
from unittest import mock

from mfsflow.cli import _resume_config_path, _resume_existing_run, build_parser, main
from mfsflow.stages import MAPPING

try:
    import yaml  # noqa: F401
    from mfsflow.config.serialization import load_run_config
except ImportError:
    yaml = None
    load_run_config = None

# Some dependency-light test modules install a temporary yaml mock so they can
# import pure helpers.  That mock is not capable of serializing a resume config.
if yaml is not None and not isinstance(yaml, ModuleType):
    yaml = None
    load_run_config = None

from mfsflow.config.builder import (
    build_base_config,
    configure_reference,
    discover_fastq_pairs,
    load_samplesheet,
    resolve_samplesheet_barcodes,
    resolve_samplesheet_fastq_groups,
)


class CliInputTests(unittest.TestCase):
    def test_manual_and_custom_runs_disable_default_well_qc(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            custom_barcode = os.path.join(tmpdir, "custom.tsv")
            with open(custom_barcode, "w", encoding="utf-8") as handle:
                handle.write("wellID\tumi_barcodes\tinternal_barcodes\n")

            common = {
                "sample": "sample",
                "threads": 1,
                "stage": "Filtering",
                "tmpRoot": None,
                "outdir": tmpdir,
                "samplesheet": None,
                "fastqs": tmpdir,
                "genomeDir": tmpdir,
                "plate": None,
            }
            cases = [
                {"manual": "20", "expectBarcode": None},
                {"manual": None, "expectBarcode": custom_barcode},
            ]
            with mock.patch("mfsflow.config.builder.load_samplesheet", return_value=[]), \
                    mock.patch("mfsflow.config.builder.discover_fastq_pairs", return_value=[("R1", "R2")]), \
                    mock.patch("mfsflow.config.builder.configure_reads"), \
                    mock.patch("mfsflow.config.builder.configure_reference"):
                for case in cases:
                    args = SimpleNamespace(**common, **case)
                    config, _records = build_base_config(args, "/toolkit")
                    self.assertEqual(
                        config["well_qc"],
                        {
                            "min_reads": None,
                            "min_mapping_ratio": None,
                            "min_genes": None,
                            "min_umis": None,
                        },
                    )

    def test_later_stage_requires_explicit_resume(self):
        with self.assertRaises(SystemExit) as error, mock.patch("sys.stderr"):
            main([
                "--fastqs", "/reads",
                "--genomeDir", "/reference",
                "--sample", "sample",
                "--stage", MAPPING,
            ])

        self.assertEqual(error.exception.code, 2)

    def test_resume_config_path_accepts_root_or_processing_directory(self):
        root = os.path.abspath("/tmp/project")
        expected = os.path.join(root, "XPRESS_PROCESSING", "config", "run_config.yaml")
        self.assertEqual(_resume_config_path(root), expected)
        self.assertEqual(_resume_config_path(os.path.join(root, "XPRESS_PROCESSING")), expected)

    @unittest.skipIf(yaml is None, "PyYAML is not installed")
    def test_resume_preserves_discovered_barcode_table_and_saved_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            processing = os.path.join(tmpdir, "XPRESS_PROCESSING")
            config_dir = os.path.join(processing, "config")
            os.makedirs(config_dir)
            config_path = os.path.join(config_dir, "run_config.yaml")
            config = {
                "project": "sample",
                "out_dir": processing,
                "num_threads": 12,
                "which_Stage": "Filtering",
                "sample": {"sample_type": "discover", "discovered_sample_type": "manual"},
                "sequence_files": {"file1": {"name": "/saved/R1.fq.gz"}, "file2": {"name": "/saved/R2.fq.gz"}},
            }
            with open(config_path, "w", encoding="utf-8") as handle:
                json.dump(config, handle)
            barcode_path = os.path.join(config_dir, "expect_id_barcode.tsv")
            barcode_content = "wellID\tumi_barcodes\tinternal_barcodes\nMANUAL1\tAAAA\tCCCC\n"
            with open(barcode_path, "w", encoding="utf-8") as handle:
                handle.write(barcode_content)

            args = SimpleNamespace(
                outdir=tmpdir,
                stage=MAPPING,
                sample=None,
                threads=None,
                tmpRoot=None,
            )
            with mock.patch("mfsflow.cli._execute_pipeline") as execute:
                _resume_existing_run(args, build_parser())

            resumed = load_run_config(config_path)
            self.assertEqual(resumed["which_Stage"], MAPPING)
            self.assertEqual(resumed["num_threads"], 12)
            self.assertEqual(resumed["sequence_files"], config["sequence_files"])
            with open(barcode_path, encoding="utf-8") as handle:
                self.assertEqual(handle.read(), barcode_content)
            execute.assert_called_once()

    def test_discover_fastq_pairs_from_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            r1 = os.path.join(tmpdir, "sample_R1.fastq.gz")
            r2 = os.path.join(tmpdir, "sample_R2.fastq.gz")
            open(r1, "w").close()
            open(r2, "w").close()

            self.assertEqual(discover_fastq_pairs(tmpdir), [(r1, r2)])

    def test_samplesheet_paths_and_barcode_resolution(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            r1 = os.path.join(tmpdir, "A_R1.fq.gz")
            r2 = os.path.join(tmpdir, "A_R2.fq.gz")
            open(r1, "w").close()
            open(r2, "w").close()

            sheet = os.path.join(tmpdir, "samplesheet.csv")
            with open(sheet, "w") as f:
                f.write("read1,read2,barcode\n")
                f.write("A_R1.fq.gz,A_R2.fq.gz,CCCC\n")

            expect = os.path.join(tmpdir, "expect_id_barcode.tsv")
            with open(expect, "w") as f:
                f.write("wellID\tumi_barcodes\tinternal_barcodes\n")
                f.write("P1A1\tAAAA,CCCC\tTTTT\n")

            records = load_samplesheet(sheet, tmpdir)
            resolved = resolve_samplesheet_barcodes(records, expect)

            self.assertEqual(resolved[0]["read1"], r1)
            self.assertEqual(resolved[0]["read2"], r2)
            self.assertEqual(resolved[0]["barcode"], "CCCC")
            self.assertEqual(resolved[0]["wellID"], "P1A1")
            self.assertEqual(resolved[0]["barcode_type"], "umi")

    def test_samplesheet_discover_defers_duplicate_barcode_resolution(self):
        records = [{
            "read1": "/tmp/A_R1.fq.gz",
            "read2": "/tmp/A_R2.fq.gz",
            "barcode": "CACTCACACGAAGCAGCCGA",
        }]
        with tempfile.TemporaryDirectory() as tmpdir:
            expect = os.path.join(tmpdir, "expect_id_barcode.tsv")
            with open(expect, "w") as f:
                f.write("wellID\tumi_barcodes\tinternal_barcodes\n")
                f.write("MANUAL1\tCACTCACACGAAGCAGCCGA\tAAAA\n")
                f.write("P1A1\tCACTCACACGAAGCAGCCGA\tTTTT\n")

            resolved = resolve_samplesheet_fastq_groups(records, expect, "discover")
            self.assertEqual(resolved, records)

            with self.assertRaisesRegex(ValueError, "duplicated"):
                resolve_samplesheet_fastq_groups(records, expect, "manual")

    def test_discover_fastq_pairs_rejects_missing_r2(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            r1 = os.path.join(tmpdir, "sample_R1.fastq.gz")
            open(r1, "w").close()

            with self.assertRaisesRegex(FileNotFoundError, "Missing matching R2 FASTQ"):
                discover_fastq_pairs(tmpdir)

    def test_samplesheet_rejects_duplicate_fastq_pairs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            r1 = os.path.join(tmpdir, "A_R1.fq.gz")
            r2 = os.path.join(tmpdir, "A_R2.fq.gz")
            open(r1, "w").close()
            open(r2, "w").close()

            sheet = os.path.join(tmpdir, "samplesheet.csv")
            with open(sheet, "w") as f:
                f.write("read1,read2,barcode\n")
                f.write("A_R1.fq.gz,A_R2.fq.gz,CCCC\n")
                f.write("A_R1.fq.gz,A_R2.fq.gz,TTTT\n")

            with self.assertRaisesRegex(ValueError, "Duplicate read1/read2 pair"):
                load_samplesheet(sheet, tmpdir)

    def test_configure_reference_accepts_gzipped_gtf(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            os.makedirs(os.path.join(tmpdir, "genes"))
            os.makedirs(os.path.join(tmpdir, "star"))
            gz_gtf = os.path.join(tmpdir, "genes", "genes.gtf.gz")
            with gzip.open(gz_gtf, "wt") as handle:
                handle.write('chr1\tT\texon\t1\t10\t.\t+\t.\tgene_id "g1";\n')

            config = {"reference": {}}
            configure_reference(config, tmpdir)

            self.assertEqual(config["reference"]["GTF_file"], gz_gtf)
            self.assertEqual(config["reference"]["STAR_index"], os.path.join(tmpdir, "star"))


if __name__ == "__main__":
    unittest.main()
