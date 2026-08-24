import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from mfsflow.report import (
    calculate_summary_metrics,
    generate_multi_report,
    _process_sequencing_quality_data,
    _infer_transcriptome_label,
    export_deliverables_to_outs,
)
from mfsflow.scripts.generate_stats import calculate_read_ratios


class ReportMetadataTests(unittest.TestCase):
    def test_genic_ratio_matches_report_definition(self):
        mapping, genic, legacy_exon_intron = calculate_read_ratios(
            exon_reads=60,
            intron_reads=20,
            intergenic_reads=10,
            ambiguity_reads=0,
            unmapped_reads=10,
            other_unassigned_reads=0,
            all_reads=100,
        )
        self.assertAlmostEqual(mapping, 0.9)
        self.assertAlmostEqual(genic, 0.8)
        self.assertAlmostEqual(legacy_exon_intron, 3.0)

    def test_mapping_ratio_includes_mapped_unassigned_fragments(self):
        mapping, genic, _ = calculate_read_ratios(
            exon_reads=50,
            intron_reads=20,
            intergenic_reads=10,
            ambiguity_reads=0,
            unmapped_reads=10,
            other_unassigned_reads=10,
            all_reads=100,
        )
        self.assertAlmostEqual(mapping, 0.9)
        self.assertAlmostEqual(genic, 0.7)

    def test_empty_partial_report_is_written_atomically(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            outdir.mkdir()
            report_path = generate_multi_report(
                "sample",
                str(outdir),
                {
                    "project": "sample",
                    "out_dir": str(outdir),
                    "sample": {"sample_type": "manual"},
                    "reference": {},
                },
            )

            self.assertTrue(report_path.exists())
            self.assertFalse(report_path.with_suffix(report_path.suffix + ".tmp").exists())

    def test_failed_analysis_report_has_incomplete_warning(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            outdir.mkdir()
            report_path = generate_multi_report(
                "sample",
                str(outdir),
                {
                    "project": "sample",
                    "out_dir": str(outdir),
                    "sample": {"sample_type": "manual"},
                    "reference": {},
                    "_analysis_failed": True,
                },
            )

            self.assertIn("Incomplete analysis:", report_path.read_text())

    def test_transcriptome_label_uses_parent_for_star_index_dir(self):
        label = _infer_transcriptome_label({"STAR_index": "/path/to/reference/star"})
        self.assertEqual(label, "reference")

    def test_transcriptome_label_uses_index_basename_when_specific(self):
        label = _infer_transcriptome_label({"STAR_index": "/path/to/GRCh38_2024"})
        self.assertEqual(label, "GRCh38_2024")

    def test_transcriptome_label_prefers_explicit_config(self):
        label = _infer_transcriptome_label({
            "STAR_index": "/path/to/reference/star",
            "transcriptome_name": "custom_v1",
        })
        self.assertEqual(label, "custom_v1")

    def test_export_deliverables_to_outs_preserves_working_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            outdir = root / "XPRESS_PROCESSING"
            outs = root / "outs"
            expr = outdir / "expression"
            stats = outdir / "stats"
            config = outdir / "config"
            mex = expr / "Sample01.inex.umi"
            for path in (mex, stats, config):
                path.mkdir(parents=True)

            (expr / "Sample01.h5ad").write_text("h5ad")
            (mex / "matrix.mtx.gz").write_text("matrix")
            (mex / "features.tsv.gz").write_text("features")
            (mex / "barcodes.tsv.gz").write_text("barcodes")
            (stats / "Sample01.stats.tsv").write_text("stats")
            (stats / "Sample01.read_stats.json").write_text("{}")
            (outdir / "Sample01.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam").write_text("bam")
            (outdir / "Sample01.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam.bai").write_text("bai")
            (config / "run_config.yaml").write_text("project: Sample01")
            (config / "expect_id_barcode.tsv").write_text("wellID\tumi_barcodes\tinternal_barcodes\n")
            (outs / "stats").mkdir(parents=True)
            (outs / "stats" / "Sample01.obsolete.tsv").write_text("stale")
            (outs / "stats" / "Other.keep.tsv").write_text("other project")

            export_deliverables_to_outs(outdir, outs, "Sample01")

            self.assertTrue((outs / "expression" / "Sample01.h5ad").exists())
            self.assertTrue((outs / "expression" / "Sample01.inex.umi" / "matrix.mtx.gz").exists())
            self.assertTrue((outs / "stats" / "Sample01.stats.tsv").exists())
            self.assertTrue((outs / "stats" / "Sample01.read_stats.json").exists())
            self.assertTrue((outs / "bam" / "Sample01.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam").exists())
            self.assertTrue((outs / "bam" / "Sample01.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam.bai").exists())
            self.assertTrue((outs / "config" / "run_config.yaml").exists())
            self.assertTrue((outs / "config" / "expect_id_barcode.tsv").exists())
            self.assertTrue((expr / "Sample01.h5ad").exists())
            self.assertTrue(mex.exists())
            self.assertTrue((stats / "Sample01.stats.tsv").exists())
            self.assertTrue((outdir / "Sample01.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam").exists())
            self.assertTrue((config / "run_config.yaml").exists())
            self.assertFalse((outs / "stats" / "Sample01.obsolete.tsv").exists())
            self.assertTrue((outs / "stats" / "Other.keep.tsv").exists())

    def test_failed_report_does_not_export_or_replace_success_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            outdir = root / "XPRESS_PROCESSING"
            expression = outdir / "expression"
            outs = root / "outs"
            expression.mkdir(parents=True)
            outs.mkdir()
            working_h5ad = expression / "sample.h5ad"
            working_h5ad.write_text("working")
            success_report = outs / "sample_manual_report.html"
            success_report.write_text("previous-success")

            partial = generate_multi_report(
                "sample",
                str(outdir),
                {
                    "project": "sample",
                    "out_dir": str(outdir),
                    "sample": {"sample_type": "manual"},
                    "reference": {},
                    "_analysis_failed": True,
                },
            )

            self.assertEqual(partial.name, "sample_partial_report.html")
            self.assertEqual(success_report.read_text(), "previous-success")
            self.assertTrue(working_h5ad.exists())
            self.assertFalse((outs / "expression" / "sample.h5ad").exists())
            self.assertTrue((outs / "partial_run_manifest.json").exists())

    def test_success_report_is_not_replaced_when_deliverable_export_fails(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            processing = Path(tmpdir) / "XPRESS_PROCESSING"
            processing.mkdir()
            outs = Path(tmpdir) / "outs"
            outs.mkdir()
            report_path = outs / "sample_manual_report.html"
            report_path.write_text("previous successful report", encoding="utf-8")
            manifest_path = outs / "run_manifest.json"
            manifest_path.write_text('{"status": "success"}\n', encoding="utf-8")
            config = {
                "project": "sample",
                "out_dir": str(processing),
                "sample": {"sample_type": "manual"},
                "sequence_files": {"file1": {"name": "reads.fq.gz"}},
                "reference": {},
            }

            with mock.patch(
                "mfsflow.report.export_deliverables_to_outs",
                side_effect=OSError("disk full"),
            ):
                with self.assertRaisesRegex(OSError, "disk full"):
                    generate_multi_report("sample", str(processing), config)

            self.assertEqual(report_path.read_text(encoding="utf-8"), "previous successful report")
            self.assertFalse(manifest_path.exists())

    def test_sequencing_quality_summary_uses_q30_and_bcstats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            stats = outdir / "stats"
            stats.mkdir(parents=True)
            (outdir / "sample.BCstats.txt").write_text("AAAA\t60\nCCCC\t40\n")
            (stats / "sample.q30_stats.tsv").write_text(
                "metric\ttotal_bases\tq30_bases\tq30_rate\n"
                "R1\t9000\t8100\t0.900000\n"
                "BC\t2000\t1800\t0.900000\n"
                "UMI\t1000\t950\t0.950000\n"
                "R1_cDNA\t8000\t7200\t0.900000\n"
                "R2_cDNA\t9000\t7650\t0.850000\n"
            )
            context = {
                "_run_config": {
                    "sequence_files": {
                        "file1": {"base_definition": ["cDNA(11-90)", "UMI(1-10)"]},
                        "file2": {"base_definition": ["cDNA(1-90)", "BC(91-110)"]},
                    }
                }
            }

            _process_sequencing_quality_data(outdir, context)

            self.assertIn('"Total sequencing reads", "value": "100"', context["sequencing_quality_summary_data"])
            self.assertIn('"Valid barcode reads", "value": "100"', context["sequencing_quality_summary_data"])
            self.assertIn('"Unused barcode reads", "value": "0"', context["sequencing_quality_summary_data"])
            self.assertIn('"Valid barcode rate", "value": "100.0%"', context["sequencing_quality_summary_data"])
            self.assertIn('"Read2 cDNA Q30", "value": "85.0%"', context["sequencing_quality_summary_data"])

    def test_sequencing_quality_summary_prefers_read_stats_counts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            stats = outdir / "stats"
            stats.mkdir(parents=True)
            (outdir / "sample.BCstats.txt").write_text("AAAA\t999\n")
            (stats / "sample.q30_stats.tsv").write_text(
                "metric\ttotal_bases\tq30_bases\tq30_rate\n"
                "R1\t9000\t8100\t0.900000\n"
            )
            (stats / "sample.read_stats.json").write_text(
                '{"read_stats": {"BC1": {"UMI_Reads": 10, "Internal_Reads": 5}, "__NO_CB__": {"Unused BC": 5}}}'
            )
            context = {
                "_run_config": {
                    "sequence_files": {
                        "file1": {"base_definition": ["cDNA(11-90)", "UMI(1-10)"]},
                    }
                }
            }

            _process_sequencing_quality_data(outdir, context)

            self.assertIn('"Total sequencing reads", "value": "20"', context["sequencing_quality_summary_data"])
            self.assertIn('"Valid barcode reads", "value": "15"', context["sequencing_quality_summary_data"])
            self.assertIn('"Unused barcode reads", "value": "5"', context["sequencing_quality_summary_data"])
            self.assertIn('"Valid barcode rate", "value": "75.0%"', context["sequencing_quality_summary_data"])

    def test_report_loaders_use_current_project_when_stale_files_exist(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            stats = outdir / "stats"
            stats.mkdir(parents=True)
            (outdir / "sample.BCstats.txt").write_text("AAAA\t100\n")
            (outdir / "stale.BCstats.txt").write_text("AAAA\t9999\n")
            for project, total in (("sample", 9000), ("stale", 9999000)):
                (stats / f"{project}.q30_stats.tsv").write_text(
                    "metric\ttotal_bases\tq30_bases\tq30_rate\n"
                    f"R1\t{total}\t{total}\t1.0\n"
                )
            context = {
                "_run_config": {
                    "project": "sample",
                    "sequence_files": {"file1": {"base_definition": ["cDNA(1-90)"]}},
                }
            }

            _process_sequencing_quality_data(outdir, context)

            self.assertIn('"Total sequencing reads", "value": "100"', context["sequencing_quality_summary_data"])
            self.assertNotIn("11,110", context["sequencing_quality_summary_data"])

    def test_report_does_not_use_another_projects_only_artifact(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            stats = outdir / "stats"
            stats.mkdir(parents=True)
            (stats / "old_sample.stats.tsv").write_text(
                "wellID\tIntron_Exon_umis\nP1A1\t10\n",
                encoding="utf-8",
            )

            metrics = calculate_summary_metrics(outdir, project="current_sample")

            self.assertEqual(metrics["rna_estm_Num_cell"], "")


if __name__ == "__main__":
    unittest.main()
