import os
import json
import re
import shutil
import sys
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from mfsflow.report import (
    calculate_summary_metrics,
    generate_multi_report,
    _process_sequencing_quality_data,
    _process_barcode_report_data,
    _infer_transcriptome_label,
    export_deliverables_to_outs,
)
from mfsflow.scripts.generate_stats import calculate_read_ratios


class ReportMetadataTests(unittest.TestCase):
    def test_well_qc_status_marks_all_failed_for_js_fallback(self):
        # well_qc_status must carry __all_failed__ so the JS summary tables can
        # use all-wells medians instead of rendering blanks for the Active set.
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            config_dir = outdir / "config"
            config_dir.mkdir(parents=True)
            (config_dir / "expect_id_barcode.tsv").write_text(
                "wellID\tumi_barcodes\tinternal_barcodes\n"
                "P1A1\tAAAA\tCCCC\n"
            )
            context = {
                "sample_type": "auto",
                "rna_stats_table_data": json.dumps([
                    {
                        "wellID": "P1A1",
                        "all_reads": 200,
                        "internal_reads": 100,
                        "umi_reads": 100,
                        "MappingRatio": 0.2,
                        "Intron_Exon_genes": 40,
                        "Intron_Exon_umis": 50,
                    },
                ]),
            }
            config = {"sample": {"sample_type": "auto"}}

            _process_barcode_report_data(outdir, context, config)

            status = json.loads(context["well_qc_status_data"])
            self.assertTrue(status.get("__all_failed__"))

    def test_expected_well_matching_is_case_insensitive(self):
        # expect_id_barcode.tsv and stats.tsv may differ in case ("p1a1" vs
        # "P1A1"); a strict match would mark every well unexpected and zero the
        # Active count, so matching must be case-insensitive.
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            config_dir = outdir / "config"
            config_dir.mkdir(parents=True)
            (config_dir / "expect_id_barcode.tsv").write_text(
                "wellID\tumi_barcodes\tinternal_barcodes\n"
                "p1a1\tAAAA\tCCCC\n"
            )
            context = {
                "sample_type": "auto",
                "rna_stats_table_data": json.dumps([
                    {
                        "wellID": "P1A1",
                        "all_reads": 2000,
                        "internal_reads": 1000,
                        "umi_reads": 1000,
                        "MappingRatio": 0.8,
                        "Intron_Exon_genes": 150,
                        "Intron_Exon_umis": 200,
                    },
                ]),
            }
            config = {"sample": {"sample_type": "auto"}}

            _process_barcode_report_data(outdir, context, config)

            summary = json.loads(context["barcode_report_summary_data"])
            self.assertEqual(summary["expected_wells"], 1)
            self.assertEqual(summary["active_wells"], 1)

    def test_active_wells_require_all_enabled_qc_thresholds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            config_dir = outdir / "config"
            config_dir.mkdir(parents=True)
            (config_dir / "expect_id_barcode.tsv").write_text(
                "wellID\tumi_barcodes\tinternal_barcodes\n"
                "P1A1\tAAAA\tCCCC\n"
                "P1A2\tGGGG\tTTTT\n"
            )
            context = {
                "sample_type": "auto",
                "rna_stats_table_data": json.dumps([
                    {
                        "wellID": "P1A1",
                        "all_reads": 2000,
                        "internal_reads": 1000,
                        "umi_reads": 1000,
                        "MappingRatio": 0.8,
                        "Intron_Exon_genes": 150,
                        "Intron_Exon_umis": 200,
                    },
                    {
                        "wellID": "P1A2",
                        "all_reads": 2000,
                        "internal_reads": 1000,
                        "umi_reads": 1000,
                        "MappingRatio": 0.8,
                        "Intron_Exon_genes": 80,
                        "Intron_Exon_umis": 99,
                    },
                    {
                        "wellID": "UNEXPECTED",
                        "all_reads": 5000,
                        "internal_reads": 2500,
                        "umi_reads": 2500,
                        "MappingRatio": 0.9,
                        "Intron_Exon_genes": 900,
                        "Intron_Exon_umis": 900,
                    },
                ]),
            }
            config = {
                "sample": {"sample_type": "auto"},
                "well_qc": {
                    "min_reads": 1000,
                    "min_mapping_ratio": 0.30,
                    "min_genes": 100,
                    "min_umis": 100,
                },
            }

            _process_barcode_report_data(outdir, context, config)

            summary = json.loads(context["barcode_report_summary_data"])
            self.assertEqual(summary["expected_wells"], 2)
            self.assertEqual(summary["active_wells"], 1)
            self.assertEqual(summary["median_reads"], 2000)
            self.assertEqual(summary["median_genes"], 150)
            self.assertEqual(summary["median_umis"], 200)
            self.assertAlmostEqual(summary["median_mapping_ratio"], 0.8)
            cards = json.loads(context["barcode_mode_cards_data"])
            active_card = next(card for card in cards if card["label"] == "Active wells")
            self.assertIn("UMIs: >= 100", active_card["help"])
            median_umis_card = next(card for card in cards if card["label"] == "Median UMIs")
            self.assertEqual(median_umis_card["value"], "200")
            self.assertIn("Active wells", median_umis_card["help"])

    def test_manual_mode_does_not_apply_default_well_qc_thresholds(self):
        # Manual runs list wells explicitly, so wells are not re-filtered by
        # the absolute thresholds unless the user configured well_qc.
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            config_dir = outdir / "config"
            config_dir.mkdir(parents=True)
            (config_dir / "expect_id_barcode.tsv").write_text(
                "wellID\tumi_barcodes\tinternal_barcodes\n"
                "P1A1\tAAAA\tCCCC\n"
                "P1A2\tGGGG\tTTTT\n"
            )
            context = {
                "sample_type": "manual",
                "rna_stats_table_data": json.dumps([
                    {
                        "wellID": "P1A1",
                        "all_reads": 200,
                        "internal_reads": 100,
                        "umi_reads": 100,
                        "MappingRatio": 0.5,
                        "Intron_Exon_genes": 40,
                        "Intron_Exon_umis": 50,
                    },
                    {
                        "wellID": "P1A2",
                        "all_reads": 150,
                        "internal_reads": 75,
                        "umi_reads": 75,
                        "MappingRatio": 0.5,
                        "Intron_Exon_genes": 30,
                        "Intron_Exon_umis": 40,
                    },
                ]),
            }
            config = {"sample": {"sample_type": "manual"}}

            _process_barcode_report_data(outdir, context, config)

            summary = json.loads(context["barcode_report_summary_data"])
            self.assertEqual(summary["expected_wells"], 2)
            self.assertEqual(summary["active_wells"], 2)
            self.assertFalse(summary["qc_all_failed"])
            self.assertEqual(summary["median_reads"], 175)
            self.assertEqual(summary["median_genes"], 35)

    def test_all_wells_failing_qc_falls_back_to_all_well_medians(self):
        # If every well fails QC, median cards must not go blank: fall back to
        # all wells and flag it so the report can show an explicit warning.
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            config_dir = outdir / "config"
            config_dir.mkdir(parents=True)
            (config_dir / "expect_id_barcode.tsv").write_text(
                "wellID\tumi_barcodes\tinternal_barcodes\n"
                "P1A1\tAAAA\tCCCC\n"
                "P1A2\tGGGG\tTTTT\n"
            )
            context = {
                "sample_type": "auto",
                "rna_stats_table_data": json.dumps([
                    {
                        "wellID": "P1A1",
                        "all_reads": 200,
                        "internal_reads": 100,
                        "umi_reads": 100,
                        "MappingRatio": 0.2,
                        "Intron_Exon_genes": 40,
                        "Intron_Exon_umis": 50,
                    },
                    {
                        "wellID": "P1A2",
                        "all_reads": 150,
                        "internal_reads": 75,
                        "umi_reads": 75,
                        "MappingRatio": 0.2,
                        "Intron_Exon_genes": 30,
                        "Intron_Exon_umis": 40,
                    },
                    {
                        # This diagnostic row must not enter the all-wells
                        # fallback because it is absent from the expected table.
                        "wellID": "UNEXPECTED",
                        "all_reads": 10000,
                        "internal_reads": 5000,
                        "umi_reads": 5000,
                        "MappingRatio": 0.2,
                        "Intron_Exon_genes": 3000,
                        "Intron_Exon_umis": 4000,
                    },
                ]),
            }
            config = {"sample": {"sample_type": "auto"}}

            _process_barcode_report_data(outdir, context, config)

            summary = json.loads(context["barcode_report_summary_data"])
            self.assertEqual(summary["active_wells"], 0)
            self.assertTrue(summary["qc_all_failed"])
            self.assertEqual(summary["median_reads"], 175)
            self.assertEqual(summary["median_genes"], 35)
            cards = json.loads(context["barcode_mode_cards_data"])
            active_card = next(card for card in cards if card["label"] == "Active wells")
            self.assertIsNotNone(active_card.get("warning"))
            median_card = next(card for card in cards if card["label"] == "Median reads")
            self.assertIn("all wells (no well passed QC)", median_card["help"])

    def test_report_contains_qc_help_and_per_well_status(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            stats_dir = outdir / "stats"
            config_dir = outdir / "config"
            stats_dir.mkdir(parents=True)
            config_dir.mkdir()
            (config_dir / "expect_id_barcode.tsv").write_text(
                "wellID\tumi_barcodes\tinternal_barcodes\nP1A1\tAAAA\tCCCC\n"
            )
            (stats_dir / "sample.stats.tsv").write_text(
                "wellID\tinternal_reads\tumi_reads\tall_reads\tMappingRatio\t"
                "Intron_Exon_genes\tIntron_Exon_umis\n"
                "P1A1\t1000\t1000\t2000\t0.8\t150\t200\n"
            )

            report_path = generate_multi_report(
                "sample",
                str(outdir),
                {
                    "project": "sample",
                    "out_dir": str(outdir),
                    "sequence_files": ["sample_R1.fastq.gz"],
                    "sample": {"sample_type": "manual"},
                    "reference": {},
                    "well_qc": {
                        "min_reads": 1000,
                        "min_mapping_ratio": 0.30,
                        "min_genes": 100,
                        "min_umis": 100,
                    },
                },
            )

            html = report_path.read_text(encoding="utf-8")
            self.assertIn("helpButton.className = 'help-button'", html)
            self.assertIn("A well is counted as Active", html)
            self.assertIn("well_status", html)

    def test_generated_report_inline_javascript_passes_node_syntax_check(self):
        node = shutil.which("node")
        if not node:
            self.skipTest("Node.js is not installed")
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            stats_dir = outdir / "stats"
            config_dir = outdir / "config"
            stats_dir.mkdir(parents=True)
            config_dir.mkdir()
            (config_dir / "expect_id_barcode.tsv").write_text(
                "wellID\tumi_barcodes\tinternal_barcodes\nP1A1\tAAAA\tCCCC\n"
            )
            (stats_dir / "sample.stats.tsv").write_text(
                "wellID\tinternal_reads\tumi_reads\tall_reads\tMappingRatio\t"
                "Intron_Exon_genes\tIntron_Exon_umis\n"
                "P1A1\t1000\t1000\t2000\t0.8\t150\t200\n"
            )

            report_path = generate_multi_report(
                "sample",
                str(outdir),
                {
                    "project": "sample",
                    "out_dir": str(outdir),
                    "sample": {"sample_type": "auto"},
                    "reference": {},
                },
            )
            html = report_path.read_text(encoding="utf-8")
            scripts = [script for script in re.findall(r"<script[^>]*>(.*?)</script>", html, re.S) if script.strip()]
            self.assertGreaterEqual(len(scripts), 1)
            for index, script in enumerate(scripts):
                script_path = Path(tmpdir) / f"report_script_{index}.js"
                script_path.write_text(script, encoding="utf-8")
                result = subprocess.run(
                    [node, "--check", str(script_path)],
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, 0, result.stderr)

    def test_report_populates_summary_cards_for_legacy_config(self):
        # Older saved configs may not contain sequence_files even though the
        # report artifacts are present. The summary area must remain visible.
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "XPRESS_PROCESSING"
            stats_dir = outdir / "stats"
            config_dir = outdir / "config"
            stats_dir.mkdir(parents=True)
            config_dir.mkdir()
            (config_dir / "expect_id_barcode.tsv").write_text(
                "wellID\tumi_barcodes\tinternal_barcodes\nP1A1\tAAAA\tCCCC\n"
            )
            (stats_dir / "sample.stats.tsv").write_text(
                "wellID\tinternal_reads\tumi_reads\tall_reads\tMappingRatio\t"
                "Intron_Exon_genes\tIntron_Exon_umis\n"
                "P1A1\t1000\t1000\t2000\t0.8\t150\t200\n"
            )

            report_path = generate_multi_report(
                "sample",
                str(outdir),
                {
                    "project": "sample",
                    "out_dir": str(outdir),
                    "sample": {"sample_type": "auto"},
                    "reference": {},
                },
            )

            html = report_path.read_text(encoding="utf-8")
            cards_line = next(line for line in html.splitlines() if "const cards = " in line)
            cards_payload = cards_line.split("const cards = ", 1)[1].rsplit(";", 1)[0]
            cards = json.loads(cards_payload)
            self.assertEqual(len(cards), 7)
            self.assertIn("Expected wells", [card["label"] for card in cards])
            self.assertIn("Active wells", [card["label"] for card in cards])

    def test_auto_and_manual_templates_expose_active_status_interactions(self):
        template_dir = Path(__file__).resolve().parents[1] / "mfsflow" / "report_assets"
        for template_name in ("template_auto.html", "template_manual.html"):
            html = (template_dir / template_name).read_text(encoding="utf-8")
            self.assertIn('id="summary-cards"', html)
            self.assertIn("Status=%{customdata[0]}", html)
            self.assertIn("sortCol === 'well_status'", html)
            self.assertIn("wellStatus(row)", html)
            self.assertIn("const cards = ${barcode_mode_cards_data};", html)
            self.assertIn("const stats = ${rna_stats_table_data};", html)
            self.assertIn("const wellQcRaw = ${well_qc_status_data};", html)

        auto_html = (template_dir / "template_auto.html").read_text(encoding="utf-8")
        manual_html = (template_dir / "template_manual.html").read_text(encoding="utf-8")
        self.assertIn("Status=%{customdata[1]}", auto_html)
        self.assertIn("return qc ? qc.active === true : !(row && row.active === false);", auto_html)
        self.assertIn("return qc ? qc.active === true : !(row && row.active === false);", manual_html)

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
            # Without read_stats.json, Unused is NA: raw-input minus BCstats would
            # also count reads fqfilter dropped for quality (inflated number).
            self.assertIn('"Unused barcode reads", "value": "NA"', context["sequencing_quality_summary_data"])
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
