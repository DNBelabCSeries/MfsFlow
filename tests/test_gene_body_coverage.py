import os
import tempfile
import unittest
from unittest import mock

from mfsflow.scripts.run_featurecounts import (
    _project_blocks_to_gene_body,
    _pysam_blocks_1based_half_open,
    coverage_sampling_fraction,
    load_gene_models,
    parse_featurecounts_mapped_reads,
    estimate_primary_mapped_reads,
    primary_mapped_for_coverage,
    should_sample_for_coverage,
)


class GeneBodyCoverageTests(unittest.TestCase):
    def test_transcript_body_projection_keeps_strand_orientation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = os.path.join(tmpdir, "genes.gtf")
            with open(gtf, "w") as f:
                f.write('chr1\tT\texon\t1\t100\t.\t+\t.\tgene_id "plus"; transcript_id "plus_t1"; gene_name "plus";\n')
                f.write('chr1\tT\texon\t1\t100\t.\t-\t.\tgene_id "minus"; transcript_id "minus_t1"; gene_name "minus";\n')

            models = load_gene_models(gtf)

            plus = models["plus"]
            minus = models["minus"]

            plus_overlap, plus_bins = _project_blocks_to_gene_body(plus, [(1, 11)])
            minus_overlap, minus_bins = _project_blocks_to_gene_body(minus, [(91, 101)])

            self.assertEqual(10, plus_overlap)
            self.assertEqual(10, minus_overlap)
            self.assertEqual([1.0] * 10, plus_bins[:10])
            self.assertEqual([1.0] * 10, minus_bins[:10])
            self.assertEqual(0.0, plus_bins[10])
            self.assertEqual(0.0, minus_bins[10])

    def test_longest_transcript_is_selected_per_gene(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = os.path.join(tmpdir, "genes.gtf")
            with open(gtf, "w") as f:
                f.write('chr1\tT\texon\t1\t100\t.\t+\t.\tgene_id "gene1"; transcript_id "short";\n')
                f.write('chr1\tT\texon\t1\t150\t.\t+\t.\tgene_id "gene1"; transcript_id "long";\n')

            models = load_gene_models(gtf)
            self.assertEqual("long", models["gene1"]["transcript_id"])
            self.assertEqual(150, models["gene1"]["length"])

    def test_pysam_blocks_are_converted_to_1based_half_open(self):
        class DummyRead:
            def get_blocks(self):
                return [(0, 100), (150, 175)]

        self.assertEqual(
            [(1, 101), (151, 176)],
            _pysam_blocks_1based_half_open(DummyRead()),
        )

    def test_coverage_sampling_fraction_targets_max_reads(self):
        self.assertEqual(1.0, coverage_sampling_fraction(100, 500))
        self.assertEqual(0.5, coverage_sampling_fraction(1000, 500))
        self.assertEqual(1.0, coverage_sampling_fraction(0, 500))
        self.assertEqual(1.0, coverage_sampling_fraction(1000, 0))

    def test_coverage_sampling_is_stable_by_read_name(self):
        first = should_sample_for_coverage("readA", 0.25, seed=42, source_label="UMI")
        second = should_sample_for_coverage("readA", 0.25, seed=42, source_label="UMI")
        self.assertEqual(first, second)
        self.assertTrue(should_sample_for_coverage("readA", 1.0, seed=42))
        self.assertFalse(should_sample_for_coverage("readA", 0.0, seed=42))


class FeatureCountsSummaryParserTests(unittest.TestCase):
    """Locks the pair-level behavior of parse_featurecounts_mapped_reads.

    It derives the denominator from the featureCounts .summary file, which
    includes supplementary alignments (FLAG 0x800) that --primary does not
    filter. The coverage loop skips supplementary, so this remains an
    approximation, but PE pairs are counted once rather than once per mate.
    """

    def _write_summary(self, tmpdir, bam_name, rows):
        bam_path = os.path.join(tmpdir, bam_name)
        summary_path = bam_path[:-4] + ".counts.txt.summary"
        with open(summary_path, "w") as fh:
            fh.write("Status\tSample\n")
            for key, val in rows.items():
                fh.write(f"{key}\t{val}\n")
        return bam_path

    def test_pe_uses_fragment_count_without_doubling(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bam = self._write_summary(tmpdir, "s1.bam", {
                "Assigned": 1000,
                "Unassigned_Unmapped": 200,
                "Unassigned_NoFeatures": 50,
            })
            # PE summary rows are already pair/fragment counts.
            self.assertEqual(1050, parse_featurecounts_mapped_reads(bam, "PE"))

    def test_se_does_not_double(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bam = self._write_summary(tmpdir, "s2.bam", {
                "Assigned": 500,
                "Unassigned_Unmapped": 100,
            })
            self.assertEqual(500, parse_featurecounts_mapped_reads(bam, "SE"))

    def test_supplementary_inflation_is_intentional(self):
        """featureCounts --primary does not filter supplementary (0x800).

        Supplementary alignments land in Assigned or Unassigned_NoFeatures,
        so they are included in the summary total. This test locks that
        behavior: the returned count includes supplementary, which is the
        known source of denominator inflation. Do NOT change this to subtract
        supplementary — the summary cannot distinguish them.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            bam = self._write_summary(tmpdir, "s3.bam", {
                "Assigned": 800,          # includes ~200 supplementary
                "Unassigned_NoFeatures": 100,
                "Unassigned_Unmapped": 50,
            })
            result = parse_featurecounts_mapped_reads(bam, "SE")
            # sum - unmapped = 800 + 100 + 50 - 50 = 900
            # This is 900, NOT the true primary-only count of ~700.
            self.assertEqual(900, result)

    def test_missing_summary_returns_none(self):
        self.assertIsNone(parse_featurecounts_mapped_reads("/nonexistent.bam", "PE"))

    def test_malformed_lines_are_skipped(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bam_path = os.path.join(tmpdir, "s4.bam")
            summary_path = bam_path[:-4] + ".counts.txt.summary"
            with open(summary_path, "w") as fh:
                fh.write("Status\tSample\n")
                fh.write("Assigned\t1000\n")
                fh.write("Garbage\tNaN\n")
                fh.write("Unassigned_Unmapped\t100\n")
            # PE summary rows are already pair/fragment counts.
            self.assertEqual(1000, parse_featurecounts_mapped_reads(bam_path, "PE"))

    def test_fallback_to_estimate_when_summary_missing(self):
        import mfsflow.scripts.run_featurecounts as rf
        orig = rf.estimate_primary_mapped_reads
        rf.estimate_primary_mapped_reads = lambda b, s, layout: 4242
        try:
            self.assertEqual(4242, primary_mapped_for_coverage("/nonexistent.bam", "samtools", "PE"))
        finally:
            rf.estimate_primary_mapped_reads = orig

    def test_fallback_estimate_counts_primary_r1_for_pe(self):
        import mfsflow.scripts.run_featurecounts as rf

        with mock.patch.object(rf.subprocess, "check_output", return_value="101\n"):
            self.assertEqual(101, estimate_primary_mapped_reads("sample.bam", "samtools", "PE"))
            self.assertEqual(101, estimate_primary_mapped_reads("sample.bam", "samtools", "SE"))
            pe_cmd = rf.subprocess.check_output.call_args_list[0].args[0]
            se_cmd = rf.subprocess.check_output.call_args_list[1].args[0]
            self.assertEqual(["-f", "64"], pe_cmd[5:7])
            self.assertNotIn("-f", se_cmd)


if __name__ == "__main__":
    unittest.main()
