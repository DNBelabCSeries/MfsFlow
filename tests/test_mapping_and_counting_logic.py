import os
import sys
import gzip
import tempfile
import unittest

from mfsflow.scripts.mapping_analysis import build_star_misc_base, setup_gtf
from mfsflow.scripts.run_featurecounts import build_featurecounts_cmd, normalize_read_category, resolve_counting_strand_modes, should_count_read


class MappingAndCountingLogicTests(unittest.TestCase):
    def test_build_star_misc_base_uses_target_specific_overhang(self):
        misc = build_star_misc_base(
            "/ref/star",
            8,
            "PE",
            False,
            "/tmp/genes.gtf",
            101,
        )
        self.assertIn("--sjdbOverhang 100", misc)

    def test_build_star_misc_base_skips_overhang_when_index_has_sjdb(self):
        misc = build_star_misc_base(
            "/ref/star",
            8,
            "PE",
            True,
            "/tmp/genes.gtf",
            0,
        )
        self.assertNotIn("--sjdbOverhang", misc)
        self.assertNotIn("--sjdbGTFfile", misc)

    def test_setup_gtf_decompresses_gzipped_gtf(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gz_gtf = os.path.join(tmpdir, "genes.gtf.gz")
            with gzip.open(gz_gtf, "wt") as handle:
                handle.write('chr1\tT\texon\t1\t10\t.\t+\t.\tgene_id "g1";\n')

            final_gtf, extra = setup_gtf(
                {"reference": {"GTF_file": gz_gtf, "STAR_index": tmpdir}},
                "Sample01",
                tmpdir,
                "samtools",
            )

            self.assertEqual(extra, "")
            self.assertTrue(final_gtf.endswith("Sample01.final_annot.gtf"))
            with open(final_gtf) as handle:
                self.assertIn('gene_id "g1"', handle.read())

    def test_featurecounts_pe_layout_enables_pair_options(self):
        cmd = build_featurecounts_cmd(
            "featureCounts",
            "in.bam",
            "features.saf",
            "out.txt",
            4,
            1,
            "PE",
        )
        self.assertIn("-p", cmd)
        self.assertIn("-C", cmd)

    def test_featurecounts_se_layout_does_not_enable_pair_options(self):
        cmd = build_featurecounts_cmd(
            "featureCounts",
            "in.bam",
            "features.saf",
            "out.txt",
            4,
            1,
            "SE",
        )
        self.assertNotIn("-p", cmd)
        self.assertNotIn("-C", cmd)

    def test_resolve_counting_strand_modes_preserves_historical_internal_default(self):
        umi_mode, internal_mode = resolve_counting_strand_modes({"strand": 2})
        self.assertEqual((umi_mode, internal_mode), (2, 0))

    def test_resolve_counting_strand_modes_supports_explicit_internal_override(self):
        umi_mode, internal_mode = resolve_counting_strand_modes({"strand": 1, "internal_strand": 0})
        self.assertEqual((umi_mode, internal_mode), (1, 0))

    def test_should_count_read_filters_secondary_and_supplementary(self):
        self.assertTrue(should_count_read("q1\t99\tchr1\t1\t255\t50M\t=\t1\t0\tACGT\tFFFF"))
        self.assertFalse(should_count_read("q2\t355\tchr1\t1\t255\t50M\t=\t1\t0\tACGT\tFFFF"))
        self.assertFalse(should_count_read("q3\t2147\tchr1\t1\t255\t50M\t=\t1\t0\tACGT\tFFFF"))

    def test_normalize_read_category_collapses_unassigned_reasons(self):
        self.assertEqual(normalize_read_category("Exon"), "Exon")
        self.assertEqual(normalize_read_category("MappingQuality"), "Other_Unassigned")
        self.assertEqual(normalize_read_category("FragmentLength"), "Other_Unassigned")


if __name__ == "__main__":
    unittest.main()
