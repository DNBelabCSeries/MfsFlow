import os
import sys
import tempfile
import unittest

from mfsflow.scripts.run_featurecounts import _best_intron_assignment, load_saf_interval_index, parse_gtf_and_create_saf


class SafBuilderTests(unittest.TestCase):
    def test_introns_are_clipped_away_from_overlapping_gene_bodies(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = os.path.join(tmpdir, "genes.gtf")
            with open(gtf, "w") as f:
                f.write('chr1\tT\texon\t100\t120\t.\t+\t.\tgene_id "geneA"; gene_name "A";\n')
                f.write('chr1\tT\texon\t200\t220\t.\t+\t.\tgene_id "geneA"; gene_name "A";\n')
                f.write('chr1\tT\texon\t150\t160\t.\t+\t.\tgene_id "geneB"; gene_name "B";\n')
                f.write('chr1\tT\texon\t400\t420\t.\t+\t.\tgene_id "geneB"; gene_name "B";\n')

            saf, gene_map = parse_gtf_and_create_saf(gtf, os.path.join(tmpdir, "sample"))
            self.assertEqual(gene_map["geneA"], "A")
            self.assertTrue(os.path.exists(os.path.join(tmpdir, "sample.exon.saf")))
            self.assertTrue(os.path.exists(os.path.join(tmpdir, "sample.intron.saf")))

            with open(saf) as f:
                lines = [line.strip().split("\t") for line in f if not line.startswith("GeneID")]

            introns = [row for row in lines if row[0] == "geneA__INTRON__"]
            self.assertEqual(
                introns,
                [
                    ["geneA__INTRON__", "chr1", "121", "149", "+"],
                ],
            )

            self.assertNotIn(["geneA__INTRON__", "chr1", "150", "160", "+"], introns)
            self.assertNotIn(["geneA__INTRON__", "chr1", "161", "199", "+"], introns)

    def test_introns_respect_r_length_limits(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = os.path.join(tmpdir, "genes.gtf")
            with open(gtf, "w") as f:
                f.write('chr1\tT\texon\t100\t120\t.\t+\t.\tgene_id "short"; gene_name "short";\n')
                f.write('chr1\tT\texon\t126\t140\t.\t+\t.\tgene_id "short"; gene_name "short";\n')
                f.write('chr1\tT\texon\t1000\t1020\t.\t+\t.\tgene_id "long"; gene_name "long";\n')
                f.write('chr1\tT\texon\t101050\t101070\t.\t+\t.\tgene_id "long"; gene_name "long";\n')

            saf, _gene_map = parse_gtf_and_create_saf(gtf, os.path.join(tmpdir, "sample"))
            with open(saf) as f:
                content = f.read()

            self.assertNotIn("short__INTRON__", content)
            self.assertNotIn("long__INTRON__", content)

    def test_hybrid_intron_assignment_uses_intron_saf(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            saf = os.path.join(tmpdir, "introns.saf")
            with open(saf, "w") as f:
                f.write("GeneID\tChr\tStart\tEnd\tStrand\n")
                f.write("geneA__INTRON__\tchr1\t121\t149\t+\n")
                f.write("geneB__INTRON__\tchr1\t200\t230\t-\n")

            index, starts = load_saf_interval_index(saf)

            gene, category = _best_intron_assignment("chr1", [(124, 130)], 0, False, index, starts)
            self.assertEqual((gene, category), ("geneA", "Intron"))

            gene, category = _best_intron_assignment("chr1", [(124, 130)], 1, True, index, starts)
            self.assertEqual((gene, category), (None, "Intergenic"))


if __name__ == "__main__":
    unittest.main()
