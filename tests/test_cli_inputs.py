import os
import sys
import gzip
import tempfile
import unittest

from mfsflow.config.builder import (
    configure_reference,
    discover_fastq_pairs,
    load_samplesheet,
    resolve_samplesheet_barcodes,
    resolve_samplesheet_fastq_groups,
)


class CliInputTests(unittest.TestCase):
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
