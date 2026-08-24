import io
import os
import sys
import unittest

from mfsflow.scripts.fqfilter import (
    extract_seq,
    fastq_iter,
    hamming_distance,
    iter_synchronized_fastq,
    parse_definition,
)


class FqfilterLogicTests(unittest.TestCase):
    def test_parse_definition_supports_multiple_ranges(self):
        self.assertEqual(
            parse_definition("BC(1-4,9-10);UMI(5-8);cDNA(11-20)"),
            {
                "BC": [(0, 4), (8, 10)],
                "UMI": [(4, 8)],
                "cDNA": [(10, 20)],
            },
        )

    def test_smartseq3_no_pattern_keeps_full_read_as_cdna(self):
        definition = {
            "UMI": [(0, 10)],
            "cDNA": [(10, 20)],
        }
        bc, bc_q, umi, umi_q, cdna, cdna_q = extract_seq(
            b"ACGTACGTACGTACGTACGT",
            b"IIIIIIIIIIIIIIIIIIII",
            definition,
            ss3_no_pattern=True,
        )
        self.assertEqual(bc, b"")
        self.assertEqual(bc_q, b"")
        self.assertEqual(umi, b"")
        self.assertEqual(umi_q, b"")
        self.assertEqual(cdna, b"ACGTACGTACGTACGTACGT")
        self.assertEqual(cdna_q, b"IIIIIIIIIIIIIIIIIIII")

    def test_smartseq3_pattern_hamming_limit(self):
        self.assertEqual(hamming_distance(b"ATTGCGCAATG", b"ATTGCGCAATA", limit=1), 1)
        self.assertEqual(hamming_distance(b"ATTGCGCAATG", b"TTTGCGCAATA", limit=1), 2)

    def test_fastq_iter_rejects_incomplete_record(self):
        handle = io.BytesIO(b"@read1\nACGT\n+\n")
        with self.assertRaisesRegex(ValueError, "Incomplete FASTQ record 1"):
            list(fastq_iter(handle, "R1.fastq"))

    def test_fastq_iter_rejects_sequence_quality_length_mismatch(self):
        handle = io.BytesIO(b"@read1\nACGT\n+\nIII\n")
        with self.assertRaisesRegex(ValueError, "SEQ/QUAL length mismatch"):
            list(fastq_iter(handle, "R1.fastq"))

    def test_synchronized_fastq_rejects_unequal_mate_counts(self):
        r1 = io.BytesIO(b"@read1/1\nACGT\n+\nIIII\n@read2/1\nACGT\n+\nIIII\n")
        r2 = io.BytesIO(b"@read1/2\nACGT\n+\nIIII\n")
        with self.assertRaisesRegex(ValueError, "FASTQ mate count mismatch at record 2"):
            list(iter_synchronized_fastq([r1, r2], ["R1.fastq", "R2.fastq"]))

    def test_synchronized_fastq_rejects_mismatched_mate_ids(self):
        r1 = io.BytesIO(b"@read1/1\nACGT\n+\nIIII\n")
        r2 = io.BytesIO(b"@other/2\nACGT\n+\nIIII\n")
        with self.assertRaisesRegex(ValueError, "FASTQ mate ID mismatch"):
            list(iter_synchronized_fastq([r1, r2], ["R1.fastq", "R2.fastq"]))


if __name__ == "__main__":
    unittest.main()
