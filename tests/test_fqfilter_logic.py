import os
import sys
import unittest

from mfsflow.scripts.fqfilter import extract_seq, hamming_distance


class FqfilterLogicTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
