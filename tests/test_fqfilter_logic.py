import os
import sys
import unittest

from mfsflow.scripts.fqfilter import extract_seq, hamming_distance, parse_definition


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


if __name__ == "__main__":
    unittest.main()
