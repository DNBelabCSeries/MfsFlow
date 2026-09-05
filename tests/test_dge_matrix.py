import gzip
import os
import tempfile
import unittest

from mfsflow.scripts.dge_analysis import write_sparse_matrix


class DgeMatrixTests(unittest.TestCase):
    def test_sparse_matrix_streaming_preserves_sorted_mex_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            write_sparse_matrix(
                {
                    "BC2": {"G2": 2, "G1": 1},
                    "BC1": {"G2": 3, "G3": 0},
                },
                ["G1", "G2", "G3"],
                {"G1": "Gene 1", "G2": "Gene 2", "G3": "Gene 3"},
                ["BC1", "BC2"],
                tmpdir,
                "sample.exon.umi",
            )

            matrix_path = os.path.join(tmpdir, "expression", "sample.exon.umi", "matrix.mtx.gz")
            with gzip.open(matrix_path, "rt", encoding="utf-8") as handle:
                lines = [line.strip() for line in handle if line.strip() and not line.startswith("%")]

            self.assertEqual(lines[0], "3 2 3")
            self.assertEqual(lines[1:], ["2 1 3", "1 2 1", "2 2 2"])


if __name__ == "__main__":
    unittest.main()
