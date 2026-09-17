import sys
import unittest
from unittest import mock

try:
    import yaml  # noqa: F401
except ModuleNotFoundError:
    sys.modules["yaml"] = mock.Mock()

from mfsflow.pipeline_modules import _parallel_pair_command, _split_job_budget


class PipelineModuleSchedulingTests(unittest.TestCase):
    def test_gnu_split_budget_accounts_for_stream_count(self):
        self.assertEqual(_split_job_budget(20, ["R1.fq"]), 20)
        self.assertEqual(_split_job_budget(20, ["R1.fq.gz"]), 10)
        self.assertEqual(
            _split_job_budget(20, ["R1.fq", "R2.fq"], paired=True),
            10,
        )
        self.assertEqual(
            _split_job_budget(20, ["R1.fq.gz", "R2.fq.gz"], paired=True),
            10,
        )

    def test_seqkit_budget_is_bounded_by_assigned_threads(self):
        self.assertEqual(
            _split_job_budget(
                20,
                ["R1.fq.gz", "R2.fq.gz"],
                paired=True,
                use_seqkit=True,
            ),
            20,
        )

    def test_parallel_pair_command_waits_for_both_mates(self):
        command = _parallel_pair_command("split-r1", "split-r2")
        self.assertIn("(split-r1) &", command)
        self.assertIn("(split-r2) &", command)
        self.assertIn("wait $mfsflow_r1", command)
        self.assertIn("wait $mfsflow_r2", command)
        self.assertIn("mfsflow_s1", command)
        self.assertIn("mfsflow_s2", command)


if __name__ == "__main__":
    unittest.main()
