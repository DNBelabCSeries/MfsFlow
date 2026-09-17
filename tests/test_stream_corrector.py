import os
import sys
import unittest
import subprocess
import tempfile
from mfsflow.scripts import stream_corrector

from mfsflow.scripts.stream_corrector import get_or_apply_correction


class FakeRead:
    def __init__(self, flag, seq, qual, tags):
        self.flag = flag
        self.query_sequence = seq
        self.query_qualities = qual
        self.tags = dict(tags)

    def has_tag(self, tag):
        return tag in self.tags

    def get_tag(self, tag):
        return self.tags[tag]

    def set_tag(self, tag, value):
        if value is None:
            self.tags.pop(tag, None)
        else:
            self.tags[tag] = value


class StreamCorrectorTests(unittest.TestCase):
    @unittest.skipIf(stream_corrector.pysam is None, "pysam not installed")
    def test_uncompressed_pipe_preserves_read_order_and_tags(self):
        pysam = stream_corrector.pysam
        with tempfile.TemporaryDirectory() as root:
            source = os.path.join(root, "input.bam")
            output = os.path.join(root, "output.bam")
            idmap = os.path.join(root, "ids.tsv")
            with open(idmap, "w") as handle:
                handle.write("wellID\tumi_barcode\tinternal_barcode\nP1A1\tACGT\tTGCA\n")
            with pysam.AlignmentFile(source, "wb", header={"HD": {"VN": "1.6"}}) as bam:
                for flag in (77, 141):
                    read = pysam.AlignedSegment()
                    read.query_name = "pair1"
                    read.flag = flag
                    read.query_sequence = "AACCGGTT"
                    read.query_qualities = [30] * 8
                    read.set_tags([("CR", "ACGT"), ("CC", "ACGT"), ("CB", "P1A1"), ("UR", "AAAA")])
                    bam.write(read)
            with open(output, "wb") as handle:
                subprocess.run(
                    [sys.executable, "-m", "mfsflow.scripts.stream_corrector",
                     "--binning", os.devnull, "--idmap", idmap, "--type", "umi", source],
                    stdout=handle, stderr=subprocess.PIPE, check=True,
                )
            with pysam.AlignmentFile(source, "rb", check_sq=False) as before, \
                 pysam.AlignmentFile(output, "rb", check_sq=False) as after:
                self.assertEqual([r.to_string() for r in before], [r.to_string() for r in after])

    def test_pre_corrected_read_is_not_sequence_adjusted_again(self):
        read = FakeRead(
            flag=77,
            seq="AAACCCCC",
            qual=[30, 31, 32, 33, 34, 35, 36, 37],
            tags={"CR": "RAWBC", "CC": "UMIBC", "CB": "P1A1", "UR": "UMI"},
        )

        correction = get_or_apply_correction(read, {}, {"UMIBC": "P1A1"}, set())

        self.assertFalse(correction.is_internal)
        self.assertEqual(read.query_sequence, "AAACCCCC")
        self.assertEqual(read.query_qualities, [30, 31, 32, 33, 34, 35, 36, 37])

    def test_missing_pre_corrected_tag_falls_back_to_correction(self):
        read = FakeRead(
            flag=77,
            seq="AAACCCCC",
            qual=[30, 31, 32, 33, 34, 35, 36, 37],
            tags={"CR": "RAWBC", "CB": "STALE"},
        )

        correction = get_or_apply_correction(read, {}, {"RAWBC": "P1A1"}, set())

        self.assertEqual(correction.well_id, "P1A1")
        self.assertEqual(read.tags["CC"], "RAWBC")
        self.assertEqual(read.tags["CB"], "P1A1")


if __name__ == "__main__":
    unittest.main()
