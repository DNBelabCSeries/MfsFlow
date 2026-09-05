import os
import sys
import unittest

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
