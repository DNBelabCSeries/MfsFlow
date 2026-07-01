import os
import sys
import tempfile
import unittest

from mfsflow.scripts.barcode_corrector import correct_read_barcode, load_bc_map, load_id_map


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


class BarcodeCorrectorTests(unittest.TestCase):
    def test_load_maps(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            binmap = os.path.join(tmpdir, "binning.txt")
            idmap = os.path.join(tmpdir, "expect_id_barcode.tsv")

            with open(binmap, "w") as f:
                f.write("falseBC,hamming,trueBC,n\n")
                f.write("aaaa,1,CCCC,7\n")

            with open(idmap, "w") as f:
                f.write("wellID\tumi_barcodes\tinternal_barcodes\n")
                f.write("P1A1\tCCCC,GGGG\tTTTT\n")

            self.assertEqual(load_bc_map(binmap), {"AAAA": "CCCC"})
            id_map, internal = load_id_map(idmap)
            self.assertEqual(id_map["CCCC"], "P1A1")
            self.assertEqual(id_map["GGGG"], "P1A1")
            self.assertEqual(id_map["TTTT"], "P1A1")
            self.assertEqual(internal, {"TTTT"})

    def test_umi_read_is_trimmed_and_tagged(self):
        read = FakeRead(
            flag=77,
            seq="AAACCCCC",
            qual=[30, 31, 32, 33, 34, 35, 36, 37],
            tags={"CR": "aaaa", "UR": "UMI", "UY": "III"},
        )

        correction = correct_read_barcode(
            read,
            bc_map={"AAAA": "CCCC"},
            id_map={"CCCC": "P1A1"},
            internal_bcs={"TTTT"},
        )

        self.assertFalse(correction.is_internal)
        self.assertEqual(read.query_sequence, "CCCCC")
        self.assertEqual(read.query_qualities, [33, 34, 35, 36, 37])
        self.assertEqual(read.tags["CR"], "AAAA")
        self.assertEqual(read.tags["CC"], "CCCC")
        self.assertEqual(read.tags["CB"], "P1A1")
        self.assertEqual(read.tags["UR"], "UMI")

    def test_internal_read_prepends_umi_and_clears_umi_tags(self):
        read = FakeRead(
            flag=77,
            seq="GGGG",
            qual=[20, 21, 22, 23],
            tags={"CR": "tttt", "UR": "ACG", "UY": "ABC"},
        )

        correction = correct_read_barcode(
            read,
            bc_map={},
            id_map={"TTTT": "P1A1"},
            internal_bcs={"TTTT"},
        )

        self.assertTrue(correction.is_internal)
        self.assertEqual(read.query_sequence, "ACGGGGG")
        self.assertEqual(read.query_qualities, [32, 33, 34, 20, 21, 22, 23])
        self.assertNotIn("UR", read.tags)
        self.assertNotIn("UY", read.tags)
        self.assertEqual(read.tags["CB"], "P1A1")

    def test_missing_barcode_returns_none(self):
        read = FakeRead(flag=77, seq="AAAA", qual=[30, 30, 30, 30], tags={})
        self.assertIsNone(correct_read_barcode(read, {}, {}, set()))

    def test_unmatched_barcode_clears_existing_cc_cb_tags(self):
        """Covers the tag-clearing branch: raw_bc not in id_map → well_id is None.

        If the read previously had CC/CB tags (e.g. from a prior correction
        pass), they must be removed, not left stale.
        """
        read = FakeRead(
            flag=77,
            seq="CCCCC",
            qual=[30, 31, 32, 33, 34],
            tags={"CR": "xxxx", "UR": "UMI", "UY": "III", "CC": "OLD", "CB": "OLD_WELL"},
        )

        correction = correct_read_barcode(
            read,
            bc_map={"XXXX": "YYYY"},  # XXXX corrects to YYYY, but YYYY not in id_map
            id_map={},
            internal_bcs=set(),
        )

        self.assertFalse(correction.is_internal)
        self.assertIsNone(correction.well_id)
        self.assertEqual(read.tags["CR"], "XXXX")
        self.assertNotIn("CC", read.tags)
        self.assertNotIn("CB", read.tags)

    def test_unmatched_barcode_with_no_prior_tags(self):
        """Clearing tags that don't exist should not raise."""
        read = FakeRead(
            flag=77,
            seq="CCCCC",
            qual=[30, 31, 32, 33, 34],
            tags={"CR": "xxxx", "UR": "UMI", "UY": "III"},
        )

        correction = correct_read_barcode(
            read,
            bc_map={"XXXX": "YYYY"},
            id_map={},
            internal_bcs=set(),
        )

        self.assertIsNone(correction.well_id)
        self.assertNotIn("CC", read.tags)
        self.assertNotIn("CB", read.tags)


if __name__ == "__main__":
    unittest.main()
